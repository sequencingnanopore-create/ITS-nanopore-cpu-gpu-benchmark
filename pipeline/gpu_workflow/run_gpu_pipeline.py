#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    GPU pipeline for Nanopore ITS amplicon analysis.

    Processes ONE trimmed FASTQ file (post-Porechop) through:
        1. Amplicon Sorter  (primary draft consensus)
        2. Rescue pathway   (minimap2 + single Racon if step 1 fails)
        3. Iterative Racon  (3 rounds)
        4. Medaka           (neural-network polishing)
        5. VSEARCH UCHIME   (de-novo chimera removal)

Usage:
    python pipeline/run_gpu_pipeline.py \
        --input  trimmed/B33_trimmed.fastq \
        --output results/B33 \
        --threads 8 \
        --medaka-model r941_min_sup_g507

    # To process all barcodes in a loop:
    for fq in trimmed/*.fastq; do
        bc=$(basename "$fq" _trimmed.fastq)
        python pipeline/run_gpu_pipeline.py \
            --input "$fq" \
            --output "results/$bc" \
            --threads 8 \
            --medaka-model r941_min_sup_g507
    done

Requirements:
    amplicon_sorter, minimap2, racon, medaka_consensus, vsearch
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

def run_cmd(cmd, step_name, stdout_file=None):
    """Run a shell command with structured logging.
    If stdout_file is given, redirect stdout there (required by Racon).
    """
    logger.info("[%s] %s", step_name, " ".join(str(c) for c in cmd))
    try:
        if stdout_file is not None:
            with open(stdout_file, "w") as fh:
                subprocess.run(cmd, stdout=fh, check=True)
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        logger.error("[%s] Failed (exit code %d).", step_name, exc.returncode)
        raise


# ─────────────────────────────────────────────────────────────
# Step 1 — Amplicon Sorter
# ─────────────────────────────────────────────────────────────

def run_amplicon_sorter(input_fastq, output_dir, similarity_cutoff):
    """
    Run Amplicon Sorter in single-sample mode.
    Returns path to draft consensus FASTA, or None if the tool fails
    or produces no output (rescue pathway is triggered in that case).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3", "amplicon_sorter/amplicon_sorter.py",
        "-i", str(input_fastq),
        "-o", str(output_dir),
        "-sc", str(similarity_cutoff),
    ]

    try:
        run_cmd(cmd, "AmpliconSorter")
    except subprocess.CalledProcessError:
        logger.warning("[AmpliconSorter] Primary step failed → rescue pathway.")
        return None

    # Amplicon Sorter names its consensus *_trimmed_consensussequences.fasta
    candidates = list(output_dir.glob("*_trimmed_consensussequences.fasta"))
    if candidates:
        return candidates[0]

    plain = output_dir / "consensus.fasta"
    if plain.exists():
        return plain

    logger.warning("[AmpliconSorter] No consensus found in output → rescue pathway.")
    return None


# ─────────────────────────────────────────────────────────────
# Rescue pathway
# ─────────────────────────────────────────────────────────────

def rescue_pathway(reads_fastq, rescue_dir, threads):
    """
    Fallback when Amplicon Sorter fails.
    Aligns reads against themselves (overlap mode) and runs one Racon
    round to produce a provisional draft.
    """
    logger.info("[Rescue] Building draft via minimap2 (ava-ont) + Racon …")
    rescue_dir.mkdir(parents=True, exist_ok=True)

    paf = rescue_dir / "rescue_overlaps.paf"
    run_cmd(
        ["minimap2", "-x", "ava-ont", "-t", str(threads),
         str(reads_fastq), str(reads_fastq)],
        "minimap2-rescue",
        stdout_file=paf,
    )

    draft = rescue_dir / "rescue_draft.fasta"
    # BUG FIX 1: Racon writes its output to stdout, not a positional argument.
    run_cmd(
        ["racon", "-t", str(threads),
         str(reads_fastq), str(paf), str(reads_fastq)],
        "Racon-rescue",
        stdout_file=draft,
    )

    if not draft.exists() or draft.stat().st_size == 0:
        raise RuntimeError("Rescue pathway produced an empty or missing draft.")

    return draft


# ─────────────────────────────────────────────────────────────
# Step 2 — Iterative Racon
# ─────────────────────────────────────────────────────────────

def _minimap2(reference, reads, paf, threads):
    run_cmd(
        ["minimap2", "-x", "map-ont", "-t", str(threads),
         str(reference), str(reads)],
        "minimap2",
        stdout_file=paf,
    )


def run_racon_iterative(reads_fastq, draft_fasta, racon_dir, threads, rounds=3):
    """
    Run `rounds` iterative Racon rounds.
    BUG FIX 2: We loop for 3 rounds (matching the Methods section).
    BUG FIX 1: Each Racon call writes to stdout captured to a file.
    """
    racon_dir.mkdir(parents=True, exist_ok=True)
    current_ref = draft_fasta

    for i in range(1, rounds + 1):
        logger.info("[Racon] Round %d / %d", i, rounds)
        paf = racon_dir / f"overlaps_round{i}.paf"
        polished = racon_dir / f"racon_round{i}.fasta"

        _minimap2(current_ref, reads_fastq, paf, threads)
        run_cmd(
            ["racon", "-t", str(threads),
             str(reads_fastq), str(paf), str(current_ref)],
            f"Racon-round{i}",
            stdout_file=polished,
        )
        current_ref = polished

    return current_ref


# ─────────────────────────────────────────────────────────────
# Step 3 — Medaka
# ─────────────────────────────────────────────────────────────

def run_medaka(input_fastq, draft_fasta, output_dir, model, threads):
    """Run Medaka neural-network polishing (requires GPU)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    run_cmd(
        ["medaka_consensus",
         "-i", str(input_fastq),
         "-d", str(draft_fasta),
         "-o", str(output_dir),
         "-m", model,
         "-t", str(threads)],
        "Medaka",
    )
    polished = output_dir / "consensus.fasta"
    if not polished.exists():
        raise FileNotFoundError("Medaka did not produce consensus.fasta.")
    return polished


# ─────────────────────────────────────────────────────────────
# Step 4 — Chimera removal
# BUG FIX 3: This step was missing from the original script.
# ─────────────────────────────────────────────────────────────

def run_chimera_removal(polished_fasta, chimera_dir, vsearch_bin="vsearch"):
    """
    Dereplication + VSEARCH UCHIME de-novo chimera removal.
    Mirrors the process in chimera_removal.py for GPU barcodes.
    """
    chimera_dir.mkdir(parents=True, exist_ok=True)
    derep = chimera_dir / "derep.fasta"
    nonchim = chimera_dir / "nonchimeras.fasta"
    chim = chimera_dir / "chimeras.fasta"

    run_cmd(
        [vsearch_bin, "--derep_fulllength", str(polished_fasta),
         "--output", str(derep), "--sizeout"],
        "VSEARCH-derep",
    )
    run_cmd(
        [vsearch_bin, "--uchime_denovo", str(derep),
         "--nonchimeras", str(nonchim),
         "--chimeras", str(chim)],
        "VSEARCH-uchime",
    )
    return nonchim


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="GPU Nanopore ITS pipeline: AmpliconSorter → Racon×3 → Medaka → UCHIME"
    )
    parser.add_argument("--input", required=True,
                        help="Trimmed input FASTQ (post-Porechop)")
    parser.add_argument("--output", required=True,
                        help="Output directory for this barcode")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--medaka-model", required=True,
                        help="Medaka model string, e.g. r941_min_sup_g507")
    parser.add_argument("--amplicon-sorter-sc", type=int, default=90,
                        help="Amplicon Sorter similarity cutoff (default: 90)")
    parser.add_argument("--racon-rounds", type=int, default=3,
                        help="Number of Racon rounds (default: 3)")
    parser.add_argument("--vsearch-bin", default="vsearch",
                        help="Path to vsearch binary (default: vsearch)")
    args = parser.parse_args()

    input_fastq = Path(args.input)
    output_dir = Path(args.output)

    if not input_fastq.exists():
        logger.error("Input FASTQ not found: %s", input_fastq)
        sys.exit(1)

    # ── 1. Amplicon Sorter ────────────────────────────────────
    sorter_dir = output_dir / "AMPLICON_SORTER"
    draft = run_amplicon_sorter(input_fastq, sorter_dir, args.amplicon_sorter_sc)

    # ── Rescue ────────────────────────────────────────────────
    # BUG FIX 4: rescue pathway was absent in the original script.
    if draft is None:
        rescue_dir = output_dir / "RESCUE"
        try:
            draft = rescue_pathway(input_fastq, rescue_dir, args.threads)
        except RuntimeError as exc:
            logger.error("Both primary and rescue pathways failed: %s", exc)
            sys.exit(1)

    # ── 2. Iterative Racon ────────────────────────────────────
    racon_dir = output_dir / "RACON"
    racon_final = run_racon_iterative(
        input_fastq, draft, racon_dir,
        threads=args.threads, rounds=args.racon_rounds,
    )

    # ── 3. Medaka ─────────────────────────────────────────────
    medaka_dir = output_dir / "MEDAKA"
    medaka_out = run_medaka(
        input_fastq, racon_final, medaka_dir,
        model=args.medaka_model, threads=args.threads,
    )

    # ── 4. Chimera removal ────────────────────────────────────
    chimera_dir = output_dir / "UCHIME"
    final_fasta = run_chimera_removal(medaka_out, chimera_dir, args.vsearch_bin)

    logger.info("Final non-chimeric consensus: %s", final_fasta)
    logger.info("GPU pipeline completed successfully.")


if __name__ == "__main__":
    main()
