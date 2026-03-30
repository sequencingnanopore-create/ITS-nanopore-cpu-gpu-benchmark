#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Generates the initial draft consensus for the GPU pipeline.
    Primary method: Amplicon Sorter (single-sample mode).
    Fallback (rescue pathway): minimap2 (ava-ont) + one round of Racon
    when Amplicon Sorter fails or produces no output.
    This script mirrors the draft step used in run_gpu_pipeline.py.

Usage:
    python pipeline/gpu_workflow/draft_generation.py \
        --input trimmed/B33_trimmed.fastq \
        --output results/B33/AMPLICON_SORTER \
        --threads 8 \
        --similarity-cutoff 90

Requirements:
    amplicon_sorter/amplicon_sorter.py, minimap2, racon
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


# ═══════════════════════════════════════════════════
# USER CONFIGURATION — Modify only this block
# ═══════════════════════════════════════════════════
# Path to Amplicon Sorter script (relative to repository root)
AMPLICON_SORTER_SCRIPT = "amplicon_sorter/amplicon_sorter.py"

# Default binary paths (change only if not in $PATH)
MINIMAP2_BIN = "minimap2"
RACON_BIN    = "racon"

# Default parameters
DEFAULT_THREADS         = 8
DEFAULT_SIMILARITY_CUTOFF = 90
# ═══════════════════════════════════════════════════


def run_cmd(cmd, step_name, stdout_file=None):
    """Run a shell command with structured logging."""
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


def run_amplicon_sorter(input_fastq: Path, output_dir: Path, similarity_cutoff: int):
    """
    Run Amplicon Sorter in single-sample mode.
    Returns path to draft consensus FASTA, or None if the tool fails
    (rescue pathway will be triggered).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3", str(AMPLICON_SORTER_SCRIPT),
        "-i", str(input_fastq),
        "-o", str(output_dir),
        "-sc", str(similarity_cutoff),
    ]

    try:
        run_cmd(cmd, "AmpliconSorter")
    except subprocess.CalledProcessError:
        logger.warning("[AmpliconSorter] Primary step failed → rescue pathway.")
        return None

    # Amplicon Sorter output naming
    candidates = list(output_dir.glob("*_trimmed_consensussequences.fasta"))
    if candidates:
        return candidates[0]

    plain = output_dir / "consensus.fasta"
    if plain.exists():
        return plain

    logger.warning("[AmpliconSorter] No consensus found → rescue pathway.")
    return None


def rescue_pathway(reads_fastq: Path, rescue_dir: Path, threads: int):
    """
    Fallback pathway: minimap2 (ava-ont) + single Racon round.
    """
    logger.info("[Rescue] Building draft via minimap2 (ava-ont) + Racon …")
    rescue_dir.mkdir(parents=True, exist_ok=True)

    paf = rescue_dir / "rescue_overlaps.paf"
    run_cmd(
        [MINIMAP2_BIN, "-x", "ava-ont", "-t", str(threads),
         str(reads_fastq), str(reads_fastq)],
        "minimap2-rescue",
        stdout_file=paf,
    )

    draft = rescue_dir / "rescue_draft.fasta"
    run_cmd(
        [RACON_BIN, "-t", str(threads),
         str(reads_fastq), str(paf), str(reads_fastq)],
        "Racon-rescue",
        stdout_file=draft,
    )

    if not draft.exists() or draft.stat().st_size == 0:
        raise RuntimeError("Rescue pathway produced an empty or missing draft.")

    return draft


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 5 — Draft consensus generation (Amplicon Sorter + rescue)"
    )
    parser.add_argument("--input", required=True,
                        help="Trimmed input FASTQ (post-Porechop)")
    parser.add_argument("--output", required=True,
                        help="Output directory for draft (e.g. results/B33/AMPLICON_SORTER)")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--similarity-cutoff", type=int, default=DEFAULT_SIMILARITY_CUTOFF,
                        help="Amplicon Sorter similarity cutoff (default: 90)")
    args = parser.parse_args()

    input_fastq = Path(args.input)
    output_dir  = Path(args.output)

    if not input_fastq.exists():
        logger.error("Input FASTQ not found: %s", input_fastq)
        sys.exit(1)

    # ── Primary: Amplicon Sorter ─────────────────────────────────────
    draft = run_amplicon_sorter(input_fastq, output_dir, args.similarity_cutoff)

    # ── Rescue pathway if needed ─────────────────────────────────────
    if draft is None:
        rescue_dir = output_dir.parent / "RESCUE"
        try:
            draft = rescue_pathway(input_fastq, rescue_dir, args.threads)
            logger.info("[Rescue] Draft successfully generated at %s", draft)
        except RuntimeError as exc:
            logger.error("Both primary and rescue pathways failed: %s", exc)
            sys.exit(1)

    logger.info("Draft generation completed successfully.")
    logger.info("Final draft consensus: %s", draft)


if __name__ == "__main__":
    main()