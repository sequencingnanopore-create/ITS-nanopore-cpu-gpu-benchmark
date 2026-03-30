#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Performs the final polishing step of the GPU pipeline.
    Runs 3 iterative Racon rounds on the Amplicon Sorter draft,
    followed by Medaka neural-network consensus polishing.
    This script mirrors the polishing steps used in run_gpu_pipeline.py.

Usage:
    python pipeline/gpu_workflow/medaka_polishing.py \
        --input trimmed/B33_trimmed.fastq \
        --draft results/B33/AMPLICON_SORTER/draft.fasta \
        --output results/B33 \
        --threads 8 \
        --medaka-model r941_min_sup_g507

Requirements:
    minimap2, racon, medaka_consensus (in activated medaka_env)
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
# Binary paths (use "command" if they are in $PATH or conda environment)
MINIMAP2_BIN = "minimap2"
RACON_BIN    = "racon"

# Medaka environment and model
MEDAKA_ENV   = "/usr/local/miniconda/envs/medaka_env"
MEDAKA_MODEL = "r941_min_sup_g507"

# Default parameters
DEFAULT_THREADS = 8
DEFAULT_RACON_ROUNDS = 3
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


def run_racon_iterative(reads_fastq: Path, draft_fasta: Path, racon_dir: Path, threads: int, rounds: int):
    """Run the requested number of iterative Racon polishing rounds."""
    racon_dir.mkdir(parents=True, exist_ok=True)
    current_ref = draft_fasta

    for i in range(1, rounds + 1):
        logger.info("[Racon] Round %d / %d", i, rounds)
        paf = racon_dir / f"overlaps_round{i}.paf"
        polished = racon_dir / f"racon_round{i}.fasta"

        # minimap2 mapping
        run_cmd(
            [MINIMAP2_BIN, "-x", "map-ont", "-t", str(threads),
             str(current_ref), str(reads_fastq)],
            f"minimap2-round{i}",
            stdout_file=paf,
        )

        # Racon polishing
        run_cmd(
            [RACON_BIN, "-t", str(threads),
             str(reads_fastq), str(paf), str(current_ref)],
            f"Racon-round{i}",
            stdout_file=polished,
        )

        current_ref = polished

    return current_ref


def run_medaka(input_fastq: Path, draft_fasta: Path, medaka_dir: Path, model: str, threads: int):
    """Run Medaka neural-network polishing (requires GPU/conda env)."""
    medaka_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "-p", MEDAKA_ENV, "medaka_consensus",
        "-i", str(input_fastq),
        "-d", str(draft_fasta),
        "-o", str(medaka_dir),
        "-m", model,
        "-t", str(threads),
    ]

    run_cmd(cmd, "Medaka")
    polished = medaka_dir / "consensus.fasta"

    if not polished.exists():
        raise FileNotFoundError("Medaka did not produce consensus.fasta.")

    return polished


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 9 — Medaka polishing (Racon×3 + Medaka)"
    )
    parser.add_argument("--input", required=True,
                        help="Trimmed input FASTQ (post-Porechop)")
    parser.add_argument("--draft", required=True,
                        help="Initial draft FASTA (Amplicon Sorter or rescue)")
    parser.add_argument("--output", required=True,
                        help="Output directory for this barcode (e.g. results/B33)")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--medaka-model", default=MEDAKA_MODEL,
                        help="Medaka model string (default: r941_min_sup_g507)")
    parser.add_argument("--racon-rounds", type=int, default=DEFAULT_RACON_ROUNDS,
                        help="Number of Racon rounds (default: 3)")
    args = parser.parse_args()

    input_fastq = Path(args.input)
    draft_fasta = Path(args.draft)
    output_dir  = Path(args.output)

    if not input_fastq.exists():
        logger.error("Input FASTQ not found: %s", input_fastq)
        sys.exit(1)
    if not draft_fasta.exists():
        logger.error("Draft FASTA not found: %s", draft_fasta)
        sys.exit(1)

    # ── 1. Iterative Racon ───────────────────────────────────────────
    racon_dir = output_dir / "RACON"
    racon_final = run_racon_iterative(
        input_fastq, draft_fasta, racon_dir,
        threads=args.threads, rounds=args.racon_rounds
    )

    # ── 2. Medaka polishing ──────────────────────────────────────────
    medaka_dir = output_dir / "MEDAKA"
    medaka_out = run_medaka(
        input_fastq, racon_final, medaka_dir,
        model=args.medaka_model, threads=args.threads
    )

    logger.info("Medaka polishing completed successfully.")
    logger.info("Final polished consensus: %s", medaka_out)


if __name__ == "__main__":
    main()