#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    CPU pipeline for Nanopore ITS amplicon analysis.

    Processes ONE trimmed FASTQ file (post-Porechop) through:
        1. Bayesian hyperparameter optimization via Optuna
           (VSEARCH clustering + MUSCLE consensus + single Racon round)
        2. Selection of the best-performing trial
        3. VSEARCH UCHIME de-novo chimera removal on the best consensus

    This script calls optimize_pipeline.py internally to run the
    Optuna study and inherits its best result.

Usage:
    python pipeline/run_cpu_pipeline.py \
        --input  trimmed/B33_trimmed.fastq \
        --output results/B33 \
        --blast-db  db/unite2024ITS \
        --threads 4 \
        --trials 12

    # To process all barcodes in a loop:
    for fq in trimmed/*.fastq; do
        bc=$(basename "$fq" _trimmed.fastq)
        python pipeline/run_cpu_pipeline.py \
            --input "$fq" \
            --output "results/$bc" \
            --blast-db db/unite2024ITS \
            --threads 4 --trials 12
    done

Requirements:
    vsearch, minimap2, racon, muscle, blastn, optuna, biopython
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


def run_cmd(cmd, step_name, stdout_file=None):
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


def run_chimera_removal(best_consensus, chimera_dir, vsearch_bin="vsearch"):
    """Dereplication + VSEARCH UCHIME de-novo chimera removal."""
    chimera_dir.mkdir(parents=True, exist_ok=True)
    derep   = chimera_dir / "derep.fasta"
    nonchim = chimera_dir / "nonchimeras.fasta"
    chim    = chimera_dir / "chimeras.fasta"

    run_cmd(
        [vsearch_bin, "--derep_fulllength", str(best_consensus),
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


def main():
    parser = argparse.ArgumentParser(
        description="CPU Nanopore ITS pipeline: Optuna(VSEARCH+Racon) → UCHIME"
    )
    parser.add_argument("--input", required=True,
                        help="Trimmed input FASTQ (post-Porechop)")
    parser.add_argument("--output", required=True,
                        help="Output directory for this barcode")
    parser.add_argument("--blast-db", required=True,
                        help="BLAST database prefix (UNITE-based)")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--trials", type=int, default=12,
                        help="Number of Optuna trials (default: 12)")
    parser.add_argument("--vsearch-bin",  default="vsearch")
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--racon-bin",    default="racon")
    parser.add_argument("--muscle-bin",   default="muscle")
    parser.add_argument("--blastn-bin",   default="blastn")
    args = parser.parse_args()

    input_fastq = Path(args.input)
    output_dir  = Path(args.output)

    if not input_fastq.exists():
        logger.error("Input FASTQ not found: %s", input_fastq)
        sys.exit(1)

    vsearch_dir = output_dir / "VSEARCH"
    vsearch_dir.mkdir(parents=True, exist_ok=True)

    # ── Optuna optimization ──────────────────────────────────
    logger.info("Starting Optuna optimization (%d trials) …", args.trials)

    # optimize_pipeline.py lives in the same directory as this script
    script_dir = Path(__file__).parent
    optimize_script = script_dir / "optimize_pipeline.py"

    run_cmd(
        [sys.executable, str(optimize_script),
         "--input",       str(input_fastq),
         "--output",      str(vsearch_dir),
         "--blast-db",    args.blast_db,
         "--trials",      str(args.trials),
         "--threads",     str(args.threads),
         "--vsearch-bin", args.vsearch_bin,
         "--minimap2-bin", args.minimap2_bin,
         "--racon-bin",   args.racon_bin,
         "--muscle-bin",  args.muscle_bin,
         "--blastn-bin",  args.blastn_bin],
        "Optuna-optimization",
    )

    best_consensus = vsearch_dir / "polished_best.fasta"
    if not best_consensus.exists():
        logger.error("optimize_pipeline.py did not produce polished_best.fasta")
        sys.exit(1)

    # ── Chimera removal ──────────────────────────────────────
    chimera_dir = output_dir / "UCHIME"
    final_fasta = run_chimera_removal(best_consensus, chimera_dir, args.vsearch_bin)

    logger.info("Final non-chimeric consensus: %s", final_fasta)
    logger.info("CPU pipeline completed successfully.")


if __name__ == "__main__":
    main()
