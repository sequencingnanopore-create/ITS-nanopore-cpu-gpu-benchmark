#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28

Description:
    Calculates substitution (S), insertion (I) and deletion (D) error rates
    by aligning Nanopore reads against reference sequences using edlib.
    Processes trimmed FASTQ files per barcode and generates a summary table
    with mean and standard error for each error type (used in stacked error
    plots and Supplementary Material).

Usage:
    python analysis/error_profile_analysis.py --mode GPU
    python analysis/error_profile_analysis.py --mode CPU

Requirements:
    edlib, biopython, pandas, numpy
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
import edlib

logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════
# USER CONFIGURATION — Modify only this block
# ═══════════════════════════════════════════════════
# Base results directory (change to CPU folder if needed)
BASE_DIR = Path("/content/drive/MyDrive/results/GPU")

# Reference FASTA files for each taxon group
REFS = {
    "aspergillus": Path("/content/drive/MyDrive/logs/references/ref_asper.fna"),
    "cladosporium": Path("/content/drive/MyDrive/logs/references/ref_clado.fna"),
    "rhizopus": Path("/content/drive/MyDrive/logs/references/ref_rhizo.fna")
}

# Maximum number of reads to process per barcode (for speed)
MAX_READS = 300

# Output directory for the TSV table
OUTPUT_DIR = Path("/content/drive/MyDrive/results")
# ═══════════════════════════════════════════════════


def identify_taxon(barcode_num: int) -> str | None:
    """Map barcode number to taxon group."""
    if 25 <= barcode_num <= 36:
        return "rhizopus"
    elif 37 <= barcode_num <= 48:
        return "aspergillus"
    elif 49 <= barcode_num <= 56:
        return "cladosporium"
    return None


def count_errors(read_seq: str, ref_seq: str):
    """Align read to reference and return S/I/D percentages."""
    aln = edlib.align(read_seq, ref_seq, task="path", mode="NW")
    cigar = aln["cigar"]
    subs = cigar.count("X")
    ins  = cigar.count("I")
    dels = cigar.count("D")
    total = subs + ins + dels
    if total == 0:
        return 0.0, 0.0, 0.0
    return subs / total * 100, ins / total * 100, dels / total * 100


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 6 — Error profile analysis (S/I/D rates)"
    )
    parser.add_argument("--base-dir", type=Path, default=BASE_DIR,
                        help="Base directory with per-barcode results")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                        help="Directory where the summary TSV will be saved")
    parser.add_argument("--mode", choices=["CPU", "GPU"], default="GPU",
                        help="Pipeline mode (for output filename prefix)")
    parser.add_argument("--max-reads", type=int, default=MAX_READS,
                        help="Maximum reads to process per barcode")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = []

    for folder in sorted(Path(args.base_dir).iterdir()):
        if not folder.is_dir():
            continue

        # Extract barcode number from folder name
        m = re.search(r"B(\d+)", folder.name)
        if not m:
            continue
        barcode_num = int(m.group(1))
        taxon = identify_taxon(barcode_num)
        if taxon is None:
            logger.warning("Unknown taxon for barcode %s — skipping", folder.name)
            continue

        ref_file = REFS.get(taxon)
        if not ref_file or not ref_file.exists():
            logger.warning("Reference FASTA not found for %s: %s", taxon, ref_file)
            continue

        # Find trimmed FASTQ
        porechop_dir = folder / "PORECHOP"
        if not porechop_dir.exists():
            continue

        fastq_files = list(porechop_dir.glob("*trimmed*.fastq"))
        if not fastq_files:
            continue
        fastq_file = fastq_files[0]

        logger.info("Processing %s (barcode %d, taxon: %s)", folder.name, barcode_num, taxon)

        # Load reference
        ref = str(next(SeqIO.parse(ref_file, "fasta")).seq)

        S_list, I_list, D_list = [], [], []

        for idx, rec in enumerate(SeqIO.parse(fastq_file, "fastq")):
            if idx >= args.max_reads:
                break
            r = str(rec.seq)
            s_pct, i_pct, d_pct = count_errors(r, ref)
            if s_pct == 0 and i_pct == 0 and d_pct == 0:
                continue
            S_list.append(s_pct)
            I_list.append(i_pct)
            D_list.append(d_pct)

        if len(S_list) == 0:
            continue

        results.append({
            "barcode": folder.name,
            "taxon": taxon,
            "n_reads": len(S_list),
            "S_mean": np.mean(S_list),
            "S_se": np.std(S_list, ddof=1) / np.sqrt(len(S_list)) if len(S_list) > 1 else 0.0,
            "I_mean": np.mean(I_list),
            "I_se": np.std(I_list, ddof=1) / np.sqrt(len(I_list)) if len(I_list) > 1 else 0.0,
            "D_mean": np.mean(D_list),
            "D_se": np.std(D_list, ddof=1) / np.sqrt(len(D_list)) if len(D_list) > 1 else 0.0,
        })

    df = pd.DataFrame(results)

    output_path = args.output_dir / f"error_profile_summary_{args.mode}.tsv"
    df.to_csv(output_path, sep="\t", index=False)

    logger.info("Error profile analysis completed successfully (mode: %s)", args.mode)
    logger.info("   → Output saved to: %s", output_path)
    print(df.head())


if __name__ == "__main__":
    main()