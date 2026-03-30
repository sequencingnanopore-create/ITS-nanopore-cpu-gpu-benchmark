#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Scans all barcode folders and selects one representative barcode per
    taxonomic group (Rhizopus, Aspergillus, Cladosporium) based on the
    median number of reads in the trimmed FASTQ file.
    The selected representatives are used for Optuna convergence analysis
    and error-profile pilot runs (Supplementary Material).

Usage:
    python analysis/representative_barcodes.py --mode GPU
    python analysis/representative_barcodes.py --mode CPU

Requirements:
    pandas, pathlib
"""

import argparse
import re
from pathlib import Path

import pandas as pd

logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════
# USER CONFIGURATION — Modify only this block
# ═══════════════════════════════════════════════════
# Base directory containing the per-barcode results
BASE_DIR = Path("/content/drive/MyDrive/results/GPU")   # ← Change to /CPU if needed

# Output directory for the representative table
OUTPUT_DIR = Path("/content/drive/MyDrive/results")
# ═══════════════════════════════════════════════════


def assign_group(barcode_num: int) -> str:
    """Map barcode number to taxonomic group."""
    if 25 <= barcode_num <= 36:
        return "Rhizopus"
    elif 37 <= barcode_num <= 48:
        return "Aspergillus"
    elif 49 <= barcode_num <= 56:
        return "Cladosporium"
    return "Unknown"


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 11 — Representative barcode selection by median read count"
    )
    parser.add_argument("--base-dir", type=Path, default=BASE_DIR,
                        help="Base directory with per-barcode results")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                        help="Directory where the representative table will be saved")
    parser.add_argument("--mode", choices=["CPU", "GPU"], default="GPU",
                        help="Pipeline mode (for output filename prefix)")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    records = []

    for folder in sorted(p for p in args.base_dir.iterdir() if p.is_dir()):
        porechop_dir = folder / "PORECHOP"
        if not porechop_dir.exists():
            continue

        # Find trimmed FASTQ
        fastq_files = list(porechop_dir.glob("*_trimmed*.fastq"))
        if not fastq_files:
            continue
        fastq_path = fastq_files[0]

        # Count reads (FASTQ = 4 lines per read)
        with open(fastq_path) as fh:
            line_count = sum(1 for _ in fh)
        read_count = line_count // 4

        # Extract barcode number
        m = re.search(r"B(\d+)", folder.name)
        if not m:
            continue
        barcode_num = int(m.group(1))

        group = assign_group(barcode_num)

        records.append({
            "barcode": folder.name,
            "barcode_number": barcode_num,
            "group": group,
            "read_count": read_count,
            "fastq_file": fastq_path.name,
        })

    if not records:
        logger.error("No barcodes found — check BASE_DIR configuration.")
        return

    df = pd.DataFrame(records)

    # Find representative per group (closest to median read count)
    representatives = {}
    for grp in df["group"].unique():
        if grp == "Unknown":
            continue
        subset = df[df["group"] == grp].copy()
        median_val = subset["read_count"].median()
        # Find row with smallest absolute deviation from median
        rep_idx = (subset["read_count"] - median_val).abs().idxmin()
        rep_barcode = subset.loc[rep_idx, "barcode"]
        representatives[grp] = rep_barcode

    # Add representative flag
    df["is_representative"] = df["barcode"].isin(representatives.values())

    # Save full table
    output_path = args.output_dir / f"representative_barcodes_{args.mode}.tsv"
    df.to_csv(output_path, sep="\t", index=False)

    # Save compact representative list for easy reference
    rep_df = df[df["is_representative"]][["group", "barcode", "read_count"]]
    rep_df.to_csv(args.output_dir / f"representatives_summary_{args.mode}.tsv", sep="\t", index=False)

    logger.info("Representative barcode selection completed successfully (mode: %s)", args.mode)
    logger.info("   → Full table: %s", output_path)
    logger.info("   → Summary of representatives: %s", args.output_dir / f"representatives_summary_{args.mode}.tsv")
    logger.info("Selected representatives:\n%s", rep_df.to_string(index=False))


if __name__ == "__main__":
    main()