#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Generates MD5 checksums for all basecalled FASTQ files (one per barcode)
    located in the BASECALLED subfolder. The resulting CSV is used for data
    integrity verification and is included in the repository for reproducibility
    (Supplementary Material).

Usage:
    python analysis/md5_checksums.py --mode GPU
    python analysis/md5_checksums.py --mode CPU

Requirements:
    hashlib, csv, pathlib
"""

import argparse
import csv
import hashlib
from pathlib import Path

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

# Output CSV file for MD5 checksums
OUTPUT_CSV = Path("/content/drive/MyDrive/results/fastq_md5s.csv")
# ═══════════════════════════════════════════════════


def compute_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Compute MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 13 — MD5 checksum generation for basecalled FASTQ files"
    )
    parser.add_argument("--base-dir", type=Path, default=BASE_DIR,
                        help="Base directory with per-barcode results")
    parser.add_argument("--output-csv", type=Path, default=OUTPUT_CSV,
                        help="Path to output CSV file")
    parser.add_argument("--mode", choices=["CPU", "GPU"], default="GPU",
                        help="Pipeline mode (for logging only)")
    args = parser.parse_args()

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["barcode", "file_name", "md5"])

        processed = 0
        for barcode_folder in sorted(p for p in args.base_dir.iterdir() if p.is_dir()):
            fastq_path = barcode_folder / "BASECALLED" / f"{barcode_folder.name}.fastq"

            if not fastq_path.is_file():
                logger.warning("FASTQ not found for %s — skipping", barcode_folder.name)
                continue

            md5_hash = compute_md5(fastq_path)
            writer.writerow([barcode_folder.name, fastq_path.name, md5_hash])
            processed += 1

            logger.info("✓ %s: %s", barcode_folder.name, md5_hash)

    logger.info("MD5 checksum generation completed successfully (mode: %s)", args.mode)
    logger.info("   → Processed %d barcodes", processed)
    logger.info("   → Output saved to: %s", args.output_csv)


if __name__ == "__main__":
    main()