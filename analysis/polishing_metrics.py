#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Computes polishing metrics by aligning Draft, Racon and Medaka
    consensus sequences using minimap2 (map-ont mode).
    For each comparison (Draft_vs_Racon, Racon_vs_Medaka, Draft_vs_Medaka)
    it reports:
        - Identity (%)
        - Number of indels
        - Indels per kb
    Results are saved as a TSV table for Supplementary Material and plotting.

Usage:
    python analysis/polishing_metrics_analysis.py --mode GPU
    python analysis/polishing_metrics_analysis.py --mode CPU

Requirements:
    minimap2, pandas
"""

import argparse
import subprocess
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
# Base results directory (GPU or CPU)
BASE_DIR = Path("/content/drive/MyDrive/results/GPU")   # ← Change to /CPU if needed

# Path to minimap2 binary
MINIMAP2_BIN = "/content/drive/MyDrive/apps/minimap2-2.30_x64-linux/minimap2"

# Number of threads for minimap2
THREADS = 6

# Output directory for the metrics TSV
OUTPUT_DIR = Path("/content/drive/MyDrive/results")
# ═══════════════════════════════════════════════════


def run_minimap2(ref: Path, query: Path, out_paf: Path):
    """Run minimap2 and capture PAF output."""
    cmd = [
        MINIMAP2_BIN, "-x", "map-ont", "-t", str(THREADS),
        str(ref), str(query)
    ]
    with open(out_paf, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def parse_paf(paf_file: Path):
    """
    Parse PAF file and return (identity, indels, aligned_length).
    Identity = matches / aligned_length.
    Indels   = aligned_length - matches.
    """
    matches = 0
    aln_len = 0

    with open(paf_file) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            matches += int(fields[9])
            aln_len += int(fields[10])

    if aln_len == 0:
        return None

    identity = matches / aln_len
    indels = aln_len - matches
    return identity, indels, aln_len


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 10 — Polishing metrics (Draft vs Racon vs Medaka)"
    )
    parser.add_argument("--base-dir", type=Path, default=BASE_DIR,
                        help="Base directory with per-barcode results")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                        help="Directory where the metrics TSV will be saved")
    parser.add_argument("--mode", choices=["CPU", "GPU"], default="GPU",
                        help="Pipeline mode (for output filename prefix)")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = []

    for barcode in sorted(p for p in args.base_dir.iterdir() if p.is_dir()):
        bdir = barcode

        draft  = bdir / "AMPLICON_SORTER" / "draft.fasta"
        racon  = bdir / "RACON" / "consensus_racon.fasta"
        medaka = bdir / "MEDAKA" / "consensus.fasta"

        if not (draft.exists() and racon.exists() and medaka.exists()):
            logger.warning("Missing polishing files for %s — skipping", barcode.name)
            continue

        logger.info("Processing polishing metrics for %s", barcode.name)

        comparisons = [
            ("Draft_vs_Racon", draft, racon),
            ("Racon_vs_Medaka", racon, medaka),
            ("Draft_vs_Medaka", draft, medaka),
        ]

        for label, ref, qry in comparisons:
            paf = bdir / f"{label}.paf"
            run_minimap2(ref, qry, paf)

            parsed = parse_paf(paf)
            if parsed is None:
                continue

            identity, indels, aln_len = parsed

            results.append({
                "Barcode": barcode.name,
                "Comparison": label,
                "Identity": round(identity, 4),
                "Indels": indels,
                "Aligned_length": aln_len,
                "Indels_per_kb": round(indels / (aln_len / 1000), 2),
            })

    df = pd.DataFrame(results)

    output_csv = args.output_dir / f"polishing_metrics_{args.mode}.tsv"
    df.to_csv(output_csv, sep="\t", index=False)

    logger.info("Polishing metrics analysis completed successfully (mode: %s)", args.mode)
    logger.info("   → Output saved to: %s", output_csv)


if __name__ == "__main__":
    main()