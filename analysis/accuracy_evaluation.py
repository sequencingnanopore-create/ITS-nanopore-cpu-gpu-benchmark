#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Computes taxonomic accuracy metrics (species-level and genus-level)
    against known ground-truth for each barcode. Produces two TSV files
    ready for Supplementary Material:
        - species_majority_per_barcode_{CPU/GPU}.tsv
        - accuracy_summary_{CPU/GPU}.tsv

Usage:
    python analysis/accuracy_evaluation.py --mode GPU
    python analysis/accuracy_evaluation.py --mode CPU

Requirements:
    pandas
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
# Directory containing the per-barcode results (CPU or GPU)
RESULTS_DIR = Path("/content/drive/MyDrive/results/GPU")   # ← Change to /CPU if needed

# Output directory for the generated TSV files
OUTPUT_DIR  = Path("/content/drive/MyDrive/results")

# Pipeline mode (used for output filename prefix)
MODE = "GPU"      # or "CPU"
# ═══════════════════════════════════════════════════


def extract_barcode_number(qseqid: str):
    """Extract barcode number from qseqid (e.g. HBAN14_B33)."""
    m = re.search(r"_B(\d+)", qseqid)
    return int(m.group(1)) if m else None


def main():
    parser = argparse.ArgumentParser(
        description="CODIGO 1 — Taxonomic accuracy evaluation (CPU/GPU)"
    )
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR,
                        help="Directory with per-barcode results")
    parser.add_argument("--output-dir",  type=Path, default=OUTPUT_DIR,
                        help="Directory where TSV files will be saved")
    parser.add_argument("--mode", choices=["CPU", "GPU"], default=MODE,
                        help="Pipeline mode (for output filename prefix)")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Load taxonomy file ────────────────────────────────────────
    taxonomy_file = args.results_dir / "work" / "final_taxonomy.tsv"
    if not taxonomy_file.exists():
        logger.error("final_taxonomy.tsv not found: %s", taxonomy_file)
        return

    df = pd.read_csv(taxonomy_file, sep="\t")
    df["barcode_num"] = df["qseqid"].apply(extract_barcode_number)

    # ── 2. Ground truth (species-level) ──────────────────────────────
    ground_truth = {}
    for b in [48, 47, 45, 44, 43, 41, 40, 38, 37, 46]:
        ground_truth[b] = "Aspergillus_niger"
    for b in [49, 50, 51, 52, 53, 54, 55, 56]:
        ground_truth[b] = "Cladosporium_cladosporioides"
    for b in [36, 34, 35, 33, 32, 31, 30, 29, 25, 26, 24]:
        ground_truth[b] = "Rhizopus_stolonifer"

    # ── 3. Majority species per barcode ──────────────────────────────
    rows = []
    for barcode, sub in df.groupby("barcode_num"):
        if barcode not in ground_truth:
            continue

        valid = sub["species"].dropna()
        valid = valid[(valid != "") & (valid != "Unassigned")]

        if valid.empty:
            rows.append({
                "barcode": barcode,
                "expected_species": ground_truth[barcode],
                "observed_species": "No valid species",
                "correct": False
            })
            continue

        species_counts = valid.value_counts()
        major_species = species_counts.idxmax()
        major_reads = species_counts.max()
        total_reads = len(valid)
        dominance = major_reads / total_reads

        rows.append({
            "barcode": barcode,
            "expected_species": ground_truth[barcode],
            "observed_species": major_species,
            "reads_supporting_majority": major_reads,
            "total_valid_reads": total_reads,
            "dominance_fraction": round(dominance, 3),
            "correct": major_species == ground_truth[barcode]
        })

    res = pd.DataFrame(rows).sort_values("barcode")

    # ── 4. Global accuracy metrics ───────────────────────────────────
    def get_genus(name):
        if not isinstance(name, str) or name == "No valid species":
            return None
        return name.split("_")[0]

    res["expected_genus"] = res["expected_species"].apply(get_genus)
    res["observed_genus"] = res["observed_species"].apply(get_genus)

    n_total = len(res)
    n_correct_species = res["correct"].sum()
    n_correct_genus = (
        (~res["correct"]) &
        (res["expected_genus"] == res["observed_genus"]) &
        (res["observed_genus"].notna())
    ).sum()
    n_incorrect = n_total - n_correct_species - n_correct_genus

    accuracy_summary = pd.DataFrame({
        "metric": ["Correct (species level)", "Correct (genus level)", "Incorrect"],
        "count": [n_correct_species, n_correct_genus, n_incorrect],
        "percentage": [
            round(100 * n_correct_species / n_total, 2),
            round(100 * n_correct_genus / n_total, 2),
            round(100 * n_incorrect / n_total, 2)
        ]
    })

    # ── 5. Export Supplementary Material files ───────────────────────
    prefix = args.mode
    res.to_csv(
        args.output_dir / f"species_majority_per_barcode_{prefix}.tsv",
        sep="\t", index=False
    )
    accuracy_summary.to_csv(
        args.output_dir / f"accuracy_summary_{prefix}.tsv",
        sep="\t", index=False
    )

    logger.info("Accuracy evaluation completed successfully (mode: %s)", prefix)
    logger.info("   → %s", args.output_dir / f"species_majority_per_barcode_{prefix}.tsv")
    logger.info("   → %s", args.output_dir / f"accuracy_summary_{prefix}.tsv")


if __name__ == "__main__":
    main()