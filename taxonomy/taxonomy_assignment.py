#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Hybrid BLAST + SINTAX taxonomic assignment pipeline.

    Collects non-chimeric consensus sequences from all processed barcodes,
    then applies a hierarchical decision framework (Rules A/B/C) that
    integrates BLAST local-alignment and SINTAX probabilistic classification.

    Decision rules:
        Rule A — BLAST species (pident ≥ 97 %, qcov ≥ 90 %) + SINTAX agreement
        Rule B — Extreme BLAST dominance (pident ≥ 99 %, qcov ≥ 95 %,
                 top-vs-second gap ≥ 1 %)
        Rule C — Genus-level consensus (BLAST pident ≥ 94 %, qcov ≥ 80 %,
                 SINTAX genus confidence ≥ 0.90)
        Fallback — LCA / weak-hit genus / "Fungi sp."

    Outputs:
        final_taxonomy.tsv         — one row per OTU
        final_taxonomy_verbose.tsv — merged BLAST+SINTAX table
        unassigned.tsv             — OTUs that fell back to Fungi sp.

Usage:
    python pipeline/taxonomy_assignment.py \
        --results-dir results/ \
        --blast-db    db/unite2024ITS \
        --unite-fasta db/UNITE_public_19.02.2025.fasta \
        --sintax-fasta db/utax_reference_dataset_19.02.2025.fasta \
        --output-dir  results/taxonomy \
        --threads 4
"""

import argparse
import glob
import logging
import math
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)

# ── Thresholds ───────────────────────────────────────────────
SPECIES_PID      = 97.0
SPECIES_QCOV     = 90.0
GENUS_PID        = 94.0
GENUS_QCOV       = 80.0
RULE_B_PID       = 99.0
RULE_B_QCOV      = 95.0
RULE_B_DELTA     = 1.0
SINTAX_SP_CUTOFF = 0.90
SINTAX_GN_CUTOFF = 0.90
BLAST_EVALUE     = 1e-5
BLAST_MAX_HITS   = 10


# ─────────────────────────────────────────────────────────────
# 1 — Gather consensus sequences from all barcodes
# ─────────────────────────────────────────────────────────────

def gather_queries(results_dir, output_fasta):
    """
    Collect UCHIME/nonchimeras.fasta from every barcode directory.
    Each sequence ID is prefixed with the barcode folder name.
    """
    count = 0
    with open(output_fasta, "w") as out:
        for bc_dir in sorted(Path(results_dir).iterdir()):
            if not bc_dir.is_dir():
                continue
            nonchim = bc_dir / "UCHIME" / "nonchimeras.fasta"
            if not nonchim.exists():
                logger.warning("Missing: %s  — skipping barcode %s",
                               nonchim, bc_dir.name)
                continue
            for rec in SeqIO.parse(nonchim, "fasta"):
                rec.id          = f"{bc_dir.name}__{rec.id}"
                rec.description = ""
                SeqIO.write(rec, out, "fasta")
                count += 1

    if count == 0:
        raise SystemExit("No OTUs collected — check that UCHIME step completed.")

    logger.info("Collected %d OTUs from %s", count, results_dir)
    return count


# ─────────────────────────────────────────────────────────────
# 2 — BLAST
# ─────────────────────────────────────────────────────────────

def run_blast(queries_fasta, blast_db, output_tsv, threads, blastn_bin="blastn"):
    fmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qcovs qcovhsp stitle"
    )
    cmd = [
        blastn_bin,
        "-query",           str(queries_fasta),
        "-db",              blast_db,
        "-outfmt",          fmt,
        "-max_target_seqs", str(BLAST_MAX_HITS),
        "-num_threads",     str(threads),
        "-evalue",          str(BLAST_EVALUE),
        "-out",             str(output_tsv),
    ]
    logger.info("[BLAST] Running …")
    subprocess.run(cmd, check=True)


def parse_blast(tsv_path, acc_map):
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qcovs", "qcovhsp", "stitle",
    ]
    df = pd.read_csv(tsv_path, sep="\t", names=cols, low_memory=False)
    df["pident"]  = pd.to_numeric(df["pident"],  errors="coerce").fillna(0.0)
    df["qcovhsp"] = pd.to_numeric(df["qcovhsp"], errors="coerce").fillna(0.0)
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce").fillna(0.0)
    df = df[df["qcovhsp"] > 1.0]

    def _map(acc):
        if pd.isna(acc):
            return ("", "")
        base = str(acc).split(".")[0]
        return acc_map.get(str(acc), acc_map.get(base, ("", "")))

    df[["blast_genus", "blast_species"]] = (
        df["sseqid"].apply(lambda x: pd.Series(_map(x)))
    )

    agg = df.groupby("qseqid").agg(
        avg_pid  =("pident",   "mean"),
        avg_qcov =("qcovhsp",  "mean"),
        n_hits   =("sseqid",   "count"),
        best_bits=("bitscore", "max"),
    ).reset_index()
    return df.merge(agg, on="qseqid", how="left")


# ─────────────────────────────────────────────────────────────
# 3 — SINTAX
# ─────────────────────────────────────────────────────────────

def run_sintax(queries_fasta, sintax_db, output_tsv, threads, vsearch_bin="vsearch"):
    cmd = [
        vsearch_bin,
        "--sintax",        str(queries_fasta),
        "--db",            sintax_db,
        "--tabbedout",     str(output_tsv),
        "--sintax_cutoff", "0.0",
        "--strand",        "both",
        "--threads",       str(threads),
    ]
    logger.info("[SINTAX] Running …")
    subprocess.run(cmd, check=True)


def parse_sintax(tsv_path):
    rows = []
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            q, annot = parts[0], parts[1]
            genus = species = ""
            gscore = sscore = 0.0
            for tok in annot.split(","):
                tok = tok.strip()
                m = re.match(r"g:([^(\s]+)(?:\(([\d.]+)\))?", tok)
                if m:
                    genus  = m.group(1)
                    gscore = float(m.group(2)) if m.group(2) else 0.0
                m = re.match(r"s:([^(\s]+)(?:\(([\d.]+)\))?", tok)
                if m:
                    species = m.group(1)
                    sscore  = float(m.group(2)) if m.group(2) else 0.0
            rows.append((q, genus, gscore, species, sscore))
    return pd.DataFrame(rows, columns=[
        "qseqid", "sintax_genus", "sintax_genus_score",
        "sintax_species", "sintax_species_score",
    ])


# ─────────────────────────────────────────────────────────────
# 4 — UNITE accession map
# ─────────────────────────────────────────────────────────────

def build_unite_map(unite_fasta):
    """Parse UNITE FASTA headers → dict {accession: (genus, species)}."""
    logger.info("[UNITE] Building accession map from %s …", unite_fasta)
    acc_map = {}
    with open(unite_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line.strip()
            genus = species = ""
            mg = re.search(r"g__([^;|]+)", header)
            ms = re.search(r"s__([^;|]+)", header)
            if mg:
                genus = mg.group(1).strip()
            if ms:
                species = ms.group(1).strip()
            for tok in re.split(r"[>\s|]+", header):
                if tok and tok not in acc_map:
                    acc_map[tok] = (genus, species)
    return acc_map


# ─────────────────────────────────────────────────────────────
# 5 — Taxonomy decision
# ─────────────────────────────────────────────────────────────

def _penalty(pid, qcov):
    try:
        return math.log1p(100 - pid) + math.log1p(100 - qcov)
    except Exception:
        return 0.0


def assign_taxonomy(merged_df):
    """Apply Rules A / B / C and fallbacks to each OTU."""
    records = []
    for q, sub in merged_df.groupby("qseqid"):
        row = dict(qseqid=q, rule="Fallback", method="Fungi_sp",
                   genus="Fungi", species="", confidence=0.0)

        # ── Rule A: BLAST species + SINTAX agreement ─────────
        sp_ok = sub[
            (sub["avg_pid"]  >= SPECIES_PID) &
            (sub["avg_qcov"] >= SPECIES_QCOV) &
            (sub["blast_species"] != "")
        ]
        if not sp_ok.empty:
            best = sp_ok.sort_values("avg_pid", ascending=False).iloc[0]
            sx   = sub.iloc[0]
            if sx["sintax_species_score"] >= SINTAX_SP_CUTOFF and \
               sx["sintax_species"] == best["blast_species"]:
                row.update(rule="A", method="BLAST_species+SINTAX",
                           genus=best["blast_genus"],
                           species=best["blast_species"],
                           confidence=best["avg_pid"])
                row["penalty"] = _penalty(best["avg_pid"], best["avg_qcov"])
                records.append(row)
                continue

        # ── Rule B: Extreme BLAST dominance ──────────────────
        hits_sorted = sub.sort_values("bitscore", ascending=False)
        if len(hits_sorted) >= 2:
            top  = hits_sorted.iloc[0]
            sec  = hits_sorted.iloc[1]
            if (top["pident"] >= RULE_B_PID and
                    top["qcovhsp"] >= RULE_B_QCOV and
                    (top["pident"] - sec["pident"]) >= RULE_B_DELTA):
                row.update(rule="B", method="BLAST_dominant",
                           genus=top["blast_genus"],
                           species=top["blast_species"],
                           confidence=top["pident"])
                row["penalty"] = _penalty(top["pident"], top["qcovhsp"])
                records.append(row)
                continue

        # ── Rule C: Genus-level consensus ────────────────────
        gn_ok = sub[
            (sub["avg_pid"]  >= GENUS_PID) &
            (sub["avg_qcov"] >= GENUS_QCOV) &
            (sub["blast_genus"] != "")
        ]
        sx = sub.iloc[0]
        if not gn_ok.empty and sx["sintax_genus_score"] >= SINTAX_GN_CUTOFF:
            best = gn_ok.sort_values("avg_pid", ascending=False).iloc[0]
            if sx["sintax_genus"] == best["blast_genus"]:
                row.update(rule="C", method="genus_consensus",
                           genus=best["blast_genus"], species="",
                           confidence=best["avg_pid"])
                row["penalty"] = _penalty(best["avg_pid"], best["avg_qcov"])
                records.append(row)
                continue

        # ── Fallback ─────────────────────────────────────────
        row["penalty"] = 99.0
        records.append(row)

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Hybrid BLAST+SINTAX ITS taxonomy assignment (Rules A/B/C)"
    )
    parser.add_argument("--results-dir",  required=True,
                        help="Directory containing per-barcode subdirectories")
    parser.add_argument("--blast-db",     required=True,
                        help="BLAST database prefix (UNITE-derived)")
    parser.add_argument("--unite-fasta",  required=True,
                        help="UNITE release FASTA (for accession → genus/species map)")
    parser.add_argument("--sintax-fasta", required=True,
                        help="UNITE UTAX reference FASTA for SINTAX")
    parser.add_argument("--output-dir",   required=True,
                        help="Output directory")
    parser.add_argument("--threads",      type=int, default=4)
    parser.add_argument("--blastn-bin",   default="blastn")
    parser.add_argument("--vsearch-bin",  default="vsearch")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    queries    = out / "otu_queries.fasta"
    blast_tsv  = out / "blast_results.tsv"
    sintax_tsv = out / "sintax_results.tsv"

    # 1. Gather
    gather_queries(args.results_dir, queries)

    # 2. BLAST
    run_blast(queries, args.blast_db, blast_tsv, args.threads, args.blastn_bin)

    # 3. SINTAX
    run_sintax(queries, args.sintax_fasta, sintax_tsv, args.threads, args.vsearch_bin)

    # 4. Parse & merge
    acc_map   = build_unite_map(args.unite_fasta)
    blast_df  = parse_blast(blast_tsv, acc_map)
    sintax_df = parse_sintax(sintax_tsv)
    merged    = blast_df.merge(sintax_df, on="qseqid", how="left")

    # 5. Assign
    final = assign_taxonomy(merged)

    # 6. Save
    merged.to_csv(out / "final_taxonomy_verbose.tsv", sep="\t", index=False)
    final.to_csv(out / "final_taxonomy.tsv",          sep="\t", index=False)

    unassigned = final[final["method"] == "Fungi_sp"]
    if not unassigned.empty:
        unassigned.to_csv(out / "unassigned.tsv", sep="\t", index=False)

    logger.info("Final taxonomy → %s", out / "final_taxonomy.tsv")
    logger.info("  Total OTUs: %d", len(final))
    logger.info("  Unassigned: %d", len(unassigned))


if __name__ == "__main__":
    main()
