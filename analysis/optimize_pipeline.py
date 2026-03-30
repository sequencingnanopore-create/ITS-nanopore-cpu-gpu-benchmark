#!/usr/bin/env python3
"""
Authors: Daniel Albuja, Peter Zambrano, Paolo Maldonado,
         J. Riczury Olmos-Tovar, Edwin Vera
Date: 2026-03-28
Description:
    Bayesian hyperparameter optimization for VSEARCH-based ITS clustering.
    Called by run_cpu_pipeline.py; can also be run standalone.

    Five hyperparameters are optimized per barcode:
        1. min_length        — minimum read length (400–1000 bp)
        2. min_quality       — minimum mean Phred quality (10–20)
        3. cluster_identity  — VSEARCH clustering identity (0.85–0.99)
        4. min_cluster_size  — minimum reads per cluster for consensus (2–5)
        5. consensus_identity — majority-vote threshold in MUSCLE alignments (0.4–0.8)

    Objective (minimized):
        Score = P_ambiguity + P_singleton + P_diversity + P_BLAST + P_structure

    Best trial's polished consensus is saved as polished_best.fasta.
    Optuna trial history is exported as optuna_trials.json.

Usage:
    python pipeline/optimize_pipeline.py \
        --input trimmed/B33_trimmed.fastq \
        --output results/B33/VSEARCH \
        --blast-db db/unite2024ITS \
        --trials 12 \
        --threads 4
"""

import argparse
import json
import logging
import math
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import optuna
from Bio import SeqIO

optuna.logging.set_verbosity(optuna.logging.WARNING)

logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────

def run_cmd(cmd, stdout_file=None, check=True):
    if stdout_file is not None:
        with open(stdout_file, "w") as fh:
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result


def count_ambiguous(fasta_path):
    """Return mean % of ambiguous (N) bases across sequences in a FASTA."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return 100.0
    pcts = []
    for rec in records:
        seq = str(rec.seq).upper()
        n_count = seq.count("N")
        pcts.append(100.0 * n_count / max(len(seq), 1))
    return float(np.mean(pcts))


def shannon_entropy(counts):
    """Shannon entropy of a list of counts (bits)."""
    total = sum(counts)
    if total == 0:
        return 0.0
    probs = [c / total for c in counts if c > 0]
    return -sum(p * math.log2(p) for p in probs)


def blast_penalty(fasta_path, blast_db, blastn_bin, threads):
    """
    BLAST-based penalty: smooth log function of deviation from perfect
    identity (100 %) and query coverage (100 %).
    Returns infinity when no hits are found.
    """
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        return float("inf")

    result = subprocess.run(
        [blastn_bin,
         "-query", str(fasta_path),
         "-db", blast_db,
         "-outfmt", "6 pident qcovhsp",
         "-max_target_seqs", "1",
         "-num_threads", str(threads),
         "-evalue", "1e-5"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    lines = result.stdout.decode().strip().splitlines()
    if not lines:
        return float("inf")

    penalties = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        try:
            pid  = float(parts[0])
            qcov = float(parts[1])
            pen  = math.log1p(100 - pid) + math.log1p(100 - qcov)
            penalties.append(pen)
        except ValueError:
            continue

    return float(np.mean(penalties)) if penalties else float("inf")


def majority_consensus(alignment_fasta, threshold):
    """
    Build a consensus from a MUSCLE alignment using a majority-vote rule.
    Positions where no base meets `threshold` are written as 'N'.
    """
    records = list(SeqIO.parse(alignment_fasta, "fasta"))
    if not records:
        return ""

    length = len(str(records[0].seq))
    consensus = []
    for i in range(length):
        col = [str(rec.seq)[i].upper() for rec in records]
        col = [b for b in col if b != "-"]
        if not col:
            consensus.append("N")
            continue
        counts = Counter(col)
        top_base, top_count = counts.most_common(1)[0]
        if top_count / len(col) >= threshold:
            consensus.append(top_base)
        else:
            consensus.append("N")
    return "".join(consensus)


# ─────────────────────────────────────────────────────────────
# Single trial execution
# ─────────────────────────────────────────────────────────────

def run_trial(trial_dir, reads_fastq, params, bins, blast_db, threads):
    """
    Execute one full clustering + polishing trial.
    Returns (score, polished_fasta_path).  score = inf when the trial fails.
    """
    trial_dir.mkdir(parents=True, exist_ok=True)

    min_length         = params["min_length"]
    min_quality        = params["min_quality"]
    cluster_identity   = params["cluster_identity"]
    min_cluster_size   = params["min_cluster_size"]
    consensus_identity = params["consensus_identity"]

    # ── 1. Read filtering ────────────────────────────────────
    filtered = trial_dir / "filtered.fastq"
    run_cmd(
        [bins["vsearch"], "--fastq_filter", str(reads_fastq),
         "--fastq_minlen", str(min_length),
         "--fastq_qmin",   str(min_quality),
         "--fastqout", str(filtered)],
    )
    if not filtered.exists() or filtered.stat().st_size == 0:
        return float("inf"), None

    # ── 2. Dereplication ─────────────────────────────────────
    derep = trial_dir / "derep.fasta"
    run_cmd(
        [bins["vsearch"], "--derep_fulllength", str(filtered),
         "--output", str(derep), "--sizeout"],
    )

    # ── 3. De-novo chimera removal ───────────────────────────
    nochim = trial_dir / "nochimeras.fasta"
    run_cmd(
        [bins["vsearch"], "--uchime_denovo", str(derep),
         "--nonchimeras", str(nochim)],
    )
    if not nochim.exists() or nochim.stat().st_size == 0:
        return float("inf"), None

    # ── 4. Clustering ────────────────────────────────────────
    centroids = trial_dir / "centroids.fasta"
    uc_file   = trial_dir / "clusters.uc"
    run_cmd(
        [bins["vsearch"], "--cluster_size", str(nochim),
         "--id", str(cluster_identity),
         "--strand", "both",
         "--centroids", str(centroids),
         "--uc", str(uc_file),
         "--sizeout"],
    )
    if not centroids.exists() or centroids.stat().st_size == 0:
        return float("inf"), None

    # Parse cluster sizes from .uc file
    cluster_sizes: dict[str, int] = {}
    with open(uc_file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if parts[0] == "C":
                cluster_id = parts[1]
                size = int(parts[2])
                cluster_sizes[cluster_id] = size

    if not cluster_sizes:
        return float("inf"), None

    # ── 5. Consensus generation ──────────────────────────────
    consensus_seqs = []
    for cluster_id, size in cluster_sizes.items():
        if size < min_cluster_size:
            continue

        # Extract reads belonging to this cluster (re-parse .uc)
        read_ids = set()
        with open(uc_file) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if parts[0] == "H" and parts[1] == cluster_id:
                    read_ids.add(parts[8])

        if not read_ids:
            continue

        cluster_fasta = trial_dir / f"cluster_{cluster_id}.fasta"
        with open(cluster_fasta, "w") as out:
            for rec in SeqIO.parse(filtered, "fastq"):
                if rec.id in read_ids:
                    out.write(f">{rec.id}\n{rec.seq}\n")

        if not cluster_fasta.exists() or cluster_fasta.stat().st_size == 0:
            continue

        # MUSCLE alignment
        aln = trial_dir / f"cluster_{cluster_id}_aln.fasta"
        run_cmd([bins["muscle"], "-in", str(cluster_fasta), "-out", str(aln)],
                check=False)
        if not aln.exists():
            continue

        consensus_seq = majority_consensus(aln, consensus_identity)
        if consensus_seq and consensus_seq.count("N") / len(consensus_seq) < 0.5:
            consensus_seqs.append((f"cluster_{cluster_id}_size{size}", consensus_seq))

    if not consensus_seqs:
        return float("inf"), None

    consensus_fasta = trial_dir / "consensus.fasta"
    with open(consensus_fasta, "w") as out:
        for name, seq in consensus_seqs:
            out.write(f">{name}\n{seq}\n")

    # ── 6. Racon polishing (single round) ────────────────────
    paf = trial_dir / "overlaps.paf"
    run_cmd(
        [bins["minimap2"], "-x", "map-ont", "-t", str(threads),
         str(consensus_fasta), str(filtered)],
        stdout_file=paf,
    )
    polished = trial_dir / "polished.fasta"
    run_cmd(
        [bins["racon"], "-t", str(threads),
         str(filtered), str(paf), str(consensus_fasta)],
        stdout_file=polished,
    )
    if not polished.exists() or polished.stat().st_size == 0:
        return float("inf"), None

    # ── 7. Objective score ───────────────────────────────────
    sizes = list(cluster_sizes.values())

    # Ambiguity penalty
    p_ambig = count_ambiguous(polished)

    # Singleton penalty
    n_singletons = sum(1 for s in sizes if s == 1)
    p_singleton = n_singletons / max(len(sizes), 1)

    # Diversity penalty (normalized Shannon entropy)
    h = shannon_entropy(sizes)
    h_max = math.log2(len(sizes)) if len(sizes) > 1 else 1.0
    p_diversity = h / h_max if h_max > 0 else 0.0

    # BLAST penalty
    p_blast = blast_penalty(polished, blast_db, bins["blastn"], threads)
    if math.isinf(p_blast):
        p_blast = 10.0  # finite fallback

    # Structural penalty: penalise if fewer than 3 clusters have > 1 read
    n_supported = sum(1 for s in sizes if s > 1)
    p_struct = 1.0 if n_supported < 3 else 0.0

    score = p_ambig + p_singleton + p_diversity + p_blast + p_struct
    return score, polished


# ─────────────────────────────────────────────────────────────
# Optuna objective
# ─────────────────────────────────────────────────────────────

def make_objective(reads_fastq, output_dir, bins, blast_db, threads):
    def objective(trial):
        params = {
            "min_length":         trial.suggest_int("min_length", 400, 1000),
            "min_quality":        trial.suggest_int("min_quality", 10, 20),
            "cluster_identity":   trial.suggest_float("cluster_identity", 0.85, 0.99),
            "min_cluster_size":   trial.suggest_int("min_cluster_size", 2, 5),
            "consensus_identity": trial.suggest_float("consensus_identity", 0.4, 0.8),
        }
        trial_dir = output_dir / f"trial_{trial.number}"
        score, _ = run_trial(trial_dir, reads_fastq, params, bins, blast_db, threads)
        return score if not math.isinf(score) else 1e9
    return objective


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Optuna hyperparameter optimization for VSEARCH ITS clustering"
    )
    parser.add_argument("--input",       required=True, help="Trimmed FASTQ")
    parser.add_argument("--output",      required=True, help="Output directory")
    parser.add_argument("--blast-db",    required=True, help="BLAST DB prefix")
    parser.add_argument("--trials",      type=int, default=12)
    parser.add_argument("--threads",     type=int, default=4)
    parser.add_argument("--seed",        type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")
    parser.add_argument("--vsearch-bin",  default="vsearch")
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--racon-bin",    default="racon")
    parser.add_argument("--muscle-bin",   default="muscle")
    parser.add_argument("--blastn-bin",   default="blastn")
    args = parser.parse_args()

    reads_fastq = Path(args.input)
    output_dir  = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not reads_fastq.exists():
        logger.error("Input FASTQ not found: %s", reads_fastq)
        sys.exit(1)

    bins = {
        "vsearch":  args.vsearch_bin,
        "minimap2": args.minimap2_bin,
        "racon":    args.racon_bin,
        "muscle":   args.muscle_bin,
        "blastn":   args.blastn_bin,
    }

    sampler = optuna.samplers.TPESampler(seed=args.seed)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    study.optimize(
        make_objective(reads_fastq, output_dir, bins, args.blast_db, args.threads),
        n_trials=args.trials,
    )

    # ── Export trial history ─────────────────────────────────
    trial_records = []
    for t in study.trials:
        trial_records.append({
            "trial":  t.number,
            "value":  t.value if t.value is not None else None,
            "params": t.params,
            "state":  str(t.state),
        })
    history_path = output_dir / "optuna_trials.json"
    with open(history_path, "w") as fh:
        json.dump(trial_records, fh, indent=2)
    logger.info("Optuna trial history → %s", history_path)

    # ── Copy best polished consensus ─────────────────────────
    best_trial = study.best_trial
    logger.info(
        "Best trial: %d  |  score: %.4f  |  params: %s",
        best_trial.number, best_trial.value, best_trial.params,
    )
    best_polished = output_dir / f"trial_{best_trial.number}" / "polished.fasta"
    if not best_polished.exists():
        logger.error("Best trial polished FASTA not found: %s", best_polished)
        sys.exit(1)

    final_path = output_dir / "polished_best.fasta"
    shutil.copy(best_polished, final_path)
    logger.info("Best polished consensus → %s", final_path)


if __name__ == "__main__":
    main()
