#!/usr/bin/env python3
"""
Collapse overlapping miniprot gene models (same locus hit by several reference
queries) down to one best model per locus.

Two models are considered the same locus if they are on the same scaffold and
their genomic intervals overlap (by at least 1 bp) or are within --merge-dist
bp of each other.

Best model = highest native miniprot alignment score (column 6 of the GFF;
this already balances query coverage, identity, gaps and introns, so it
correctly prefers a full-length, slightly-lower-identity hit over a
truncated near-identical self-hit), ties broken by longer CDS, then fewer
internal stop codons.

Usage:
    dedup_loci.py <meta.tsv> <cds.fa> <prot.fa> <out_prefix> [--merge-dist N]
"""
import sys
import argparse
from collections import defaultdict

from Bio import SeqIO


def load_fasta_dict(path):
    return SeqIO.to_dict(SeqIO.parse(path, "fasta"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta_tsv")
    ap.add_argument("cds_fa")
    ap.add_argument("prot_fa")
    ap.add_argument("out_prefix")
    ap.add_argument("--merge-dist", type=int, default=5000,
                     help="merge gene models within this many bp of each other on the same scaffold")
    ap.add_argument("--min-identity", type=float, default=0.25,
                     help="discard gene models below this miniprot identity before clustering loci")
    args = ap.parse_args()

    rows = []
    with open(args.meta_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["identity"] = float(row["identity"]) if row["identity"] not in ("NA", "") else 0.0
            row["score"] = float(row["score"]) if row.get("score", "NA") not in ("NA", "") else 0.0
            row["cds_len_nt"] = int(row["cds_len_nt"])
            row["internal_stops"] = int(row["internal_stops"])
            rows.append(row)

    print(f"[dedup] {len(rows)} raw gene models loaded", file=sys.stderr)

    rows = [r for r in rows if r["identity"] >= args.min_identity]
    print(f"[dedup] {len(rows)} gene models retained after identity >= {args.min_identity} filter", file=sys.stderr)

    # group by scaffold, sort by start, merge overlapping/nearby intervals
    by_chrom = defaultdict(list)
    for row in rows:
        by_chrom[row["chrom"]].append(row)

    loci = []  # list of list-of-rows
    for chrom, chrom_rows in by_chrom.items():
        chrom_rows.sort(key=lambda r: r["start"])
        current_cluster = []
        current_end = None
        for row in chrom_rows:
            if current_cluster and row["start"] <= current_end + args.merge_dist:
                current_cluster.append(row)
                current_end = max(current_end, row["end"])
            else:
                if current_cluster:
                    loci.append(current_cluster)
                current_cluster = [row]
                current_end = row["end"]
        if current_cluster:
            loci.append(current_cluster)

    print(f"[dedup] {len(loci)} distinct loci after merging", file=sys.stderr)

    def sort_key(r):
        return (r["score"], r["cds_len_nt"], -r["internal_stops"], r["identity"])

    best_rows = []
    for cluster in loci:
        best = max(cluster, key=sort_key)
        best["n_models_at_locus"] = len(cluster)
        best_rows.append(best)

    cds_dict = load_fasta_dict(args.cds_fa)
    prot_dict = load_fasta_dict(args.prot_fa)

    with open(f"{args.out_prefix}.cds.fa", "w") as cds_out, \
         open(f"{args.out_prefix}.prot.fa", "w") as prot_out, \
         open(f"{args.out_prefix}.meta.tsv", "w") as meta_out:

        extra_cols = ["n_models_at_locus"]
        meta_out.write("\t".join(header + extra_cols) + "\n")
        for row in sorted(best_rows, key=lambda r: (r["chrom"], r["start"])):
            mrna_id = row["mrna_id"]
            cds_out.write(f">{mrna_id}\n{cds_dict[mrna_id].seq}\n")
            prot_out.write(f">{mrna_id}\n{prot_dict[mrna_id].seq}\n")
            meta_out.write("\t".join(str(row[c]) for c in header + extra_cols) + "\n")

    print(f"[dedup] Wrote {len(best_rows)} best-per-locus models to {args.out_prefix}.*", file=sys.stderr)


if __name__ == "__main__":
    main()
