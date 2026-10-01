#!/usr/bin/env python3
"""
Combine the deduplicated per-locus gene models with the tree-based ingroup /
outgroup classification into one final per-species report, and write final
CDS/protein FASTA files with informative headers.

A locus is called:
    pseudogene            if it has any internal stop codon, frameshift, or
                           in-frame stop reported by miniprot
    functional_complete   LoF-free and protein length >= --complete-frac of
                           the median reference protein length (computed from
                           --reference-proteins, so this generalizes to any
                           gene without hardcoding a reference length)
    functional_partial    LoF-free but shorter than the completeness cutoff
                           (likely a scaffold-edge / assembly-gap truncation)

Usage:
    finalize_report.py <dedup.meta.tsv> <dedup.cds.fa> <dedup.prot.fa> \
                        <tree_classification.tsv> <reference_proteins.fa> \
                        <gene_name> <species> <accession> <out_prefix> \
                        [--complete-frac 0.8]
"""
import sys
import argparse
import statistics

from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dedup_meta")
    ap.add_argument("dedup_cds")
    ap.add_argument("dedup_prot")
    ap.add_argument("tree_tsv")
    ap.add_argument("reference_proteins")
    ap.add_argument("gene_name")
    ap.add_argument("species")
    ap.add_argument("accession")
    ap.add_argument("out_prefix")
    ap.add_argument("--complete-frac", type=float, default=0.8)
    args = ap.parse_args()

    meta_rows = {}
    with open(args.dedup_meta) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            meta_rows[row["mrna_id"]] = row

    tree_rows = {}
    if args.tree_tsv and args.tree_tsv != "NONE":
        with open(args.tree_tsv) as fh:
            thead = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                vals = line.rstrip("\n").split("\t")
                row = dict(zip(thead, vals))
                tree_rows[row["candidate_id"]] = row

    cds_dict = SeqIO.to_dict(SeqIO.parse(args.dedup_cds, "fasta"))
    prot_dict = SeqIO.to_dict(SeqIO.parse(args.dedup_prot, "fasta"))

    ref_lens = [len(rec.seq) for rec in SeqIO.parse(args.reference_proteins, "fasta")]
    median_ref_len = statistics.median(ref_lens) if ref_lens else 0
    complete_thresh = median_ref_len * args.complete_frac
    print(f"[finalize] median reference protein length = {median_ref_len:.0f} aa "
          f"(n={len(ref_lens)}), completeness threshold = {complete_thresh:.0f} aa", file=sys.stderr)

    out_rows = []
    for mrna_id, row in meta_rows.items():
        prot_len = int(row["prot_len_aa"])
        internal_stops = int(row["internal_stops"])
        frameshift = int(row.get("gff_frameshift", 0))
        stopcodon = int(row.get("gff_stopcodon", 0))

        tree_row = tree_rows.get(mrna_id, {})
        clade_call = tree_row.get("clade_call", "NOT_IN_TREE")

        is_lof = (internal_stops > 0) or (frameshift > 0) or (stopcodon > 0)
        if is_lof:
            functional_call = "pseudogene"
        elif prot_len >= complete_thresh:
            functional_call = "functional_complete"
        else:
            functional_call = "functional_partial"

        gene_id = f"{args.species}|{args.accession}|{row['chrom']}:{row['start']}-{row['end']}({row['strand']})"

        out_rows.append(
            {
                "gene_id": gene_id,
                "mrna_id": mrna_id,
                "species": args.species,
                "accession": args.accession,
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "best_query": row["query"],
                "identity": row["identity"],
                "n_exons": row["n_exons"],
                "cds_len_nt": row["cds_len_nt"],
                "prot_len_aa": prot_len,
                "internal_stops": internal_stops,
                "gff_frameshift": frameshift,
                "gff_stopcodon": stopcodon,
                "functional_call": functional_call,
                "clade_call": clade_call,
                "terminal_branch_length": tree_row.get("terminal_branch_length", "NA"),
                "long_branch_suspect": tree_row.get("long_branch_suspect", "NA"),
                "dist_to_ingroup": tree_row.get("dist_to_ingroup", "NA"),
                "dist_to_outgroup": tree_row.get("dist_to_outgroup", "NA"),
            }
        )

    cols = [
        "gene_id", "mrna_id", "species", "accession", "chrom", "start", "end", "strand",
        "best_query", "identity", "n_exons", "cds_len_nt", "prot_len_aa",
        "internal_stops", "gff_frameshift", "gff_stopcodon",
        "functional_call", "clade_call", "terminal_branch_length", "long_branch_suspect",
        "dist_to_ingroup", "dist_to_outgroup",
    ]

    with open(f"{args.out_prefix}.summary.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(out_rows, key=lambda r: (r["chrom"], int(r["start"]))):
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # Only a confident, unambiguous ingroup call goes in the main FASTA files.
    # AMBIGUOUS(...)/OUTGROUP_SUSPECT/NOT_IN_TREE calls are real tree-placement
    # failures, not just low-confidence versions of a real gene -- they go to
    # a separate file so they don't pollute the main candidate set. They stay
    # in summary.tsv either way, so nothing is hidden, just kept out of the
    # sequence files you'd otherwise use directly.
    confident_rows = [r for r in out_rows if r["clade_call"] == args.gene_name]
    ambiguous_rows = [r for r in out_rows if r["clade_call"] != args.gene_name]

    def write_fasta(rows, cds_path, prot_path):
        with open(cds_path, "w") as cds_out, open(prot_path, "w") as prot_out:
            for r in rows:
                mrna_id = r["mrna_id"]
                header = f">{r['gene_id']}|{r['clade_call']}|{r['functional_call']}"
                cds_out.write(f"{header}\n{cds_dict[mrna_id].seq}\n")
                prot_out.write(f"{header}\n{prot_dict[mrna_id].seq}\n")

    write_fasta(confident_rows, f"{args.out_prefix}.cds.fa", f"{args.out_prefix}.prot.fa")
    write_fasta(ambiguous_rows, f"{args.out_prefix}.ambiguous.cds.fa", f"{args.out_prefix}.ambiguous.prot.fa")

    print(f"[finalize] Wrote {len(confident_rows)} confident loci to {args.out_prefix}.{{summary.tsv,cds.fa,prot.fa}} "
          f"and {len(ambiguous_rows)} ambiguous/outgroup loci to {args.out_prefix}.ambiguous.{{cds.fa,prot.fa}}",
          file=sys.stderr)
    for r in out_rows:
        print(f"[finalize]   {r['gene_id']}  clade={r['clade_call']}  call={r['functional_call']}  "
              f"prot_len={r['prot_len_aa']}  suspect={r['long_branch_suspect']}", file=sys.stderr)


if __name__ == "__main__":
    main()
