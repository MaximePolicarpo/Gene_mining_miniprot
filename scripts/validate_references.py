#!/usr/bin/env python3
"""
Sanity-check that --reference_proteins, --reference_alignment and
--reference_tree are mutually consistent before spending any time on
miniprot/mafft/FastTree:
  - every reference_proteins id must be present in reference_alignment
    (reference_alignment = reference_proteins + optional outgroup sequences)
  - reference_tree's leaf set must exactly match reference_alignment's
    sequence set (it should be the tree already built from that alignment)

Usage:
    validate_references.py <reference_proteins.fa> <reference_alignment.fa> <reference_tree>
Exits non-zero with a clear message on any mismatch.
"""
import sys

from Bio import SeqIO, Phylo


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <reference_proteins.fa> <reference_alignment.fa> <reference_tree>")

    ref_prot_path, ref_aln_path, ref_tree_path = sys.argv[1:4]

    ref_prot_ids = {rec.id for rec in SeqIO.parse(ref_prot_path, "fasta")}
    ref_aln_ids = {rec.id for rec in SeqIO.parse(ref_aln_path, "fasta")}
    tree = Phylo.read(ref_tree_path, "newick")
    tree_leaf_names = {t.name for t in tree.get_terminals()}

    errors = []

    missing_from_aln = ref_prot_ids - ref_aln_ids
    if missing_from_aln:
        errors.append(
            f"{len(missing_from_aln)} reference_proteins id(s) are missing from reference_alignment, "
            f"e.g.: {sorted(missing_from_aln)[:5]}"
        )

    only_in_aln_not_tree = ref_aln_ids - tree_leaf_names
    only_in_tree_not_aln = tree_leaf_names - ref_aln_ids
    if only_in_aln_not_tree:
        errors.append(
            f"{len(only_in_aln_not_tree)} reference_alignment sequence(s) are missing from reference_tree "
            f"(stale tree?), e.g.: {sorted(only_in_aln_not_tree)[:5]}"
        )
    if only_in_tree_not_aln:
        errors.append(
            f"{len(only_in_tree_not_aln)} reference_tree leaf(ves) are missing from reference_alignment "
            f"(stale alignment?), e.g.: {sorted(only_in_tree_not_aln)[:5]}"
        )

    n_outgroup = len(ref_aln_ids - ref_prot_ids)
    print(
        f"[validate_references] {len(ref_prot_ids)} reference_proteins (ingroup), "
        f"{len(ref_aln_ids)} reference_alignment sequences ({n_outgroup} outgroup-only), "
        f"{len(tree_leaf_names)} reference_tree leaves",
        file=sys.stderr,
    )

    if errors:
        print("[validate_references] ERROR: reference files are inconsistent:", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    print("[validate_references] OK: reference_proteins / reference_alignment / reference_tree are consistent",
          file=sys.stderr)


if __name__ == "__main__":
    main()
