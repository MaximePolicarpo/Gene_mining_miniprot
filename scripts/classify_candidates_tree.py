#!/usr/bin/env python3
"""
Given a FastTree tree built from a reference alignment (reference_proteins +
some extra outgroup sequences used only for rooting) with candidate sequences
added (mafft --add --keeplength), call each candidate as nested in the target
gene's reference set ("ingroup") or as an outgroup/suspect placement, and flag
candidates with suspiciously long terminal branches.

Ingroup / outgroup membership requires no header-naming convention: a tree
leaf is "ingroup" if its name is found in --reference-proteins, "candidate" if
found in --candidate-ids, and "outgroup" otherwise (i.e. any extra sequence
present in the reference alignment/tree but not in reference_proteins -- for
example a handful of distant paralogues added purely to root the tree).

Usage:
    classify_candidates_tree.py <treefile> <candidate_ids.txt> <reference_proteins.fa> \
                                 <gene_name> <out_tsv> [--k 3] [--branch-outlier-factor 5]
"""
import sys
import argparse
import statistics

from Bio import Phylo
from Bio import SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("treefile")
    ap.add_argument("candidate_ids")
    ap.add_argument("reference_proteins")
    ap.add_argument("gene_name")
    ap.add_argument("out_tsv")
    ap.add_argument("--k", type=int, default=3, help="number of nearest references averaged per group")
    ap.add_argument("--branch-outlier-factor", type=float, default=5.0,
                     help="flag candidate as long-branch if its terminal branch length exceeds "
                          "factor * median ingroup terminal branch length")
    ap.add_argument("--ambiguous-ratio", type=float, default=1.3,
                     help="if the closer group's mean distance is within this ratio of the "
                          "second-closest group's mean distance, call it ambiguous")
    args = ap.parse_args()

    with open(args.candidate_ids) as fh:
        candidate_ids = set(line.strip() for line in fh if line.strip())

    ingroup_ids = {rec.id for rec in SeqIO.parse(args.reference_proteins, "fasta")}

    tree = Phylo.read(args.treefile, "newick")
    terminals = tree.get_terminals()
    name_to_term = {t.name: t for t in terminals}

    ingroup_refs, outgroup_refs, candidates_in_tree = [], [], []
    for t in terminals:
        name = t.name
        if name in candidate_ids:
            candidates_in_tree.append(name)
        elif name in ingroup_ids:
            ingroup_refs.append(name)
        else:
            outgroup_refs.append(name)

    print(f"[classify] tree has {len(terminals)} leaves: "
          f"{len(ingroup_refs)} ingroup ({args.gene_name}) refs, "
          f"{len(outgroup_refs)} outgroup refs, {len(candidates_in_tree)} candidates",
          file=sys.stderr)

    ref_branch_lengths = [
        name_to_term[n].branch_length for n in ingroup_refs
        if name_to_term[n].branch_length is not None
    ]
    median_ref_bl = statistics.median(ref_branch_lengths) if ref_branch_lengths else 0.0
    bl_threshold = median_ref_bl * args.branch_outlier_factor
    print(f"[classify] median ingroup terminal branch length = {median_ref_bl:.4f}, "
          f"long-branch threshold = {bl_threshold:.4f}", file=sys.stderr)

    results = []
    for cand_name in candidates_in_tree:
        cand_term = name_to_term[cand_name]
        cand_bl = cand_term.branch_length if cand_term.branch_length is not None else 0.0

        def nearest_mean(ref_list):
            if not ref_list:
                return None
            dists = sorted(tree.distance(cand_term, name_to_term[r]) for r in ref_list)
            k = min(args.k, len(dists))
            return sum(dists[:k]) / k

        in_mean = nearest_mean(ingroup_refs)
        out_mean = nearest_mean(outgroup_refs)

        group_means = {args.gene_name: in_mean, "OUTGROUP": out_mean}
        group_means = {k: v for k, v in group_means.items() if v is not None}
        ranked = sorted(group_means.items(), key=lambda kv: kv[1])

        best_group, best_dist = ranked[0]
        second_group, second_dist = ranked[1] if len(ranked) > 1 else (None, None)

        ambiguous = False
        if second_dist is not None and best_dist > 0:
            if (second_dist / best_dist) < args.ambiguous_ratio:
                ambiguous = True

        call = best_group
        if best_group == "OUTGROUP":
            call = "OUTGROUP_SUSPECT"
        elif ambiguous:
            call = f"AMBIGUOUS({best_group}~{second_group})"

        long_branch = (bl_threshold > 0) and (cand_bl > bl_threshold)

        results.append(
            {
                "candidate_id": cand_name,
                "clade_call": call,
                "terminal_branch_length": round(cand_bl, 5),
                "long_branch_suspect": long_branch,
                "dist_to_ingroup": round(in_mean, 5) if in_mean is not None else "NA",
                "dist_to_outgroup": round(out_mean, 5) if out_mean is not None else "NA",
            }
        )

    missing = candidate_ids - set(candidates_in_tree)
    if missing:
        print(f"[classify] WARNING: {len(missing)} candidate ids not found as leaves in the tree: "
              f"{sorted(missing)}", file=sys.stderr)

    with open(args.out_tsv, "w") as out:
        cols = ["candidate_id", "clade_call", "terminal_branch_length", "long_branch_suspect",
                "dist_to_ingroup", "dist_to_outgroup"]
        out.write("\t".join(cols) + "\n")
        for r in results:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"[classify] Wrote {len(results)} candidate calls to {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
