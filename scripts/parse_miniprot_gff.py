#!/usr/bin/env python3
"""
Parse a miniprot --gff output file, extract CDS (nucleotide) and protein
sequences for every predicted gene model (mRNA), and write a metadata table.

miniprot --gff output encodes, per mRNA, attributes we rely on:
  - Target=<query_protein> <qstart> <qend>   (1-based, inclusive, aa coords)
  - Rank=<n>
  - Identity=<float 0-1>
  - Positive=<float 0-1>
  - Frameshift events are stop codons / frameshifts already spliced out by
    miniprot *unless* they are real frameshifts in the genome, which miniprot
    reports as a stop codon ('*') or a frameshift marker in the translated
    protein it would produce; here we translate the extracted CDS ourselves
    and record internal stop codons / non-3n length as evidence of a
    pseudogene / frameshift.

Usage:
    parse_miniprot_gff.py <gff> <genome.fa> <out_prefix>

Writes:
    <out_prefix>.cds.fa    nucleotide CDS per mRNA (spliced, in mRNA orientation)
    <out_prefix>.prot.fa   translated protein per mRNA (stops kept as '*')
    <out_prefix>.meta.tsv  one row per mRNA with coordinates + miniprot stats
"""
import sys
import re
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


def parse_attributes(attr_str):
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k] = v
    return attrs


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <gff> <genome.fa> <out_prefix>")

    gff_path, genome_path, out_prefix = sys.argv[1:4]

    print(f"[parse_gff] Loading genome {genome_path} ...", file=sys.stderr)
    genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    print(f"[parse_gff] {len(genome)} scaffolds loaded.", file=sys.stderr)

    mrna_info = {}          # mrna_id -> dict of attrs/coords
    mrna_cds_exons = defaultdict(list)  # mrna_id -> list of (start0, end, strand, chrom)

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_str = cols
            attrs = parse_attributes(attr_str)

            if feature == "mRNA":
                mrna_id = attrs.get("ID", attrs.get("Parent", f"mRNA_{len(mrna_info)+1}"))
                target = attrs.get("Target", "")
                target_fields = target.split()
                query_name = target_fields[0] if target_fields else "NA"
                mrna_info[mrna_id] = {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "query": query_name,
                    "identity": attrs.get("Identity", "NA"),
                    "positive": attrs.get("Positive", "NA"),
                    "rank": attrs.get("Rank", "NA"),
                    "score": float(score) if score not in (".", "") else 0.0,
                    "gff_frameshift": int(attrs.get("Frameshift", 0)),
                    "gff_stopcodon": int(attrs.get("StopCodon", 0)),
                }
            elif feature == "CDS":
                mrna_id = attrs.get("Parent", None)
                if mrna_id is None:
                    continue
                mrna_cds_exons[mrna_id].append((int(start) - 1, int(end), strand, chrom))

    print(f"[parse_gff] {len(mrna_info)} mRNA models found.", file=sys.stderr)

    cds_out = open(f"{out_prefix}.cds.fa", "w")
    prot_out = open(f"{out_prefix}.prot.fa", "w")
    meta_out = open(f"{out_prefix}.meta.tsv", "w")
    meta_out.write(
        "\t".join(
            [
                "mrna_id", "chrom", "start", "end", "strand", "query",
                "identity", "positive", "rank", "score", "n_exons", "cds_len_nt",
                "prot_len_aa", "internal_stops", "ends_with_stop", "multiple_of_3",
                "gff_frameshift", "gff_stopcodon",
            ]
        )
        + "\n"
    )

    n_written = 0
    for mrna_id, info in mrna_info.items():
        exons = mrna_cds_exons.get(mrna_id)
        if not exons:
            continue
        chrom = info["chrom"]
        if chrom not in genome:
            print(f"[parse_gff] WARNING: {chrom} not found in genome fasta, skipping {mrna_id}", file=sys.stderr)
            continue
        strand = info["strand"]
        # sort exons by genomic start, concatenate respecting strand
        exons_sorted = sorted(exons, key=lambda e: e[0])
        seq_parts = []
        for (s0, e, strand_exon, _chrom) in exons_sorted:
            seq_parts.append(genome[chrom].seq[s0:e])
        cds_seq = sum(seq_parts[1:], seq_parts[0]) if seq_parts else Seq("")
        if strand == "-":
            cds_seq = cds_seq.reverse_complement()

        cds_len = len(cds_seq)
        multiple_of_3 = (cds_len % 3 == 0)

        # translate (trim to full codons to avoid Biopython warning)
        trim_len = cds_len - (cds_len % 3)
        prot_seq = str(cds_seq[:trim_len].translate())

        ends_with_stop = prot_seq.endswith("*")
        body = prot_seq[:-1] if ends_with_stop else prot_seq
        internal_stops = body.count("*")

        # miniprot's CDS features give the correct genomic span for a locus
        # even when it contains a frameshift, but a frameshift re-synchronises
        # the reading frame *within* a single reported CDS block (see
        # Frameshift= in the GFF) rather than splitting it in two -- naively
        # translating straight through one therefore produces real amino
        # acids up to the frameshift point and out-of-frame garbage after it,
        # almost always hitting a spurious stop codon within a few codons.
        # Whatever we write here is what later gets aligned onto the
        # reference alignment and placed on the tree (see dedup_loci.py /
        # mine_gene.sh step 6), so garbage downstream of a stop would corrupt
        # that alignment for no benefit -- the classification step already
        # calls this locus a pseudogene from internal_stops/frameshift alone,
        # it doesn't need the garbled tail. Keep only the clean prefix up to
        # (not including) the first stop codon, whether that stop is the
        # natural end of a normal gene or a premature one in a pseudogene.
        clean_prot_seq = prot_seq.split("*")[0]

        header = f"{mrna_id}"
        cds_out.write(f">{header}\n{cds_seq}\n")
        prot_out.write(f">{header}\n{clean_prot_seq}\n")

        meta_out.write(
            "\t".join(
                map(
                    str,
                    [
                        mrna_id, chrom, info["start"], info["end"], strand, info["query"],
                        info["identity"], info["positive"], info["rank"], info["score"],
                        len(exons_sorted), cds_len, len(clean_prot_seq),
                        internal_stops, ends_with_stop, multiple_of_3,
                        info["gff_frameshift"], info["gff_stopcodon"],
                    ],
                )
            )
            + "\n"
        )
        n_written += 1

    cds_out.close()
    prot_out.close()
    meta_out.close()
    print(f"[parse_gff] Wrote {n_written} gene models to {out_prefix}.{{cds,prot,meta}}.*", file=sys.stderr)


if __name__ == "__main__":
    main()
