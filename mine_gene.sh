#!/usr/bin/env bash
# mine_gene.sh - mine a single-copy (or small-paralogue-family) gene from a
# bird genome assembly using miniprot, then place candidates on a reference
# tree to confirm they are nested in the target gene's clade (vs. an
# outgroup/suspect placement).
#
# Generalized from the original SI/ADAG-specific pipeline: point it at any
# gene by supplying your own reference protein set + alignment + tree.
#
# Usage:
#   mine_gene.sh --reference_proteins REF.prot --reference_alignment REF.aln \
#                --reference_tree REF.aln.treefile --gene_name GENE \
#                --ncbi_accession GCA_xxxxxxxxx.x -o OUTDIR [OPTIONS]
#
#   mine_gene.sh --reference_proteins REF.prot --reference_alignment REF.aln \
#                --reference_tree REF.aln.treefile --gene_name GENE \
#                -g GENOME.FA -s SPECIES -o OUTDIR [OPTIONS]
#
# Required:
#   --reference_proteins FILE   unaligned FASTA of the target gene's known
#                                homologues; used both as the miniprot query
#                                database and to define the "ingroup" leaf
#                                set for tree placement (any extra sequence in
#                                reference_alignment/reference_tree that is
#                                NOT in this file is treated as outgroup)
#   --reference_alignment FILE  protein alignment of reference_proteins PLUS
#                                any outgroup sequences (mafft --add target)
#   --reference_tree FILE       FastTree newick tree built from
#                                reference_alignment (used only to confirm
#                                gene_name/labeling context; placement itself
#                                rebuilds a tree each run after adding
#                                candidates, see step 6)
#   --gene_name NAME             short label used in output file/dir names
#                                and as the ingroup clade label
#   one of:
#     --ncbi_accession ACC       GCA/GCF accession; genome + species name are
#                                resolved and downloaded automatically
#     -g/--genome FILE           pre-downloaded genome FASTA (then -s/--species
#                                is required unless --ncbi_accession is also
#                                given just to resolve the species name)
#   -o/--outdir DIR
#
# Optional:
#   -s/--species NAME            override the auto-resolved species name
#   -t/--threads INT              default 4
#   --min-identity FLOAT          default 0.25
#   --maxintron INT                default 200000
#   --merge-dist INT               default 5000
#   --complete-frac FLOAT          default 0.8
#   --clean-tmp                    delete tmp/ and the genome index when done
#   --skip-permissive-pass         run only the fast, default-settings miniprot
#                                   pass. Faster, but will miss a secondary
#                                   paralogous locus whenever a stronger-
#                                   scoring copy of the gene exists elsewhere
#                                   in the same genome (see README, "Known
#                                   limitations")
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v miniprot >/dev/null 2>&1; then
    echo "ERROR: miniprot not found on PATH. Did you 'conda activate gene_miner_env'?" >&2
    exit 1
fi

REFERENCE_PROTEINS=""
REFERENCE_ALIGNMENT=""
REFERENCE_TREE=""
GENE_NAME=""
ACCESSION=""
GENOME_FA=""
SPECIES=""
OUTDIR=""
THREADS=4
MIN_IDENTITY=0.25
MAX_INTRON=200000
MERGE_DIST=5000
COMPLETE_FRAC=0.8
CLEAN_TMP=0
SKIP_PERMISSIVE_PASS=0

usage() { grep '^# ' "$0" | sed 's/^# \?//'; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference_proteins) REFERENCE_PROTEINS="$2"; shift 2 ;;
        --reference_alignment) REFERENCE_ALIGNMENT="$2"; shift 2 ;;
        --reference_tree) REFERENCE_TREE="$2"; shift 2 ;;
        --gene_name) GENE_NAME="$2"; shift 2 ;;
        --ncbi_accession|-a) ACCESSION="$2"; shift 2 ;;
        -g|--genome) GENOME_FA="$2"; shift 2 ;;
        -s|--species) SPECIES="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --min-identity) MIN_IDENTITY="$2"; shift 2 ;;
        --maxintron) MAX_INTRON="$2"; shift 2 ;;
        --merge-dist) MERGE_DIST="$2"; shift 2 ;;
        --complete-frac) COMPLETE_FRAC="$2"; shift 2 ;;
        --clean-tmp) CLEAN_TMP=1; shift 1 ;;
        --skip-permissive-pass) SKIP_PERMISSIVE_PASS=1; shift 1 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

for req in REFERENCE_PROTEINS REFERENCE_ALIGNMENT REFERENCE_TREE GENE_NAME OUTDIR; do
    if [[ -z "${!req}" ]]; then
        echo "ERROR: --${req,,} is required" >&2
        usage
    fi
done
if [[ -z "$GENOME_FA" && -z "$ACCESSION" ]]; then
    echo "ERROR: provide either -g GENOME.FA or --ncbi_accession ACC" >&2
    usage
fi
if [[ -n "$GENOME_FA" && -z "$ACCESSION" && -z "$SPECIES" ]]; then
    echo "ERROR: -s/--species is required when using -g without --ncbi_accession" >&2
    usage
fi

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

############################################
# 0. sanity-check that the three reference files agree with each other
############################################
log "Validating reference_proteins / reference_alignment / reference_tree"
python3 "${SCRIPT_DIR}/scripts/validate_references.py" \
    "$REFERENCE_PROTEINS" "$REFERENCE_ALIGNMENT" "$REFERENCE_TREE"

############################################
# 1. genome + species resolution (cached per-accession, shared across genes)
############################################
GENOME_CACHE="${OUTDIR}/genomes_cache"
mkdir -p "$GENOME_CACHE"

if [[ -n "$ACCESSION" ]]; then
    RESOLVED=$(bash "${SCRIPT_DIR}/scripts/download_genome.sh" "$ACCESSION" "$GENOME_CACHE")
    RESOLVED_SPECIES=$(echo "$RESOLVED" | sed -n '1p')
    RESOLVED_GENOME=$(echo "$RESOLVED" | sed -n '2p')
    [[ -z "$SPECIES" ]] && SPECIES="$RESOLVED_SPECIES"
    [[ -z "$GENOME_FA" ]] && GENOME_FA="$RESOLVED_GENOME"
fi

ACCESSION_LABEL="${ACCESSION:-$(basename "${GENOME_FA%.*}")}"
LABEL="${SPECIES}.${ACCESSION_LABEL}"
SPDIR="${OUTDIR}/${GENE_NAME}/${LABEL}"
TMPDIR="${SPDIR}/tmp"
mkdir -p "$SPDIR" "$TMPDIR"

############################################
# 2. miniprot index (cached per-accession, shared across genes)
############################################
if [[ -n "$ACCESSION" ]]; then
    MPI="${GENOME_CACHE}/${ACCESSION}.mpi"
else
    MPI="${SPDIR}/genome.mpi"
fi
if [[ ! -s "$MPI" ]]; then
    log "Building miniprot index"
    miniprot -t "$THREADS" -d "$MPI" "$GENOME_FA" >/dev/null
else
    log "miniprot index already present: $MPI"
fi

############################################
# 3. miniprot alignment (reference_proteins as queries), two passes
############################################
# Pass 1 (default miniprot settings) reports, for each query, only its single
# best hit genome-wide. If a genome carries more than one paralogous copy of
# the gene at different divergence levels (e.g. a lineage-specific
# duplication alongside a degraded/truncated copy), every query's best hit
# lands on the strongest copy and the weaker one is never reported at all --
# not filtered downstream, simply never emitted by miniprot. Pass 2 relaxes
# both the chaining secondary/primary ratio (-p) and the output score
# fraction (--outs) to recover those secondary loci; dedup_loci.py (step 5)
# then collapses the redundant hits from both passes down to one best model
# per locus, so merging them here is safe.
RAW_GFF="${TMPDIR}/raw.gff"
if [[ ! -s "$RAW_GFF" ]]; then
    log "Running miniprot pass 1/2 (default settings, max-intron=${MAX_INTRON})"
    miniprot -t "$THREADS" --gff -G "$MAX_INTRON" \
        "$MPI" "$REFERENCE_PROTEINS" > "${TMPDIR}/raw_pass1.gff" 2> "${TMPDIR}/miniprot_pass1.log"

    if [[ "$SKIP_PERMISSIVE_PASS" -eq 0 ]]; then
        log "Running miniprot pass 2/2 (permissive: recovers secondary paralogous loci)"
        miniprot -t "$THREADS" --gff -G "$MAX_INTRON" -p 0.1 --outs=0.1 -N 50 -P P2 \
            "$MPI" "$REFERENCE_PROTEINS" > "${TMPDIR}/raw_pass2.gff" 2> "${TMPDIR}/miniprot_pass2.log"
        cat "${TMPDIR}/raw_pass1.gff" "${TMPDIR}/raw_pass2.gff" > "$RAW_GFF"
    else
        cp "${TMPDIR}/raw_pass1.gff" "$RAW_GFF"
    fi
else
    log "Raw GFF already present: $RAW_GFF"
fi

############################################
# 4. parse GFF -> CDS / protein / metadata
############################################
log "Parsing GFF"
python3 "${SCRIPT_DIR}/scripts/parse_miniprot_gff.py" \
    "$RAW_GFF" "$GENOME_FA" "${TMPDIR}/raw"

############################################
# 5. dedup: one best model per locus
############################################
log "Deduplicating overlapping gene models"
python3 "${SCRIPT_DIR}/scripts/dedup_loci.py" \
    "${TMPDIR}/raw.meta.tsv" "${TMPDIR}/raw.cds.fa" "${TMPDIR}/raw.prot.fa" \
    "${TMPDIR}/dedup" --min-identity "$MIN_IDENTITY" --merge-dist "$MERGE_DIST"

N_LOCI=$(($(wc -l < "${TMPDIR}/dedup.meta.tsv") - 1))
log "$N_LOCI candidate locus/loci found"

if [[ "$N_LOCI" -le 0 ]]; then
    log "No candidate loci found for $SPECIES ($ACCESSION) -- writing empty outputs"
    : > "${SPDIR}/${LABEL}.${GENE_NAME}.summary.tsv"
    : > "${SPDIR}/${LABEL}.${GENE_NAME}.cds.fa"
    : > "${SPDIR}/${LABEL}.${GENE_NAME}.prot.fa"
    exit 0
fi

############################################
# 6. phylogenetic placement: mafft --add --keeplength, then FastTree
############################################
log "Aligning candidates to the reference alignment (mafft --add --keeplength)"
mafft --thread "$THREADS" --add "${TMPDIR}/dedup.prot.fa" --keeplength \
    "$REFERENCE_ALIGNMENT" > "${TMPDIR}/placed.aln" 2> "${TMPDIR}/mafft.log"

log "Building FastTree"
FastTree < "${TMPDIR}/placed.aln" > "${TMPDIR}/placed.aln.treefile" 2> "${TMPDIR}/fasttree.log"

grep '^>' "${TMPDIR}/dedup.prot.fa" | sed 's/^>//' > "${TMPDIR}/candidate_ids.txt"

log "Classifying candidates on the tree (ingroup/outgroup call + long-branch suspects)"
python3 "${SCRIPT_DIR}/scripts/classify_candidates_tree.py" \
    "${TMPDIR}/placed.aln.treefile" "${TMPDIR}/candidate_ids.txt" \
    "$REFERENCE_PROTEINS" "$GENE_NAME" "${TMPDIR}/tree_classification.tsv"

############################################
# 7. final report + output FASTA
############################################
log "Writing final report"
python3 "${SCRIPT_DIR}/scripts/finalize_report.py" \
    "${TMPDIR}/dedup.meta.tsv" "${TMPDIR}/dedup.cds.fa" "${TMPDIR}/dedup.prot.fa" \
    "${TMPDIR}/tree_classification.tsv" "$REFERENCE_PROTEINS" "$GENE_NAME" \
    "$SPECIES" "$ACCESSION_LABEL" "${SPDIR}/${LABEL}.${GENE_NAME}" --complete-frac "$COMPLETE_FRAC"

if [[ "$CLEAN_TMP" -eq 1 ]]; then
    log "Cleaning up tmp/ (genome index is kept in genomes_cache/ for reuse across genes)"
    rm -rf "$TMPDIR"
fi

log "Done. Results in ${SPDIR}/"
