#!/usr/bin/env bash
# run_batch.sh - run mine_gene.sh for every accession listed in an NCBI-style
# assembly_summary.txt (as produced by genome_updater.sh, or any file whose
# first tab-separated column is a GCA/GCF accession).
#
# Usage:
#   run_batch.sh --accession_list assembly_summary.txt \
#                --reference_proteins REF.prot --reference_alignment REF.aln \
#                --reference_tree REF.aln.treefile --gene_name GENE \
#                -o OUTDIR [-t THREADS] [extra mine_gene.sh options...]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ACCESSION_LIST=""
GENE_NAME=""
OUTDIR=""
THREADS=8
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --accession_list) ACCESSION_LIST="$2"; shift 2 ;;
        --gene_name) GENE_NAME="$2"; EXTRA_ARGS+=("--gene_name" "$2"); shift 2 ;;
        -o|--outdir) OUTDIR="$2"; EXTRA_ARGS+=("-o" "$2"); shift 2 ;;
        -t|--threads) THREADS="$2"; EXTRA_ARGS+=("-t" "$2"); shift 2 ;;
        *) EXTRA_ARGS+=("$1"); shift 1 ;;
    esac
done

if [[ -z "$ACCESSION_LIST" || -z "$GENE_NAME" || -z "$OUTDIR" ]]; then
    echo "ERROR: --accession_list, --gene_name and -o/--outdir are required" >&2
    echo "Usage: run_batch.sh --accession_list assembly_summary.txt --reference_proteins REF.prot \\" >&2
    echo "                    --reference_alignment REF.aln --reference_tree REF.aln.treefile \\" >&2
    echo "                    --gene_name GENE -o OUTDIR [-t THREADS] [extra mine_gene.sh options...]" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r accession _rest; do
    [[ -z "$accession" || "$accession" == "#"* ]] && continue
    summary_glob=("${OUTDIR}/${GENE_NAME}/"*".${accession}/"*".${GENE_NAME}.summary.tsv")
    if [[ -s "${summary_glob[0]}" ]]; then
        echo "[batch] $accession already done, skipping" >&2
        continue
    fi
    echo "[batch] === $accession ===" >&2
    bash "${SCRIPT_DIR}/mine_gene.sh" --ncbi_accession "$accession" "${EXTRA_ARGS[@]}" \
        2>&1 | sed "s/^/[${accession}] /"
done < <(awk -F'\t' '{print $1}' "$ACCESSION_LIST")

echo "[batch] All done. Combining summaries..." >&2
OUT_COMBINED="${OUTDIR}/all_species_${GENE_NAME}_summary.tsv"
first=1
for f in "$OUTDIR/$GENE_NAME"/*/*."$GENE_NAME".summary.tsv; do
    [[ -s "$f" ]] || continue
    if [[ "$first" -eq 1 ]]; then
        cat "$f" > "$OUT_COMBINED"
        first=0
    else
        tail -n +2 "$f" >> "$OUT_COMBINED"
    fi
done
echo "[batch] Combined summary: $OUT_COMBINED" >&2
