#!/usr/bin/env bash
# download_genome.sh - resolve an NCBI genome accession (GCA or GCF) to its FTP
# directory, download the genome FASTA and the organism name, with no
# dependency on a pre-fetched assembly_summary.txt.
#
# Usage: download_genome.sh <accession> <outdir>
#
# On success prints two lines to stdout:
#   line 1: resolved species name (spaces -> underscores, parenthetical
#           common name / strain info stripped)
#   line 2: path to the downloaded genome FASTA
set -euo pipefail

ACCESSION="$1"
OUTDIR="$2"

mkdir -p "$OUTDIR"
GENOME_FA="${OUTDIR}/${ACCESSION}.fa"
REPORT="${OUTDIR}/${ACCESSION}.assembly_report.txt"

log() { echo "[download_genome] $*" >&2; }

PREFIX="${ACCESSION%%_*}"                      # GCA or GCF
REST="${ACCESSION#*_}"                         # 013397245.1
DIGITS="${REST%%.*}"                           # 013397245
if [[ "$PREFIX" != "GCA" && "$PREFIX" != "GCF" ]]; then
    echo "[download_genome] ERROR: accession must start with GCA_ or GCF_, got: $ACCESSION" >&2
    exit 1
fi
P1="${DIGITS:0:3}"
P2="${DIGITS:3:3}"
P3="${DIGITS:6:3}"
PARENT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${PREFIX}/${P1}/${P2}/${P3}/"

if [[ ! -s "$REPORT" ]]; then
    log "Resolving FTP directory for $ACCESSION from $PARENT_URL"
    LISTING=$(curl -sL --retry 3 --max-time 60 "$PARENT_URL")
    ASM_DIR=$(echo "$LISTING" | grep -oP "href=\"\K${ACCESSION}_[^\"/]+(?=/)" | head -1)
    if [[ -z "$ASM_DIR" ]]; then
        echo "[download_genome] ERROR: could not find a directory for $ACCESSION under $PARENT_URL" >&2
        exit 1
    fi
    BASE_URL="${PARENT_URL}${ASM_DIR}/"
    log "Resolved assembly directory: $BASE_URL"

    curl -sL --retry 3 --max-time 60 "${BASE_URL}${ASM_DIR}_assembly_report.txt" -o "$REPORT"
    echo "$BASE_URL" > "${OUTDIR}/${ACCESSION}.ftp_dir.txt"
else
    log "Reusing cached assembly report: $REPORT"
    BASE_URL=$(cat "${OUTDIR}/${ACCESSION}.ftp_dir.txt")
    ASM_DIR=$(basename "$BASE_URL")
fi

SPECIES=$(grep -m1 "^# Organism name:" "$REPORT" | sed -E 's/^# Organism name:\s*//; s/\s*\([^)]*\)\s*$//; s/ /_/g')
if [[ -z "$SPECIES" ]]; then
    echo "[download_genome] ERROR: could not parse organism name from $REPORT" >&2
    exit 1
fi

if [[ -s "$GENOME_FA" ]]; then
    log "$GENOME_FA already exists, skipping download"
else
    GENOME_URL="${BASE_URL}${ASM_DIR}_genomic.fna.gz"
    GENOME_GZ="${OUTDIR}/${ACCESSION}.fna.gz"
    log "Downloading $GENOME_URL"
    curl -sL --retry 3 --max-time 1800 "$GENOME_URL" -o "$GENOME_GZ"
    if [[ ! -s "$GENOME_GZ" ]]; then
        echo "[download_genome] ERROR: download failed or empty file for $ACCESSION" >&2
        exit 1
    fi
    gunzip -c "$GENOME_GZ" > "$GENOME_FA"
    rm -f "$GENOME_GZ"
    log "Done: $GENOME_FA"
fi

echo "$SPECIES"
echo "$GENOME_FA"
