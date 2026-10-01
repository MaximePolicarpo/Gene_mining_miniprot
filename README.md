# Gene_mining_miniprot

A lightweight pipeline for mining a single-copy gene (or a small family of
close paralogues, like a gene and its rare in-group duplicate) out of raw
genome assemblies, built on **miniprot**, **mafft** and **FastTree**.

> **Everything this pipeline outputs is a candidate, not a finished
> annotation.** It tells you where a plausible gene model is and whether it
> places inside your target gene's clade or looks suspect — it does not
> replace looking at the alignment and the tree yourself before trusting a
> sequence, especially anything the pipeline flagged `AMBIGUOUS`,
> `OUTGROUP_SUSPECT`, or `long_branch_suspect=True`. See
> [§9](#9-reading-the-output-candidates-not-final-annotations) before using
> any sequence downstream.

---

## 1. Installation

```bash
git clone <this-repo-url>
cd Gene_mining_miniprot
conda env create -f environment.yml
conda activate gene_miner_env
```

This installs Python 3.10, Biopython, pandas, miniprot, samtools, mafft,
FastTree, seqkit and curl — everything the pipeline needs, nothing else.

---

## 2. What you need to provide

Three files describe the gene you want to mine, plus a gene name:

| Flag | What it is |
|---|---|
| `--reference_proteins FILE` | Unaligned FASTA of the target gene's known homologues (e.g. one ortholog per species, across the clade you care about). Used as the **miniprot query database**, and defines the **"ingroup"** leaf set for tree placement. |
| `--reference_alignment FILE` | A protein alignment containing `reference_proteins` **plus** a handful of outgroup sequences (e.g. a related paralogue family) used purely to root the tree. Any sequence in this file that is *not* in `reference_proteins` is automatically treated as outgroup — no header-naming convention required. |
| `--reference_tree FILE` | The FastTree newick tree built from `reference_alignment`. The pipeline doesn't actually place candidates onto this fixed tree (see step 6 below) — it exists so the pipeline can **sanity-check that your three reference files agree with each other** before doing any real work (same leaf set in the alignment and tree, `reference_proteins` fully contained in the alignment). |
| `--gene_name NAME` | A short label used in output file/directory names and as the ingroup clade label in reports. |

See [§8](#8-building-your-own-reference-set-for-a-new-gene) for how the
bundled SI/ADAG example was built — the same recipe generalizes to any gene.

---

## 3. Quick start

### Mine one genome by NCBI accession

```bash
bash mine_gene.sh \
    --reference_proteins example_data/SI_ADAG/SI_ADAG_references.prot \
    --reference_alignment example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln \
    --reference_tree example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln.treefile \
    --gene_name SI_ADAG \
    --ncbi_accession GCA_982554685.1 \
    -o results/ \
    -t 8
```

This downloads the genome straight from NCBI (resolving the FTP path and
organism name from the accession alone — no separate species flag needed),
mines it, and writes results to
`results/SI_ADAG/Certhia_brachydactyla.GCA_982554685.1/`.

### Mine a pre-downloaded genome

```bash
bash mine_gene.sh \
    --reference_proteins example_data/SI_ADAG/SI_ADAG_references.prot \
    --reference_alignment example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln \
    --reference_tree example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln.treefile \
    --gene_name SI_ADAG \
    -g /path/to/my_genome.fa -s My_species \
    -o results/ -t 8
```

### Batch mode: mine every accession in a list

```bash
bash run_batch.sh \
    --accession_list assembly_summary.txt \
    --reference_proteins example_data/SI_ADAG/SI_ADAG_references.prot \
    --reference_alignment example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln \
    --reference_tree example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln.treefile \
    --gene_name SI_ADAG \
    -o results/ -t 8
```

`assembly_summary.txt` is any tab-separated file whose first column is a
GCA/GCF accession (exactly the format NCBI's `assembly_summary.txt` or
[`genome_updater`](https://github.com/pirovc/genome_updater) produce — other
columns are ignored, species names are resolved automatically). The batch
runner skips any accession it has already finished, so an interrupted run
can simply be relaunched, and writes a combined
`results/all_species_<gene_name>_summary.tsv` at the end.

---

## 4. Running on a cluster with no internet on compute/login nodes (e.g. GWDG)

Many HPC clusters only allow outbound internet access (needed to download a
genome, or to let `mine_gene.sh` auto-resolve `--ncbi_accession`) from the
**login node**, not from compute nodes reached via `srun`/`sbatch`. On these
clusters you cannot use `--ncbi_accession` directly in a compute job — download
the genome by hand on the login node first, then point `mine_gene.sh` at the
local file with `-g` (which never needs internet).

### Step 1 — on the login node: download the genome

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/791/375/GCF_027791375.1_UM_Iind_1.1/GCF_027791375.1_UM_Iind_1.1_genomic.fna.gz
gzip -d GCF_027791375.1_UM_Iind_1.1_genomic.fna.gz
mkdir -p results/genomes_cache
mv GCF_027791375.1_UM_Iind_1.1_genomic.fna results/genomes_cache/
```

(Build the URL from the accession the same way NCBI's own FTP layout does:
`genomes/all/<GCA|GCF>/<first 3 digits>/<next 3>/<next 3>/<accession>_<assembly name>/`.
The assembly name is in the directory listing, or in `*_assembly_report.txt`.)

### Step 2 — request a compute node

Interactively:

```bash
srun \
  --account=<your_slurm_account> \
  --partition=<a_cpu_partition> \
  --nodes=1 \
  --cpus-per-task=8 \
  --mem=50G \
  --pty bash
```

(`--mem=50G` is comfortable headroom for indexing a ~1 GB bird genome with
both miniprot passes running — see §6 step 2's memory note. Scale up for
larger genomes.)

### Step 3 — on the compute node: run the pipeline against the local genome

```bash
conda activate gene_miner_env

bash mine_gene.sh \
    --reference_proteins example_data/SI_ADAG/SI_ADAG_references.prot \
    --reference_alignment example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln \
    --reference_tree example_data/SI_ADAG/SI_ADAG.templatealignment.outgroup.prot.aln.treefile \
    --gene_name SI_ADAG \
    -g results/genomes_cache/GCF_027791375.1_UM_Iind_1.1_genomic.fna -s Indicator_indicator \
    -o results/ -t 8
```

`-s`/`--species` is required here since there's no accession for the
pipeline to resolve the organism name from.

**Alternative: submit as a batch job instead of an interactive `srun`.**
Wrap the same `mine_gene.sh` call (step 3) in a job script and `sbatch` it —
just don't forget to `conda activate gene_miner_env` inside the script
itself, before the `bash mine_gene.sh ...` line, since a batch job doesn't
inherit your interactive shell's environment.

### Step 4 — clean up the genome index afterward

The miniprot index (`genome.mpi`, several GB) is kept after the run so a
second gene can reuse it without rebuilding. Once you're done mining every
gene you need from this genome, delete it to reclaim disk space:

```bash
# when using -g without --ncbi_accession, the index lives per-species:
rm results/SI_ADAG/Indicator_indicator.GCF_027791375.1_UM_Iind_1.1_genomic/genome.mpi

# when using --ncbi_accession, it's shared in the cache instead:
rm results/genomes_cache/<ACCESSION>.mpi
```

---

## 5. Full option reference

```
mine_gene.sh --reference_proteins FILE --reference_alignment FILE \
             --reference_tree FILE --gene_name NAME \
             (--ncbi_accession ACC | -g GENOME.FA -s SPECIES) \
             -o OUTDIR [OPTIONS]

  -s, --species NAME        override the auto-resolved species name
  -t, --threads INT         CPU threads for miniprot/mafft [4]
  --min-identity FLOAT      discard miniprot hits below this identity
                            before clustering into loci [0.25]
  --maxintron INT           max intron size passed to miniprot (-G) [200000]
  --merge-dist INT          merge gene models within this many bp into
                            one locus [5000]
  --complete-frac FLOAT     a locus is "functional_complete" if its protein
                            is >= this fraction of the median reference
                            protein length [0.8]
  --clean-tmp               delete tmp/ (GFF, raw candidates, alignment,
                            tree) after a successful run; the genome index
                            is always kept in genomes_cache/ for reuse
```

---

## 6. Pipeline steps

```
  reference_proteins / reference_alignment / reference_tree
                      │
                      ▼
  0. validate_references.py
     Sanity-check the three reference files are mutually consistent
     (reference_proteins ⊆ reference_alignment leaf set == reference_tree
     leaf set) before spending any compute.
                      │
  NCBI accession ─────┤                    pre-downloaded genome ──────┐
       │              │                                                │
       ▼              │                                                │
  1. download_genome.sh                                                │
     Resolve the FTP directory straight from the accession (no local   │
     assembly_summary.txt lookup needed), download + gunzip the        │
     genome FASTA, and parse the organism name from the assembly       │
     report. Cached under outdir/genomes_cache/ by accession, shared   │
     across every gene you later mine against the same genome.         │
       └──────────────┬─────────────────────────────────────┬─────────┘
                       ▼
  2. miniprot --index
     Build the genome's miniprot index once (also cached in
     genomes_cache/ by accession). Needs roughly 8-10x the genome's
     uncompressed size in RAM (a ~1.1 GB bird genome peaked at ~10 GB RSS
     in testing) -- request memory accordingly on a cluster (see §4).
                       │
                       ▼
  3. miniprot (protein-to-genome spliced alignment), two passes
     Pass 1 (default settings) aligns every sequence in reference_proteins
     against the genome. Because the reference set is usually many closely
     related homologues (one per species), most of them will hit the *same*
     true locus — that redundancy is resolved in step 5, not here. But
     miniprot's default mode reports, for each query, only its single best
     hit genome-wide (chains scoring < 70% of that query's top chain are
     discarded before alignment, and only near-top alignments are output).
     If a genome carries more than one paralogous copy of the gene at
     different divergence levels (e.g. a lineage-specific duplication
     alongside an older, more degraded copy), every query's best hit lands
     on the strongest copy and a real but weaker-scoring copy is never
     reported at all — not filtered downstream, simply never emitted.
     Pass 2 relaxes the chaining ratio and output threshold
     (-p 0.1 --outs 0.1 -N 50) specifically to recover those secondary
     loci; both passes' hits are pooled before step 5, which collapses the
     redundancy between them. Skip pass 2 with --skip-permissive-pass if
     you know your gene is strictly single-copy and want the faster run —
     expect roughly 2-5x longer runtime with it on, and more candidate
     loci to review (see §10).
                       │
                       ▼
  4. parse_miniprot_gff.py
     Splice each reported gene model's CDS out of the genome, translate
     it, and record miniprot's own Identity / Frameshift / StopCodon
     tags plus its raw alignment score.
                       │
                       ▼
  5. dedup_loci.py
     Drop hits below --min-identity, merge overlapping/nearby (within
     --merge-dist) gene models into one locus, and keep the single
     best-scoring model per locus — ranked by miniprot's native
     alignment score (not raw identity: a self-hit might be 99%
     identical but truncated, while a slightly-lower-identity hit from
     a close relative can cover the full gene; the native score already
     balances identity against coverage, so it picks the right one).
     Every surviving locus goes on to the next step; nothing is
     filtered out based on gene identity beyond this point.
                       │
                       ▼
  6. mafft --add --keeplength  +  FastTree
     Add every candidate protein to reference_alignment (without
     disturbing its existing columns) and build a fresh tree including
     the candidates. A new tree is built every run rather than placing
     onto the fixed reference_tree, since --add --keeplength only needs
     an alignment, and rebuilding is cheap at this scale (a few hundred
     reference leaves, a handful of candidates).
                       │
                       ▼
  7. classify_candidates_tree.py
     For each candidate, compute its tree (patristic) distance to the
     nearest ingroup references (reference_proteins) and the nearest
     outgroup references (whatever is in reference_alignment but not in
     reference_proteins). Call it "<gene_name>" if it's closer to the
     ingroup, "OUTGROUP_SUSPECT" if closer to the outgroup, or
     "AMBIGUOUS(...)" if the two are too close to call confidently.
     Also flag "long_branch_suspect" if its terminal branch is more than
     5x the median ingroup terminal branch length (a sign of a bad
     alignment, a contaminant, or a badly degraded sequence).
                       │
                       ▼
  8. finalize_report.py
     Merge the locus metadata with the tree call into one per-species
     report, call each locus:
       pseudogene            any internal stop / frameshift / in-frame
                              stop reported by miniprot
       functional_complete   LoF-free and protein length >= complete-frac
                              x median reference_proteins length
       functional_partial    LoF-free but shorter (usually a scaffold-
                              edge or assembly-gap truncation)
     and writes the final <species>.<accession>.<gene_name>.{summary.tsv,
     cds.fa,prot.fa}.
```

Nothing is silently discarded after step 5's identity floor — every
surviving locus is reported with its functional call, clade call, and
suspect flags, so you can review borderline cases yourself rather than have
them disappear.

---

## 7. Output layout

```
OUTDIR/
├── genomes_cache/                       # shared across every gene_name run
│   ├── GCA_xxxxxxxxx.x.fa               # downloaded genome
│   ├── GCA_xxxxxxxxx.x.mpi              # miniprot index (reused across genes)
│   └── GCA_xxxxxxxxx.x.assembly_report.txt
│
├── <GENE_NAME>/
│   └── <Species_name>.<accession>/
│       ├── tmp/                         # raw GFF, candidates, alignment, tree
│       │   └── ...                      # (removed by --clean-tmp)
│       ├── <label>.<gene_name>.summary.tsv
│       ├── <label>.<gene_name>.cds.fa            # confident ingroup calls only
│       ├── <label>.<gene_name>.prot.fa           # confident ingroup calls only
│       ├── <label>.<gene_name>.ambiguous.cds.fa  # AMBIGUOUS/OUTGROUP_SUSPECT calls
│       └── <label>.<gene_name>.ambiguous.prot.fa # AMBIGUOUS/OUTGROUP_SUSPECT calls
│
└── all_species_<gene_name>_summary.tsv  # written by run_batch.sh
```

`summary.tsv` columns: `gene_id, mrna_id, species, accession, chrom, start,
end, strand, best_query, identity, n_exons, cds_len_nt, prot_len_aa,
internal_stops, gff_frameshift, gff_stopcodon, functional_call, clade_call,
terminal_branch_length, long_branch_suspect, dist_to_ingroup,
dist_to_outgroup`. It lists every locus found, confident or not — only the
FASTA files are split (see below).

**Why two sets of FASTA files:** `clade_call` is `<gene_name>` for a
confident ingroup placement, or `AMBIGUOUS(...)` / `OUTGROUP_SUSPECT` when
the tree placement step couldn't confirm it (see §9). Those aren't
lower-confidence versions of a real gene — they're sequences the pipeline
explicitly failed to confirm, usually a distant unrelated homologue picked
up by the permissive second miniprot pass (§6 step 3). Mixing them into the
main `cds.fa`/`prot.fa` would make every downstream use of those files (an
alignment, a tree, a BLAST database) need re-filtering first, so they're
routed to `*.ambiguous.cds.fa`/`*.ambiguous.prot.fa` instead. Nothing is
deleted — they're still in `summary.tsv` and in their own FASTA pair — just
kept out of the files you'd use directly.

---

## 8. Building your own reference set for a new gene

This is the exact recipe used to build the bundled SI/ADAG example, with the
actual commands — swap in your gene, your clade, and your own sequence
sources. It assumes `samtools`, `emboss` (for `transeq`), `muscle`, `trimal`
and `mafft` are available (`muscle`/`trimal` aren't in `environment.yml`
since they're only needed once, to build a reference set, not at mining
time — `conda install -c bioconda muscle trimal` or similar).

### 8.1 Collect candidate reference CDS sequences

Gather CDS sequences for your gene across the clade you're mining, plus any
close in-group paralogue you also want to detect (for SI/ADAG, both genes
went through this pipeline together). RefSeq annotations are the easiest
source — for birds, that meant pulling the RefSeq `rna.fna`/CDS sequences
annotated for this gene across every available RefSeq bird genome (NCBI's
[Orthologs](https://www.ncbi.nlm.nih.gov/gene/) view for a model species'
gene, or a bulk per-gene CDS download across an Entrez/Datasets query, both
work) into one multi-FASTA, one sequence per species:

```bash
cat *_rna.fna > all_candidates.cds
samtools faidx all_candidates.cds
```

### 8.2 Filter out truncated annotations

RefSeq/automated annotations for a gene this large (SI spans ~45 exons) are
frequently truncated in draft assemblies. Pick a well-annotated reference
species for your clade (chicken for birds) and keep only candidates that are
at least 80% of its CDS length:

```bash
# extract the reference species' own CDS and get its length
samtools faidx all_candidates.cds Gallus_gallus_SI_transcript_id > ref.cds
samtools faidx ref.cds
REF_LEN=$(cut -f2 ref.cds.fai)
THRESH=$(awk -v l="$REF_LEN" 'BEGIN{print l*0.8}')

# list and extract every candidate at or above that length
awk -v t="$THRESH" '$2>=t {print $1}' all_candidates.cds.fai > pass_ids.txt
samtools faidx all_candidates.cds -r pass_ids.txt > filtered.cds
```

If you're mining an in-group paralogue alongside the main gene (like ADAG
alongside SI) and already have a rough phylogenetic split of which candidate
belongs to which, repeat this filter separately per gene (each against its
own reference length), then keep the two filtered sets separate for now —
you'll label them when building the alignment in §8.3.

### 8.3 Translate to protein

```bash
transeq filtered.cds filtered.prot
sed -i 's/_1$//' filtered.prot   # transeq appends _1 to every header; strip it
```

If you have two gene sets (e.g. SI and ADAG), prefix each header so they
stay identifiable once merged — `sed 's/^>/>SI---/' SI_filtered.prot`,
`sed 's/^>/>ADAG---/' ADAG_filtered.prot` — then `cat` them together into
one `reference_proteins.fa`. This labeling is just for your own bookkeeping
at this stage; the pipeline itself never looks at header prefixes (see §2).

### 8.4 Build a trimmed template alignment, then add everything else

Aligning hundreds of full-length sequences directly tends to produce messy,
gappy alignments and wastes time re-aligning near-identical sequences from
closely related species. Instead, align and trim a small representative
subset first, then add the rest onto that fixed, trimmed column set:

```bash
# pick ~15 representative sequences per gene, plus your reference species
grep ">" reference_proteins.fa | grep "^>SI---"   | shuf -n15 | sed 's/>//' > subset.id
grep ">" reference_proteins.fa | grep "^>ADAG---" | shuf -n15 | sed 's/>//' >> subset.id
echo "Gallus_gallus_SI_id"   >> subset.id   # always keep your reference species
echo "Gallus_gallus_ADAG_id" >> subset.id
sort -u subset.id -o subset.id

xargs samtools faidx reference_proteins.fa < subset.id > subset.prot
muscle -align subset.prot -output subset.prot.aln
trimal -in subset.prot.aln -out subset.trimmed.prot.aln -gt 0.7   # drop columns >30% gaps

# add every remaining sequence onto the trimmed alignment without
# disturbing its columns
grep ">" reference_proteins.fa | sed 's/>//' | sort -u > all.id
comm -23 all.id subset.id > remaining.id
xargs samtools faidx reference_proteins.fa < remaining.id > remaining.prot
mafft --add remaining.prot --keeplength subset.trimmed.prot.aln > reference_alignment.fa
```

### 8.5 Add an outgroup

A handful of related-but-clearly-distinct sequences, added purely to root
the tree — without an outgroup, the tree placement step in §6 has no fixed
point to tell "nested in your gene's clade" from "nested somewhere else."
For SI/ADAG this meant the orthologous SI and MGAM (a related paralogue)
protein sequences from three mammals with good annotations — *Felis catus*,
*Mus musculus* and *Homo sapiens* — fetched from RefSeq, concatenated into
one small FASTA, and added the same way as §8.4:

```bash
mafft --add outgroup.prot --keeplength reference_alignment.fa > reference_alignment.outgroup.fa
```

These outgroup sequences must **not** be added to `reference_proteins.fa`
(§8.1-8.3) — the pipeline tells ingroup from outgroup purely by whether a
sequence is in `reference_proteins` or not (see §2), so putting them in both
files would make every outgroup sequence register as ingroup instead.

Three or four distantly-related species is normally enough; you don't need
exhaustive outgroup sampling; you need it reliably on the other side of the
root, further from your gene's clade than any of your real ingroup
sequences.

### 8.6 Build the tree

```bash
FastTree < reference_alignment.outgroup.fa > reference_alignment.outgroup.fa.treefile
```

You now have your three pipeline inputs: `reference_proteins.fa` (from §8.3,
pre-outgroup), `reference_alignment.outgroup.fa` (→ `--reference_alignment`)
and `reference_alignment.outgroup.fa.treefile` (→ `--reference_tree`). Run
`validate_references.py` (step 0 of §6, or just run `mine_gene.sh` — it's
called automatically) to confirm the three agree before using them for real.

---

## 9. Reading the output: candidates, not final annotations

Every sequence in `*.summary.tsv` / `*.cds.fa` / `*.prot.fa` is a **candidate
gene model**, not a confirmed, publication-ready gene. `functional_call`,
`clade_call` and the suspect flags are there to triage *which* candidates
need a closer look and roughly how urgently — they are not a substitute for
actually looking.

**Before trusting a sequence, at minimum:**

1. **Look at where it falls in the tree**, not just the summary table's
   `clade_call` column. The tree itself is the real evidence; the column is
   just our attempt to summarize it into one word. It's written to
   `tmp/placed.aln.treefile` (Newick — open it in
   [FigTree](https://github.com/rambaut/figtree/), `ete3`, iTOL, or any
   tree viewer) alongside `tmp/placed.aln`, the alignment it was built from.
   This is exactly why `--clean-tmp` deletes `tmp/` only after you're done
   inspecting a run — don't use it until you've checked the tree.
2. **`AMBIGUOUS(...)`/`OUTGROUP_SUSPECT` calls are already routed to
   `*.ambiguous.cds.fa`/`*.ambiguous.prot.fa` (see §7)** rather than the main
   FASTA files, most often because the permissive second miniprot pass
   surfaced a distant, unrelated homologue (§6 step 3) — they're kept for
   transparency, not because they're likely to be real. A confident
   ingroup call with `long_branch_suspect=True` stays in the main file but
   still deserves a look: it's a real placement, just an unusually divergent
   one.
3. **Treat `functional_partial` as "probably real but incomplete"**, not
   "broken." It usually means the locus continues past the end of whichever
   reference protein produced the best-scoring alignment there (see §10),
   not that the gene itself is truncated in the genome — worth extending
   the search region or checking the raw genomic sequence by hand.

**Manually curating a candidate into a finished gene model** (extending
exon boundaries, resolving an ambiguous splice site, correcting a probable
assembly error, etc.) is a separate step this pipeline doesn't attempt —
if you're unsure how to do this, ask Jinseok, since the manual curation
conventions used for this project aren't written down here (yet).

---

## 10. Known limitations

- **Tandem duplicates close together.** `dedup_loci.py` merges gene models
  within `--merge-dist` (default 5 kb) into a single locus and keeps only
  the best one. This is the right call for a gene that's normally
  single-copy (like SI), but if you point this pipeline at a gene family
  with real tandem duplicates sitting within that distance, the second copy
  would be silently collapsed away. The upstream
  [Chemoreceptor_mining_miniprot](https://github.com/MaximePolicarpo/Chemoreceptor_mining_miniprot)
  pipeline handles that case with dedicated tandem-rescue logic; it wasn't
  ported here since it's unnecessary complexity for a single/few-copy gene.
- **`--min-identity` is a noise floor, not a correctness filter.** Keep it
  low rather than raise it if you're worried about missing divergent
  pseudogenes — the tree placement step is what actually confirms or
  rejects a candidate, and it's cheap to run on a few extra loci.
- **The permissive second pass trades extra sensitivity for extra noise.**
  It's what makes a genuine secondary paralogue (see §6, step 3) detectable
  at all, but on a genome with no such paralogue it will often also surface
  a handful of distant, unrelated homologues (e.g. a different gene in the
  same broader enzyme family, or a processed pseudogene fragment) that
  wouldn't otherwise be reported. These are not false "complete" calls —
  they reliably come back `AMBIGUOUS(...)` or `OUTGROUP_SUSPECT` with
  `long_branch_suspect=True` from the tree placement step, so budget for a
  quick manual look at any locus flagged that way rather than treating every
  row in `summary.tsv` as confirmed. Use `--skip-permissive-pass` if you'd
  rather not deal with this and are confident the gene is single-copy in
  your species set.
- No stop-codon / start-codon extension is attempted beyond what miniprot's
  own spliced alignment reports; a `functional_partial` call often just
  means the gene continues past the end of whichever reference protein
  produced the best-scoring alignment at that locus, not necessarily that
  the genome assembly itself is incomplete there.
