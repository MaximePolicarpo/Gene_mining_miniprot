# Gene_mining_miniprot

A lightweight pipeline for mining a single-copy gene (or a small family of
close paralogues, like a gene and its rare in-group duplicate) out of raw
genome assemblies, built on **miniprot**, **mafft** and **FastTree**.

It was originally written to recover the sucrase-isomaltase gene (**SI**) and
its rare avian paralogue (**ADAG**) across bird genome assemblies, but the
pipeline itself is gene-agnostic: point it at any gene by supplying your own
reference protein set, alignment and tree. The SI/ADAG files are bundled
under `example_data/` as a worked example.

Design is loosely inspired by
[Chemoreceptor_mining_miniprot](https://github.com/MaximePolicarpo/Chemoreceptor_mining_miniprot)
(Policarpo et al.), a much larger pipeline built for multi-hundred-copy
chemoreceptor repertoires (OR, V1R, V2R, TAAR, T1R, T2R). This pipeline
reuses the same core idea — miniprot for protein-to-genome spliced alignment,
phylogenetic placement to confirm gene identity — simplified for genes that
exist in only one or two copies per genome, where you don't need tandem-array
rescue logic, dedicated multi-pass miniprot runs, or family-specific filters.

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

See [§7](#7-building-your-own-reference-set-for-a-new-gene) for how the
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

## 4. Full option reference

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

## 5. Pipeline steps

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
     genomes_cache/ by accession).
                       │
                       ▼
  3. miniprot (protein-to-genome spliced alignment)
     Align every sequence in reference_proteins against the genome.
     Because the reference set is usually many closely related homologues
     (one per species), most of them will hit the *same* true locus —
     that redundancy is resolved in step 5, not here.
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

## 6. Output layout

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
│       ├── <label>.<gene_name>.cds.fa
│       └── <label>.<gene_name>.prot.fa
│
└── all_species_<gene_name>_summary.tsv  # written by run_batch.sh
```

`summary.tsv` columns: `gene_id, mrna_id, species, accession, chrom, start,
end, strand, best_query, identity, n_exons, cds_len_nt, prot_len_aa,
internal_stops, gff_frameshift, gff_stopcodon, functional_call, clade_call,
terminal_branch_length, long_branch_suspect, dist_to_ingroup,
dist_to_outgroup`.

---

## 7. Building your own reference set for a new gene

This is how the bundled SI/ADAG example was built — the same recipe applies
to any other gene:

1. **Collect a reference CDS set.** Gather CDS sequences for your gene (and
   any close in-group paralogue you want to detect alongside it) across the
   clade you're mining — ideally one annotated ortholog per species from
   existing genome annotations.
2. **Filter out truncated annotations.** Pick a well-annotated reference
   species (e.g. chicken for birds) and keep only sequences that are at
   least ~80% of its CDS length, to drop obviously fragmentary gene models
   before they contaminate your miniprot query set.
3. **Translate** the filtered CDS to protein (e.g. EMBOSS `transeq`).
4. **Build the template alignment.** Align a representative subset with
   `muscle`/`mafft`, trim it (e.g. `trimal -gt 0.7`), then add the remaining
   sequences on top with `mafft --add --keeplength` so the trimmed column
   set is preserved.
5. **Add an outgroup.** `mafft --add` a handful of related but clearly
   distinct sequences (e.g. a different paralogue, or orthologs from a
   distantly related outgroup clade) purely to root the tree. These do
   *not* go into `reference_proteins` — only into the alignment.
6. **Build the tree:** `FastTree < alignment.aln > alignment.aln.treefile`.

You now have your `--reference_proteins` (pre-outgroup FASTA),
`--reference_alignment` and `--reference_tree`.

---

## 8. Known limitations

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
- No stop-codon / start-codon extension is attempted beyond what miniprot's
  own spliced alignment reports; a `functional_partial` call often just
  means the gene continues past the end of whichever reference protein
  produced the best-scoring alignment at that locus, not necessarily that
  the genome assembly itself is incomplete there.
