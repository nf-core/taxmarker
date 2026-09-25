# AGENTS.md

This file provides guidance to coding agents (e.g. Claude Code) when working with code in this repository.

## Project Overview

**nf-core/taxmarker** is a Nextflow bioinformatics pipeline, a re-implementation of the Sativa algorithm (Kozlov et al. 2016), that identifies taxonomically mislabelled sequences by evolutionary placement: it builds a phylogeny from a declared taxonomy, places each sequence back into it after removing it (leave-one-out), and flags sequences whose phylogenetic signal doesn't agree with their declared taxonomy. It is built from the nf-core template (currently synced to v4.1.0) and uses Nextflow DSL2. The pipeline is currently in early development (v1.0.0dev, not yet released or transferred to the nf-core org).

Requires Nextflow ≥ 25.10.4.

## Commands

**Run the pipeline (with Docker):**

```bash
nextflow run main.nf -profile docker --sequences sequences.fasta --taxonomy taxonomy.tsv --outdir results
```

**Run minimal test suite:**

```bash
nextflow run main.nf -profile test,docker --outdir results
```

**Run nf-test (unit/integration tests):**

```bash
nf-test test tests/default.nf.test
# Run all tests (respects nf-test.config ignore rules):
nf-test test
```

**Lint with nf-core tools:**

```bash
nf-core pipelines lint
```

**Format code (Prettier + Nextflow lint via prek/pre-commit):**

```bash
prek run -a
```

**Update nf-core modules:**

```bash
nf-core modules update <module-name>
```

## Architecture

### Execution flow

```
main.nf
  └── PIPELINE_INITIALISATION   (subworkflows/local/utils_nfcore_taxmarker_pipeline/main.nf)
        validates params, resolves --taxonomy/--sequences into channels
  └── NFCORE_TAXMARKER
        └── TAXMARKER            (workflows/taxmarker.nf)  ← main logic lives here
              ├── RESOLVETAXONOMY       (modules/local/resolvetaxonomy/) -- from --taxonomy,
              │     or derived from --sequences record headers if omitted (GTDB-style)
              ├── CHECKNAMECONSISTENCY  (modules/local/checknameconsistency/) -- validates
              │     taxonomy/sequences names match, rewrites problematic characters
              ├── EMBOSS_SEQRET         (modules/nf-core/emboss/seqret/) -- normalises to FASTA
              ├── GTDBWEIGHTTRANSLATE   (modules/local/gtdbweighttranslate/) -- optional,
              │     --upstream 'GTDB' to enable; computes --sequence_weights automatically
              │     from GTDB genome metadata, runs on the pristine pre-RESOLVETAXONOMY
              │     sequences (needs the raw accession/header, not yet sanitised/stripped)
              ├── WEIGHTED_CLUSTERING   (subworkflows/local/weighted_clustering/) -- reduce
              │     the input to one representative per (cluster, taxon) pair before
              │     raxtax/alignment/placement all see it (see nf-core/taxmarker#15)
              ├── RAXTAX_PREFILTER      (subworkflows/local/raxtax_prefilter/) -- optional,
              │     --skip_raxtax to disable; fast raxtax self-classification triage on
              │     unaligned sequences that flags severely mislabeled sequences before
              │     alignment and the expensive placement step
              ├── ENSURE_ALIGNED        (subworkflows/local/ensure_aligned/) -- transparently
              │     aligns unaligned input via hmmalign (--hmm/--hmm_name); already-aligned
              │     input passes through unchanged
              ├── GAPFILTER / PROFILECOVER (modules/local/{gapfilter,profilecover}/) --
              │     drop sequences too short/incomplete to place reliably (whichever of
              │     ENSURE_ALIGNED's two branches ran); each has its own skip flag
              ├── SATIVA (subworkflows/local/sativa/) -- optional, --skip_sativa to disable;
              │     builds the reference tree with RAxML-NG, delegates leave-one-out
              │     placement/scoring to Auguste Gardette's sativa-epang fork (EPA-ng-based)
              └── MULTIQC               (modules/nf-core/multiqc/)
  └── PIPELINE_COMPLETION        (subworkflows/local/utils_nfcore_taxmarker_pipeline/main.nf)
        sends email / completion summary
```

Skipping `--skip_sativa` turns the rest of the pipeline into a general-purpose
taxonomy-resolution/alignment/prefilter QC tool; the raxtax prefilter and gap/profile-cover
filters still run as configured. See `README.md` for the full parameter list.

### Key conventions

- **Module arguments**: Pass extra CLI flags to tools via `ext.args` in `conf/modules.config`, not in the module itself.
- **Output paths**: Default publish rule in `conf/modules.config` derives directory from the process name (e.g., `RESOLVETAXONOMY` → `outdir/resolvetaxonomy/`). Override per-process with a `publishDir` block.
- **Input**: No samplesheet -- `--sequences` (phylip/clustal/fasta, aligned or not) and an optional `--taxonomy` (TSV; derived from `--sequences` headers if omitted, GTDB-style). Validated implicitly by `RESOLVETAXONOMY`/`CHECKNAMECONSISTENCY`, not a JSON schema file.
- **Parameter schema**: `nextflow_schema.json` defines all pipeline parameters and is used for CLI validation (via nf-schema plugin) and help text generation.
- **Software versions**: Collected via a `channel.topic("versions")` stream and written to `pipeline_info/nf_core_taxmarker_software_mqc_versions.yml` for MultiQC.
- **nf-core modules**: Modules under `modules/nf-core/` and subworkflows under `subworkflows/nf-core/` are managed by nf-core tools — do not edit them directly. Custom/local code goes in `subworkflows/local/`.
- **Tool identity vs. pipeline identity**: `subworkflows/local/sativa/`, `modules/local/sativascore`, `modules/local/sativaloosplit`, and `--skip_sativa` all name the wrapped SATIVA placement algorithm specifically (planned to be proposed upstream as its own nf-core subworkflow) — they are not renamed when the pipeline itself was renamed from `sativa` to `taxmarker`.

### Container registries

All container profiles (`docker`, `singularity`, `apptainer`, etc.) default to `quay.io` as the registry. The `wave` profile enables on-demand container building via Seqera Wave, required for ARM64.

### Test infrastructure

- `nf-test.config` defines test directories and triggers (files that force a full test run when changed).
- Tests run with `-profile test` by default; other pipeline-level profiles exercise specific input shapes: `test_fasta`, `test_clustal`, `test_gtdb`, `test_gtdb_unaligned`, `test_gtdb_embedded`, `test_full`.
- Test fixtures are fetched remotely from the `taxmarker` branch of `erikrikarddaniel/test-datasets` (`params.pipelines_testdata_base_path`) — no test data is committed to this repo.
- Snapshot files (`*.snap`) track expected outputs — update them with `nf-test test --update-snapshot` after intentional output changes.
