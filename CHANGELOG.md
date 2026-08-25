# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-08-19
### Added
- samtools tasks cap their thread count at the CPUs actually allocated, and log a note when
  the two differ. A backend that does not map the `cpu` runtime attribute otherwise leaves
  samtools oversubscribed on one core, which is slower than running single-threaded and gives
  no indication why.
- The `@PG` header is now checked rather than assumed. `merge_bams` strips it only when it
  is actually malformed -- an `@PG` line carrying more than one `ID:` field, which is what an
  aligner embedding tabs in `CL:` produces and what htsjdk rejects -- so a well-formed header
  keeps its provenance. `sanitize_header` is replaced by `repair_flags` (the unpaired-flag
  repair that only the fixmate path needs) and `strip_bad_pg` (run the check at all), which
  were previously conflated.
- `probe_mc_tags`, which decides whether fixmate is needed from the alignments rather than
  from the input alone. `doFixmate = false` is honoured when the tags are really present; when
  they are not, fixmate runs anyway with a warning, rather than letting REDUX mark duplicates
  silently wrong. `doFixmate = true` always fixmates, so the probe only runs when its answer
  can change the decision.
- `doFixmate`. Consulted only when `run_redux` is true and the sample is Illumina: leave it
  true for alignments that lack mate CIGAR tags, set it false when the aligner already wrote
  them, as bwa-mem2 does, to skip the per-chromosome fixmate scatter. The inputs are still
  merged, since REDUX takes one alignment file per sample.
- `scheduler` (`sge` or `slurm`, validated) and `nextflow_config`. The submit wrapper is
  installed only for `sge`, where the executor embeds `h_rss`/`mem_free` directly in
  `.command.run`; the Slurm executor emits `--mem` and `--cpus-per-task` from the memory and
  cpus directives. `nextflow_config` is a list, each entry passed with its own `-c` in order,
  so a shared overlay and a per-site one compose. A `slurm` run must supply one, because the
  config the module ships selects the SGE executor.
- `pre_filtering`, which filters the primary's somatic VCF in WG and WG_PE mode and reports
  what each filter cost. Two filter sets, both primary-side: the **primary filters** WISP
  uses to decide which sites can carry MRD signal (mappability, repeat count, SNV only,
  tier, nearby indel, subclonal), and the **germline filter** which WISP
  documents but never applies.
  The result is written as `<sample>.purple.somatic.prefiltered.vcf.gz` **beside the
  untouched original**, and both go into the WG archive. Only the PASS step drops records;
  the rest add their name to the VCF FILTER column, so `FILTER="PASS"` selects the usable
  sites.
- `use_primary_filters`. In PE and WG_PE mode this chooses whether the plasma stage works
  from the prefiltered VCF or the full call set; no filtering happens in the PE stage.
- `pack_wgts`, which archives the WG results. Split out of `run_wgts` so that it runs after
  `pre_filtering` and can therefore carry both somatic VCFs. It stages symlinks and
  dereferences them into the archive, so the WG output is not copied a second time.

## [1.0.0] - 2026-08-10
### Added
- First version. Wraps nf-core/oncoanalyser 3.0.0-rc.3 for tumour-informed MRD
  detection, replacing the 2.3.0-based `purityEstimate` workflow. A single WDL covers
  three modes selected by the `mode` input: `WG` (primary tumour/normal WGTS calling),
  `PE` (WISP purity estimation of a longitudinal sample against a pre-existing WG
  tarball) and `WG_PE` (both, in sequence).
- `normal_alignments` is MANDATORY for WG and WG_PE, with no override. A tumour-only
  primary does not merely degrade the MRD result, it inverts it: with no matched normal
  SAGE cannot subtract germline variants, the somatic call set fills with germline sites,
  and WISP measures those in the patient's own cfDNA and reports a large spurious tumour
  fraction. Measured, and confirmed by re-running the same primaries with their normals.
- Resource paths come from the `oncoanalyser/3.0.0-rc.3` environment module rather than a
  hardcoded root: `$IMAGES_DIR`, `$ONCOANALYSER_FOLDER`, `$NEXTFLOW_HOME`, `$QSUB_WRAPPER`,
  `$ONCOANALYSER_OICR_CONFIG`, and `nextflow` from PATH (the module ships 25.10.4). Reference
  data comes from the separate `oncoanalyser-data/3.0.0` module via `$ONCOANALYSER_DATA_ROOT`.
  The `resources_root` input is gone and no absolute paths remain in the workflow.
- Accepts BAM or CRAM for every sample. CRAMs are decoded against hs38DH, the reference
  the Ultima vendor encodes with; multiple alignments per sample are merged.
- Per-sample sequencing platform detection from the `@RG PL` tag, overridable with
  `sequencing_platform`. Because oncoanalyser applies one `--sequencing_platform` per
  pipeline run, a WG run whose tumour and normal disagree fails before Nextflow starts
  unless `allow_mixed_platforms` is set.
- Illumina-only per-chromosome `fixmate` before REDUX, skipped for single-end platforms.
- `nextflow_stub` runs oncoanalyser with `-stub --create_stub_placeholders`, exercising
  the samplesheet, output layout and tarring in minutes.
- Germline variant calls are excluded from the WG tarball by default
  (`include_germline_outputs`). They play no part in MRD: WISP reads the PURPLE somatic
  VCF, and oncoanalyser's own `purity_estimate.nf` disables germline calling. They are
  still generated, because oncoanalyser hardcodes germline calling on.
- The PE stage does not take the primary normal. It reaches WISP's LOH path only as a
  `redux_dir`, never as a BAM, so supplying it cost a second REDUX plus AMBER (about
  5.3 hours) and produced no LOH output. As a side effect COBALT runs tumour-only
  against the diploid-regions BED in the PE stage; this moves the CNV number by about
  2% and does not touch the SNV path, which is what the assay reports.
