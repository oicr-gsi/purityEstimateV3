# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-08-19
### Added
- `pre_filtering`, which filters the primary's somatic VCF in WG and WG_PE mode and reports
  what each filter cost. Two filter sets, both primary-side: the **primary filters** WISP
  uses to decide which sites can carry MRD signal (mappability, repeat count, SNV only,
  tier, nearby indel, subclonal), and the **germline filter** WISP
  documents but never applies, because oncoanalyser does not pass the reference sample id
  through and WISP therefore has no genotype to test.
  The result is written as `<sample>.purple.somatic.prefiltered.vcf.gz` **beside the
  untouched original**, and both go into the WG archive. A reader of that archive therefore
  sees the full call set and the MRD-usable subset, which is what a clinical report needs in
  order to say whether a plasma sample is worth taking.
  Note that WISP applies the primary filters itself, on whatever sites it is given, and
  records the reason per site. The prefiltered VCF is a deliverable, not a correction; only
  the germline filter changes what WISP can see.
  Two of WISP's conditions are measurements on the cfDNA sample and so cannot be applied to
  the primary at all: average edge distance (`AED[1] >= 0.06`) and the quality ratio
  `(RC_QUAL[0]+[1]+[3])/(RC_CNT[0]+[1]+[3]) >= 18`. WISP applies both itself, so the reported
  count is an upper bound on the sites it will actually use, and the report says so.
- `primary_site_report`, the per-filter table: how many candidate sites survive each filter
  in turn, and the final count usable for MRD.
- `use_primary_filters`. In PE and WG_PE mode this chooses whether the plasma stage works
  from the prefiltered VCF or the full call set; no filtering happens in the PE stage.
  oncoanalyser resolves the somatic VCF by exact filename, so using the prefiltered one means
  staging a copy of the PURPLE directory in which that name holds the prefiltered content.
  That staging is done in `run_purity_estimate`, beside the samplesheet it feeds. A WG tarball
  produced before pre-filtering existed contains no prefiltered VCF; the run falls back to
  the full call set with a warning rather than failing.
- `pack_wgts`, which archives the WG results. Split out of `run_wgts` so that it runs after
  `pre_filtering` and can therefore carry both somatic VCFs. It stages symlinks and
  dereferences them into the archive, so the WG output is not copied a second time.

### Changed
- `run_purity_estimate` takes `primary_purple_dir` in place of `wgts_outdir`: the PURPLE
  directory itself rather than its parent, so a staged copy can be substituted.
- Comments and `parameter_meta` no longer carry measured numbers or past error strings.

### Removed
- `min_usable_sites`. The workflow does not decide whether purity estimation is worth
  running; that is a clinical decision taken from the WG deliverable. The site counts are
  still reported, they just do not gate anything.
- `apply_germline_correction`. Both filter sets always run in WG mode. What used to be
  optional is now the `use_primary_filters` choice on the PE side.

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
