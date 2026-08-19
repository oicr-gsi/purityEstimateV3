# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-08-17
### Added
- `assess_primary_variants`, which does two things with the primary's somatic VCF.
  - **Triage.** Reports how many candidate sites survive each WISP filter in turn, so a
    primary too weak to support MRD is identifiable before a plasma run is committed.
    Provisioned as `primary_site_report`. `min_usable_sites` skips purity estimation below a
    threshold; the run still succeeds and provisions whatever was produced, so the report
    explaining the decision is always delivered. Defaults to 0, which reports without gating.
  - **Germline correction.** Removes normal-supported sites from the VCF that feeds
    SAGE_APPEND. WISP documents this filter but never applies it, because oncoanalyser does
    not pass the reference sample id through and WISP therefore has no genotype to test.
    The test is normal VAF above 1%, far below heterozygous frequency, so alongside germline
    variants it also catches low-level artefacts shared by tumour and normal. Either way the
    support does not come from the tumour, so counting them as tumour signal errs in the
    false-positive direction. Disable with `apply_germline_correction = false`. The quality
    field is read as `ARCBQ`, falling back to `ABQ` for VCFs from older versions; `RABQ` is
    rejected as the raw, pre-recalibration value. The uncorrected VCF is retained as
    `<sample>.purple.somatic.prefilter.vcf.gz`.
    The report gives both counts for this filter, with their denominators: what it removes
    from the whole VCF, and how much of that was still in the usable set. The two differ by
    a lot, because the triage table applies the filter last.

### Changed
- `pipeline_info` is now OPTIONAL. A run whose purity estimation is skipped launches no
  Nextflow stage at all, so there is no `pipeline_info` to provision. `primary_site_report`
  is the output that is always present.
- Two inputs renamed to drop site-specific names: `oicr_config` to `site_config`, and
  `qsub_wrapper` to `submit_wrapper`.
- Comments and `parameter_meta` no longer carry measured timings, past error strings or
  site-specific detail.

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
