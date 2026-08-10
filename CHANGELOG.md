# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-08-10
### Added
- First version. Wraps nf-core/oncoanalyser 3.0.0-rc.3 for tumour-informed MRD
  detection, replacing the 2.3.0-based `purityEstimate` workflow. A single WDL covers
  three modes selected by the `mode` input: `WG` (primary tumour/normal WGTS calling),
  `PE` (WISP purity estimation of a longitudinal sample against a pre-existing WG
  tarball) and `WG_PE` (both, in sequence).
- `normal_alignments` is MANDATORY for WG and WG_PE, with no override. A tumour-only
  primary does not merely degrade the MRD result, it inverts it: with no matched normal
  SAGE cannot subtract germline variants, the primary somatic call set fills with germline
  sites, and WISP measures those in the patient's own cfDNA at heterozygous frequency and
  reports a large spurious tumour fraction. Measured on subject OCT_011303, one plasma
  against three primaries: NEGATIVE against the tumour/normal primary, 0.2255 and 0.2451
  against two tumour-only primaries.
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
