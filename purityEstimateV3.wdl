version 1.0

# Alignment file plus its index. Accepts BAM or CRAM: the extension is load-bearing
# downstream (see cram_to_bam and the samplesheet staging in the head-job tasks).
struct Alignment {
    File aln
    File idx
}

# purityEstimateV3 - oncoanalyser 3.0.0-rc.3 wrapper, Illumina / SBX / Ultima.
#
# User-facing docs (modes, valid input combinations, the one-platform-per-run rule,
# recommended production shapes) live in meta.description at the bottom of this file, which is
# what generates README.md. Keep them there, not here. Notes below are for maintainers.
#
# Forced by oncoanalyser 3.x, and the reason this is not a copy of 2.x purityEstimate.wdl
# (devlog.txt sections 4, 9 and 14 have the detail):
#   * --skip_msi_jitter is gone. run_redux=false now means filetype bam_redux plus info
#     generate_redux_tsvs_only, with redux STILL in --processes_manual, so REDUX regenerates
#     the TSVs instead of us faking them with touch.
#   * --sequencing_platform must be passed for SBX / Ultima.
#   * CRAM is converted to BAM against the reference it was ENCODED with. For Ultima that is
#     hs38DH, which differs from oncoanalyser's masked GRCh38 on chr9, chr13, chr16, chr21 and
#     chrX, so decoding against the wrong one silently corrupts bases.
#   * SAGE_APPEND is directly selectable, so 'orange' is no longer needed to trigger it.
#
# Why platform is resolved per sample rather than per run: fixmate_chr keeps only records with
# $7 == "=" and fixmate_discordant requires the 0x1 paired flag, so running it on single-end
# Ultima reads drops EVERY record and leaves a header-only BAM. COBALT then dies inside
# ReadDepthStatisticsNormaliser with a bare IllegalStateException and no useful message. Each
# head job is therefore passed the platform of the sample it processes: the primary for WG, the
# longitudinal for PE, the control BAM for a control run.
#
# Upstream limitation behind allow_mixed_platforms: oncoanalyser has one --sequencing_platform
# per run, so a PE run that also processes a normal of a different platform cannot be expressed
# correctly. See devlog Q6.
#
# Note on run_redux and controls: control BAMs are legacy GATK-processed (deduped + BQSR) and
# cannot be redux-processed. Controls always run with run_redux=false regardless of the
# workflow-level run_redux setting. Only subject alignments are affected by run_redux.

workflow purityEstimateV3 {
    input {
        String        mode                    # "WG", "PE", or "WG_PE"
        String        group_id
        String        subject_id
        Array[Alignment]? tumor_alignments        # required for WG mode
        Array[Alignment]? normal_alignments       # optional throughout
        Array[Alignment]? longitudinal_alignments # required for PE mode
        File?         wgts_tarball            # required for PE mode; produced by a prior WG run
        String?       tumor_sample_id         # required for PE mode; overrides @RG SM in WG modes
        String?       sequencing_platform     # "illumina" | "sbx" | "ultima"; null = detect from @RG PL
        Boolean       run_redux = false       # see note above; semantics changed in 3.x
        Boolean       run_control = false     # also estimate purity for each control (PE/WG_PE only)
        Array[Pair[String, String]]? controls # PE mode: (control_id, bam_path) pairs; BAI at bam_path+".bai"
        String        outputFileNamePrefix
        Boolean       allow_mixed_platforms = false  # production guardrail; see parameter_meta
        Boolean       nextflow_stub = false   # run oncoanalyser with -stub; see parameter_meta
        String        resources_root = "/.mounts/labs/gsi/src/hmftools"
        String        modules = "java/17 singularity/3.9.4 samtools/1.16.1"
    }

    parameter_meta {
        mode:                   "Run mode: WG=WGTS only; PE=purity estimation against a pre-existing WG tarball; WG_PE=run both sequentially"
        group_id:               "Sample group identifier used as the output subdirectory name and samplesheet group_id"
        subject_id:             "Subject/patient identifier"
        tumor_alignments:       "Primary tumour BAM/CRAM(s) with indices; required for WG and WG_PE mode. Multiple entries (e.g. per flowcell) are merged"
        normal_alignments:      "Matched normal BAM/CRAM(s) with indices; optional. Omitting it gives tumour-only WGTS, and LOH-free purity estimation"
        longitudinal_alignments: "ctDNA BAM/CRAM(s) with indices; required for PE and WG_PE mode. Multiple entries are merged"
        wgts_tarball:           "Tarball produced by a prior WG run (amber/, cobalt/, purple/, pave/, sage/); required for PE mode"
        tumor_sample_id:        "Primary tumour sample ID. Normally leave unset: in WG modes it is read from the tumour @RG SM tag, and in PE mode it is derived from the WG tarball's purple/ filenames. If given, it overrides the @RG SM tag in WG modes and is cross-checked against the tarball in PE mode"
        sequencing_platform:    "Sequencing platform passed to --sequencing_platform: illumina, sbx or ultima. Leave unset to detect it from the @RG PL tag"
        run_redux:              "When true, REDUX processes the alignments normally. When false, inputs are treated as already REDUX-processed and REDUX only regenerates its TSVs (-bqr_jitter_msi_only). Does not affect control BAMs"
        run_control:            "When true, also run purity estimation for each control BAM in the controls array (PE/WG_PE mode only)"
        controls:               "PE mode: array of (control_id, bam_path) pairs; BAI assumed at bam_path+'.bai'. Controls are BAM only and always run with run_redux=false; used only when run_control=true"
        outputFileNamePrefix:   "Output directory prefix; the pipeline writes to outputFileNamePrefix/group_id/"
        allow_mixed_platforms:  "Guardrail. oncoanalyser applies ONE --sequencing_platform per pipeline run, so two samples of different platforms in the same run means one of them gets the wrong error model. Left false (the production default) such a combination fails before Nextflow starts. Set true only for deliberate experiments: the run proceeds with a loud warning"
        nextflow_stub:          "When true, oncoanalyser runs with -stub --create_stub_placeholders: every process writes placeholder outputs instead of doing real work. This exercises the whole wrapper (samplesheet, samplesheet validation, output layout, tarring, Vidarr outputs) in minutes. Real input alignments are still required, but they can be tiny, because the pipeline never reads them. Not a Cromwell dry run"
        resources_root:         "Root of the staged oncoanalyser 3.x resources. Every resource path below is derived from it, so a future oncoanalyser/3.x module needs only these defaults changed"
        modules:                "Environment modules to load. Deliberately excludes the oncoanalyser 2.x module: its lib/ shadows the system libcurl and its bundled Nextflow (25.04.3) is too old for this pipeline"
    }

    # Staged resources, all under resources_root. String (not File) so Cromwell reads them
    # from the shared filesystem instead of localizing gigabytes per task.
    String ref_data_dir   = resources_root + "/ref_data_3.0.0"
    String images_dir     = resources_root + "/oncoanalyser3_images"
    String pipeline_dir   = resources_root + "/oncoanalyser_3.0.0-rc.3"
    String nextflow_bin   = resources_root + "/nextflow_25.10.4/nextflow"
    String nextflow_home  = resources_root + "/nextflow_home_3.0.0"
    String oicr_config    = resources_root + "/oncoanalyser_oicr.config"
    String qsub_wrapper   = resources_root + "/qsub_wrapper.sh"
    String cram_reference = resources_root + "/hs38DH/GRCh38_full_analysis_set_plus_decoy_hla.fa"

    # GRCh38 chromosomes processed individually to keep fixmate memory small.
    Array[String] chromosomes = [
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
        "chr21","chr22","chrX","chrY","chrM"
    ]

    # Fail fast and legibly on a mode/input mismatch. Without this the first symptom is
    # "select_first was called with 1 empty values" from a downstream call, which says
    # nothing about which input was missing or why.
    call validate_inputs {
        input:
            mode                = mode,
            has_tumor           = defined(tumor_alignments),
            has_normal          = defined(normal_alignments),
            has_longitudinal    = defined(longitudinal_alignments),
            has_wgts_tarball    = defined(wgts_tarball),
            has_tumor_sample_id = defined(tumor_sample_id),
            run_control         = run_control,
            has_controls        = defined(controls)
    }

    # Read sample ID and platform from one header per sample. Cheap: the header is read
    # straight off shared storage, nothing is localized.
    if (defined(tumor_alignments)) {
        call get_alignment_info as tumor_info {
            input: aln_path = select_first([tumor_alignments])[0].aln
        }
    }
    if (defined(normal_alignments)) {
        call get_alignment_info as normal_info {
            input: aln_path = select_first([normal_alignments])[0].aln
        }
    }
    if (defined(longitudinal_alignments)) {
        call get_alignment_info as longitudinal_info {
            input: aln_path = select_first([longitudinal_alignments])[0].aln
        }
    }

    # PLATFORM IS PER SAMPLE, NOT PER RUN. A run can legitimately mix technologies -- an
    # Illumina primary tumour with an Ultima ctDNA sample is the whole point of this project --
    # and fixmate must be decided from each sample's OWN @RG PL.
    #
    # Getting this wrong is silent and destructive: fixmate_chr keeps only records with
    # $7 == "=" and fixmate_discordant requires the 0x1 paired flag, but Ultima reads are
    # single-end with RNEXT "*" and no 0x1, so EVERY record is dropped and the merged BAM ends
    # up header-only. COBALT then dies inside ReadDepthStatisticsNormaliser with a bare
    # IllegalStateException and no useful message. An explicit sequencing_platform input still
    # overrides everything, for data with a missing or wrong PL tag.
    if (defined(tumor_alignments)) {
        String tumor_platform  = select_first([sequencing_platform, tumor_info.platform])
        Boolean tumor_fixmate  = run_redux && tumor_platform == "illumina"
    }
    if (defined(normal_alignments)) {
        String normal_platform = select_first([sequencing_platform, normal_info.platform])
        Boolean normal_fixmate = run_redux && normal_platform == "illumina"
    }
    if (defined(longitudinal_alignments)) {
        String longitudinal_platform = select_first([sequencing_platform, longitudinal_info.platform])
        Boolean longitudinal_fixmate = run_redux && longitudinal_platform == "illumina"
    }

    # The tumour sample ID is resolved per mode, NOT here. A workflow-level declaration is
    # evaluated eagerly by Cromwell even when only used inside a conditional, so
    # select_first([tumor_sample_id, tumor_info.sample_id]) at this scope blows up in PE mode
    # where both are legitimately absent. WG modes derive it in the WG block below (@RG SM or
    # an explicit override); PE derives it from the tarball in extract_wgts.

    # ---------------------------------------------------------------------------------
    # Stage each sample: CRAM -> BAM, then optional fixmate, then merge when needed.
    # ---------------------------------------------------------------------------------

    if (defined(tumor_alignments)) {
        scatter (a in select_first([tumor_alignments])) {
            if (basename(a.aln, ".cram") != basename(a.aln)) {
                call cram_to_bam as cram_tumor {
                    input: aln = a.aln, idx = a.idx, cram_reference = cram_reference
                }
            }
            File tumor_staged_bam = select_first([cram_tumor.bam, a.aln])
            File tumor_staged_bai = select_first([cram_tumor.bai, a.idx])
        }

        if (select_first([tumor_fixmate])) {
            scatter (p in zip(tumor_staged_bam, tumor_staged_bai)) {
                scatter (chr in chromosomes) {
                    call fixmate_chr as fixmate_tumor_chr {
                        input: bam = p.left, bai = p.right, chr = chr
                    }
                }
                call fixmate_discordant as fixmate_tumor_disc {
                    input: bam = p.left, bai = p.right
                }
            }
            call merge_bams as merge_tumor_fixmate {
                input:
                    bams = flatten([flatten(fixmate_tumor_chr.fixed_bam), fixmate_tumor_disc.fixed_bam]),
                    bais = flatten([flatten(fixmate_tumor_chr.fixed_bai), fixmate_tumor_disc.fixed_bai]),
                    sanitize_header = true
            }
        }

        # Nothing to merge for a single already-fixmate-free alignment; use it as is.
        if (!select_first([tumor_fixmate]) && length(tumor_staged_bam) > 1) {
            call merge_bams as merge_tumor_plain {
                input:
                    bams = tumor_staged_bam,
                    bais = tumor_staged_bai,
                    sanitize_header = false
            }
        }

        File tumor_bam = select_first([merge_tumor_fixmate.bam, merge_tumor_plain.bam, tumor_staged_bam[0]])
        File tumor_bai = select_first([merge_tumor_fixmate.bai, merge_tumor_plain.bai, tumor_staged_bai[0]])
    }

    if (defined(normal_alignments)) {
        scatter (a in select_first([normal_alignments])) {
            if (basename(a.aln, ".cram") != basename(a.aln)) {
                call cram_to_bam as cram_normal {
                    input: aln = a.aln, idx = a.idx, cram_reference = cram_reference
                }
            }
            File normal_staged_bam = select_first([cram_normal.bam, a.aln])
            File normal_staged_bai = select_first([cram_normal.bai, a.idx])
        }

        if (select_first([normal_fixmate])) {
            scatter (p in zip(normal_staged_bam, normal_staged_bai)) {
                scatter (chr in chromosomes) {
                    call fixmate_chr as fixmate_normal_chr {
                        input: bam = p.left, bai = p.right, chr = chr
                    }
                }
                call fixmate_discordant as fixmate_normal_disc {
                    input: bam = p.left, bai = p.right
                }
            }
            call merge_bams as merge_normal_fixmate {
                input:
                    bams = flatten([flatten(fixmate_normal_chr.fixed_bam), fixmate_normal_disc.fixed_bam]),
                    bais = flatten([flatten(fixmate_normal_chr.fixed_bai), fixmate_normal_disc.fixed_bai]),
                    sanitize_header = true
            }
        }

        if (!select_first([normal_fixmate]) && length(normal_staged_bam) > 1) {
            call merge_bams as merge_normal_plain {
                input:
                    bams = normal_staged_bam,
                    bais = normal_staged_bai,
                    sanitize_header = false
            }
        }

        File normal_bam = select_first([merge_normal_fixmate.bam, merge_normal_plain.bam, normal_staged_bam[0]])
        File normal_bai = select_first([merge_normal_fixmate.bai, merge_normal_plain.bai, normal_staged_bai[0]])
    }

    if (defined(longitudinal_alignments)) {
        scatter (a in select_first([longitudinal_alignments])) {
            if (basename(a.aln, ".cram") != basename(a.aln)) {
                call cram_to_bam as cram_longitudinal {
                    input: aln = a.aln, idx = a.idx, cram_reference = cram_reference
                }
            }
            File longitudinal_staged_bam = select_first([cram_longitudinal.bam, a.aln])
            File longitudinal_staged_bai = select_first([cram_longitudinal.bai, a.idx])
        }

        if (select_first([longitudinal_fixmate])) {
            scatter (p in zip(longitudinal_staged_bam, longitudinal_staged_bai)) {
                scatter (chr in chromosomes) {
                    call fixmate_chr as fixmate_longitudinal_chr {
                        input: bam = p.left, bai = p.right, chr = chr
                    }
                }
                call fixmate_discordant as fixmate_longitudinal_disc {
                    input: bam = p.left, bai = p.right
                }
            }
            call merge_bams as merge_longitudinal_fixmate {
                input:
                    bams = flatten([flatten(fixmate_longitudinal_chr.fixed_bam), fixmate_longitudinal_disc.fixed_bam]),
                    bais = flatten([flatten(fixmate_longitudinal_chr.fixed_bai), fixmate_longitudinal_disc.fixed_bai]),
                    sanitize_header = true
            }
        }

        if (!select_first([longitudinal_fixmate]) && length(longitudinal_staged_bam) > 1) {
            call merge_bams as merge_longitudinal_plain {
                input:
                    bams = longitudinal_staged_bam,
                    bais = longitudinal_staged_bai,
                    sanitize_header = false
            }
        }

        File longitudinal_bam = select_first([merge_longitudinal_fixmate.bam, merge_longitudinal_plain.bam, longitudinal_staged_bam[0]])
        File longitudinal_bai = select_first([merge_longitudinal_fixmate.bai, merge_longitudinal_plain.bai, longitudinal_staged_bai[0]])
    }

    # ---------------------------------------------------------------------------------
    # Head jobs
    # ---------------------------------------------------------------------------------

    if (mode == "WG" || mode == "WG_PE") {
        # Safe at this scope: WG and WG_PE both require tumor_alignments, so tumor_info ran.
        String wg_tumor_sample_id = select_first([tumor_sample_id, tumor_info.sample_id])

        call run_wgts {
            input:
                group_id            = group_id,
                subject_id          = subject_id,
                tumor_bam           = select_first([tumor_bam]),
                tumor_bai           = select_first([tumor_bai]),
                tumor_sample_id     = wg_tumor_sample_id,
                normal_bam          = normal_bam,
                normal_bai          = normal_bai,
                normal_sample_id    = normal_info.sample_id,
                normal_platform     = normal_platform,
                run_redux           = run_redux,
                nextflow_stub       = nextflow_stub,
                allow_mixed_platforms = allow_mixed_platforms,
                sequencing_platform = select_first([tumor_platform]),
                outdir              = outputFileNamePrefix,
                ref_data_dir        = ref_data_dir,
                images_dir          = images_dir,
                pipeline_dir        = pipeline_dir,
                nextflow_bin        = nextflow_bin,
                nextflow_home       = nextflow_home,
                oicr_config         = oicr_config,
                qsub_wrapper        = qsub_wrapper,
                modules             = modules
        }
    }

    # PE mode: extract the supplied tarball to get a local wgts directory.
    # WG_PE mode: the directory is already on the filesystem from run_wgts.
    if (mode == "PE") {
        call extract_wgts {
            input:
                tarball                  = select_first([wgts_tarball]),
                expected_tumor_sample_id = tumor_sample_id
        }
    }

    if (mode == "PE" || mode == "WG_PE") {
        String effective_wgts_dir = select_first([run_wgts.output_dir, extract_wgts.output_dir])

        # PE: authoritative value read off the tarball's purple/ filenames.
        # WG_PE: the WG step just wrote those files named for wg_tumor_sample_id.
        String pe_tumor_sample_id = select_first([extract_wgts.tumor_sample_id, wg_tumor_sample_id])

        call run_purity_estimate as subject_purity {
            input:
                group_id               = group_id,
                subject_id             = subject_id,
                tumor_sample_id        = pe_tumor_sample_id,
                normal_bam             = normal_bam,
                normal_bai             = normal_bai,
                normal_sample_id       = normal_info.sample_id,
                normal_platform        = normal_platform,
                longitudinal_bam       = select_first([longitudinal_bam]),
                longitudinal_bai       = select_first([longitudinal_bai]),
                longitudinal_sample_id = select_first([longitudinal_info.sample_id]),
                run_redux              = run_redux,
                nextflow_stub          = nextflow_stub,
                allow_mixed_platforms  = allow_mixed_platforms,
                sequencing_platform    = select_first([longitudinal_platform]),
                wgts_outdir            = effective_wgts_dir,
                outdir                 = outputFileNamePrefix,
                ref_data_dir           = ref_data_dir,
                images_dir             = images_dir,
                pipeline_dir           = pipeline_dir,
                nextflow_bin           = nextflow_bin,
                nextflow_home          = nextflow_home,
                oicr_config            = oicr_config,
                qsub_wrapper           = qsub_wrapper,
                modules                = modules
        }

        # Controls are legacy GATK-processed BAMs; always skip redux regardless of run_redux.
        if (run_control) {
            scatter (control in select_first([controls, []])) {
                call get_alignment_info as control_info {
                    input: aln_path = control.right
                }
                call run_purity_estimate as control_purity {
                    input:
                        group_id               = group_id,
                        subject_id             = control.left,
                        tumor_sample_id        = pe_tumor_sample_id,
                        normal_bam             = normal_bam,
                        normal_bai             = normal_bai,
                        normal_sample_id       = normal_info.sample_id,
                        normal_platform        = normal_platform,
                        longitudinal_bam       = control.right,
                        longitudinal_bai       = control.right + ".bai",
                        longitudinal_sample_id = control_info.sample_id,
                        run_redux              = false,
                        nextflow_stub          = nextflow_stub,
                        allow_mixed_platforms  = allow_mixed_platforms,
                        sequencing_platform    = select_first([sequencing_platform, control_info.platform]),
                        wgts_outdir            = effective_wgts_dir,
                        outdir                 = outputFileNamePrefix,
                        ref_data_dir           = ref_data_dir,
                        images_dir             = images_dir,
                        pipeline_dir           = pipeline_dir,
                        nextflow_bin           = nextflow_bin,
                        nextflow_home          = nextflow_home,
                        oicr_config            = oicr_config,
                        qsub_wrapper           = qsub_wrapper,
                        modules                = modules
                }
            }
        }

        call collect_results {
            input:
                outputFileNamePrefix = outputFileNamePrefix,
                subject_summary   = subject_purity.wisp_summary,
                subject_tarball   = subject_purity.wisp_tarball,
                control_summaries = select_first([control_purity.wisp_summary, []]),
                control_tarballs  = select_first([control_purity.wisp_tarball, []])
        }
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "Runs HMF oncoanalyser 3.0.0-rc.3 to estimate tumour purity in longitudinal ctDNA samples, for Illumina, Roche SBX or Ultima Genomics data. In WG mode it runs WGTS (REDUX, AMBER, COBALT, SAGE, PAVE, PURPLE) on a primary tumour with an optional matched normal and produces a tarball of the results. In PE mode it runs WISP against a pre-existing WG tarball to report the ctDNA fraction of a longitudinal sample. WG_PE does both in sequence. BAM and CRAM are both accepted.\n\n### One sequencing platform per pipeline run\n\noncoanalyser applies a single --sequencing_platform to a whole pipeline run and never checks it against the BAM headers, so two samples of different platforms in one run means one of them is analysed with the wrong error model, silently.\n\nNote that a run here means one oncoanalyser (Nextflow) invocation, not one WDL job. WG_PE launches two runs, so it can legitimately span platforms: the WG run uses the primary's platform and the PE run uses the longitudinal sample's.\n\nPlatform is read per sample from the @RG PL tag; `sequencing_platform` overrides it for data with a missing or wrong tag. Two consequences are enforced:\n\n* fixmate is applied only to Illumina samples. Ultima reads are single-end, and fixmate would drop every record and leave a header-only BAM.\n* if samples inside one pipeline run disagree on platform, the run is refused before Nextflow starts, unless `allow_mixed_platforms` is true.\n\n### Valid input combinations\n\n| mode | tumor_alignments | normal_alignments | longitudinal_alignments | wgts_tarball |\n|---|---|---|---|---|\n| WG | required | optional | - | - |\n| PE | - | optional | required | required |\n| WG_PE | required | optional | required | - |\n\nLeave `tumor_sample_id` unset. In WG modes it comes from the tumour @RG SM tag; in PE mode it is derived from the tarball's purple/ filenames, and a supplied value that disagrees is rejected. A wrong value otherwise fails deep inside SAGE_APPEND, after REDUX and COBALT have already run.\n\nOmitting `normal_alignments` costs LOH: purity then comes from COBALT and SNVs only, with no amber_dir input.\n\n### Recommended production shapes\n\n* Single platform with matched normal: WG_PE with tumour, normal and longitudinal samples, giving the full LOH-aware estimate.\n* Illumina primary with Ultima ctDNA: run WG on the Illumina pair, then PE against that tarball with the Ultima longitudinal sample and no normal_alignments, so each pipeline run stays single-platform. This is the validated configuration."
        dependencies: [
            {
                name: "oncoanalyser/3.0.0-rc.3",
                url: "https://github.com/nf-core/oncoanalyser"
            },
            {
                name: "nextflow/25.10.4",
                url: "https://www.nextflow.io"
            },
            {
                name: "samtools/1.16.1",
                url: "https://github.com/samtools/samtools"
            },
            {
                name: "hs38DH (GRCh38 full analysis set plus decoy plus HLA)",
                url: "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
            }
        ]
        output_meta: {
            wg_tarball: {
                description: "Tarball of oncoanalyser WGTS outputs (amber/, cobalt/, purple/, pave/, sage/) for the primary tumour sample; produced in WG and WG_PE mode. Used as the wgts_tarball input for a subsequent PE run. REDUX alignments are deliberately excluded to keep the archive small.",
                vidarr_label: "wgTarball"
            },
            wisp_tarballs: {
                description: "Combined tarball of WISP output directories for the subject longitudinal sample and all control samples; produced in PE and WG_PE mode.",
                vidarr_label: "wispTarballs"
            },
            wisp_summary: {
                description: "TSV file with one header row and one data row per sample (subject + controls) showing the WISP-estimated ctDNA purity fraction; produced in PE and WG_PE mode.",
                vidarr_label: "wispSummary"
            },
            pipeline_info: {
                description: "Tarball of the Nextflow pipeline_info/ directory (execution report, timeline, trace, DAG, params JSON, software versions) from the primary nextflow run; always produced. In WG_PE mode the PE run's pipeline_info is used.",
                vidarr_label: "pipelineInfo"
            }
        }
    }

    output {
        File? wg_tarball    = run_wgts.wgts_tarball    # produced in WG and WG_PE mode only
        File? wisp_tarballs = collect_results.wisp_tarballs
        File? wisp_summary  = collect_results.wisp_summary
        File  pipeline_info = select_first([subject_purity.pipeline_info_tarball, run_wgts.pipeline_info_tarball])
    }
}

# Check that the supplied inputs match the requested mode, and say precisely what is
# missing if they do not. Runs first and costs nothing.
task validate_inputs {
    input {
        String  mode
        Boolean has_tumor
        Boolean has_normal
        Boolean has_longitudinal
        Boolean has_wgts_tarball
        Boolean has_tumor_sample_id
        Boolean run_control
        Boolean has_controls
        Int memory  = 1
        Int timeout = 1
    }

    parameter_meta {
        mode:                "Requested run mode"
        has_tumor:           "Whether tumor_alignments was supplied"
        has_normal:          "Whether normal_alignments was supplied"
        has_longitudinal:    "Whether longitudinal_alignments was supplied"
        has_wgts_tarball:    "Whether wgts_tarball was supplied"
        has_tumor_sample_id: "Whether tumor_sample_id was supplied"
        run_control:         "Whether control purity estimation was requested"
        has_controls:        "Whether the controls array was supplied"
        memory:              "Memory in GB"
        timeout:             "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # Bind interpolated values to shell variables: WDL substitutes them at generation
      # time, so testing them directly produces a script full of constant comparisons.
      mode="~{mode}"
      errors=()

      case "${mode}" in
        WG)
          ~{if has_tumor then "true" else "false"} || errors+=("mode WG requires tumor_alignments")
          ;;
        WG_PE)
          ~{if has_tumor then "true" else "false"} || errors+=("mode WG_PE requires tumor_alignments")
          ~{if has_longitudinal then "true" else "false"} || errors+=("mode WG_PE requires longitudinal_alignments")
          ;;
        PE)
          ~{if has_longitudinal then "true" else "false"} || errors+=("mode PE requires longitudinal_alignments")
          ~{if has_wgts_tarball then "true" else "false"} || errors+=("mode PE requires wgts_tarball from a prior WG run")
          ;;
        *)
          errors+=("mode must be WG, PE or WG_PE, got '${mode}'")
          ;;
      esac

      if ~{run_control}; then
        ~{if has_controls then "true" else "false"} || errors+=("run_control=true requires the controls array")
        case "${mode}" in
          PE|WG_PE) ;;
          *) errors+=("run_control is only meaningful in PE or WG_PE mode") ;;
        esac
      fi

      if [ "${#errors[@]}" -gt 0 ]; then
        echo "ERROR: inputs do not satisfy mode ${mode}:" >&2
        for e in "${errors[@]}"; do
          echo "  - ${e}" >&2
        done
        exit 1
      fi

      # Advisory only: these combinations run, but not the way people usually expect.
      if [ "${mode}" != "PE" ] && ~{if has_normal then "false" else "true"}; then
        echo "note: no normal_alignments, so WGTS runs tumour-only and purity estimation has no LOH input" >&2
      fi

      echo "inputs OK for mode ${mode}"
    >>>

    output {
        String checked = read_string(stdout())
    }

    runtime {
        cpu:     1
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
    }
}

# Read the sample ID and sequencing platform out of an alignment header.
#
# aln_path is a String, not a File, so Cromwell does not localize a multi-hundred-GB
# BAM/CRAM just to read a few hundred header lines. samtools reads the header of a CRAM
# without needing the reference.
task get_alignment_info {
    input {
        String aln_path
        String modules = "samtools/1.16.1"
        Int memory  = 4
        Int timeout = 1
    }

    parameter_meta {
        aln_path: "Path to a BAM or CRAM on shared storage. Passed as a String so the file is read in place rather than localized"
        modules:  "Environment modules to load (samtools required)"
        memory:   "Memory in GB"
        timeout:  "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      samtools view -H ~{aln_path} > header.txt

      # Sample ID from the first @RG SM tag.
      awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/){sub(/^SM:/,"",$i); print $i; exit}}' \
        header.txt > sample_id.txt
      [ -s sample_id.txt ] || { echo "ERROR: no @RG SM tag in ~{aln_path}" >&2; exit 1; }

      # Refuse a file carrying several different samples: taking the first SM would be wrong.
      distinct=$(awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^SM:/) print $i}' header.txt | sort -u | wc -l)
      if [ "${distinct}" -ne 1 ]; then
        echo "ERROR: ${distinct} distinct @RG SM values in ~{aln_path}; expected exactly 1" >&2
        exit 1
      fi

      # Platform from @RG PL, lowercased to match --sequencing_platform (ULTIMA -> ultima).
      # Absent PL is treated as Illumina, which is also oncoanalyser's own default.
      awk '/^@RG/{for(i=1;i<=NF;i++) if($i~/^PL:/){sub(/^PL:/,"",$i); print tolower($i); exit}}' \
        header.txt > platform.txt
      [ -s platform.txt ] || echo illumina > platform.txt
    >>>

    output {
        String sample_id = read_string("sample_id.txt")
        String platform  = read_string("platform.txt")
    }

    runtime {
        cpu:     1
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}

# Convert a CRAM to BAM using the reference the CRAM was ENCODED against.
#
# This is mandatory, not defensive. CRAM stores read bases as differences from the
# reference, so decoding with a different reference silently yields wrong sequences with no
# error raised. Ultima vendor CRAMs are encoded against hs38DH, which differs from
# oncoanalyser's GRCh38_masked_exclusions_alts_hlas.fasta on chr9, chr13, chr16, chr21 and
# chrX (20 of 25 primary chromosomes are identical; those five are not). Converting to BAM
# makes the sequences explicit, after which the reference difference no longer matters.
task cram_to_bam {
    input {
        File   aln
        File   idx
        String cram_reference
        String modules = "samtools/1.16.1"
        Int threads = 8
        Int memory  = 16
        Int timeout = 24
    }

    parameter_meta {
        aln:            "Input CRAM file"
        idx:            "CRAM index (.crai) for aln"
        cram_reference: "Path to the FASTA the CRAM was encoded against, e.g. hs38DH for Ultima vendor CRAMs. A String so it is read from shared storage without localization"
        modules:        "Environment modules to load (samtools required)"
        threads:        "Number of samtools threads"
        memory:         "Memory in GB"
        timeout:        "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # Cromwell localizes the CRAM and its index into different directories, so symlink
      # both here with matching basenames for samtools to find the index.
      ln -s ~{aln} ./input.cram
      ln -s ~{idx} ./input.cram.crai

      [ -s "~{cram_reference}" ] || {
        echo "ERROR: CRAM reference not readable: ~{cram_reference}" >&2
        exit 1
      }

      samtools view -@ ~{threads} -T ~{cram_reference} -b -o converted.bam ./input.cram
      samtools index -@ ~{threads} converted.bam
    >>>

    output {
        File bam = "converted.bam"
        File bai = "converted.bam.bai"
    }

    runtime {
        cpu:     threads
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}

# Add mate CIGAR (MC) tags to one chromosome slice of one lane BAM.
# Handles concordant pairs only (RNEXT == "="); discordant pairs are handled
# by fixmate_discordant, which runs in parallel on the same lane BAM.
task fixmate_chr {
    input {
        File        bam
        File        bai
        String      chr
        String      modules = "samtools/1.16.1"
        Int threads = 4
        Int memory  = 16
        Int timeout = 4
    }

    parameter_meta {
        bam:       "BAM file for one lane; the original filename is preserved via symlink so samtools can locate the index automatically"
        bai:       "Index for bam"
        chr:       "Chromosome name (e.g. chr1) to extract and fixmate; one task is run per chromosome per lane"
        modules:   "Environment modules to load (samtools required)"
        threads:   "Number of samtools threads"
        memory:    "Memory in GB"
        timeout:   "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # Symlink BAM and BAI preserving the original filenames so samtools can
      # locate the index automatically (it derives the index path from the BAM
      # path; renaming would break that lookup).
      ln -s ~{bam} .
      ln -s ~{bai} .
      bam_name=$(basename ~{bam})
      samtools view -h -@ 2 "${bam_name}" ~{chr} \
        | awk '/^@/{print;next} $7=="="' \
        | samtools sort -n -u -@ 2 -m 1G - \
        | samtools fixmate -m -u -@ 2 - - \
        | samtools sort -@ ~{threads} -m 2G -o ~{chr}.fixedmate.bam -
      samtools index ~{chr}.fixedmate.bam
    >>>

    output {
        File fixed_bam = "~{chr}.fixedmate.bam"
        File fixed_bai = "~{chr}.fixedmate.bam.bai"
    }

    runtime {
        cpu:     threads
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}

# Add mate CIGAR (MC) tags to discordant pairs (mates on different chromosomes)
# from one lane BAM.  Extracts only reads where RNEXT != "=" (~1-2% of total),
# so the intermediate file is small and memory requirements are low.
task fixmate_discordant {
    input {
        File        bam
        File        bai
        String      modules = "samtools/1.16.1"
        Int threads = 4
        Int memory  = 16
        Int timeout = 3
    }

    parameter_meta {
        bam:       "BAM file for one lane; only discordant read pairs (RNEXT != '=') are processed, via a full sequential scan"
        bai:       "Index for bam; not used by this task but kept so callers can pass the pair together"
        modules:   "Environment modules to load (samtools required)"
        threads:   "Number of samtools threads"
        memory:    "Memory in GB"
        timeout:   "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # -f 1:  paired
      # -F 12: neither read nor mate unmapped
      # awk:   keep header lines and reads where mate is on a different chromosome
      # No index needed: full-file sequential scan.
      samtools view -h -@ 2 -f 1 -F 12 ~{bam} \
        | awk '/^@/{print;next} $7!="="' \
        | samtools sort -n -u -@ 2 -m 1G - \
        | samtools fixmate -m -u -@ 2 - - \
        | samtools sort -@ ~{threads} -m 2G -o discordant.fixedmate.bam -
      samtools index discordant.fixedmate.bam
    >>>

    output {
        File fixed_bam = "discordant.fixedmate.bam"
        File fixed_bai = "discordant.fixedmate.bam.bai"
    }

    runtime {
        cpu:     threads
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}

# Merge coordinate-sorted BAMs. Used for two different things:
#   sanitize_header=true   fixmate shards from one sample (chromosome slices + discordant)
#   sanitize_header=false  whole alignments from several flowcells/lanes
task merge_bams {
    input {
        Array[File] bams
        Array[File] bais
        Boolean     sanitize_header = true
        String      modules = "samtools/1.16.1"
        Int threads = 8
        Int memory  = 16
        Int timeout = 24
    }

    parameter_meta {
        bams:            "Coordinate-sorted BAMs to merge"
        bais:            "Index files corresponding to each entry in bams (same order)"
        sanitize_header: "When true, strip @PG lines and repair unpaired-read flags. Needed after per-chromosome fixmate, and harmful for vendor alignments because it discards their @PG provenance"
        modules:         "Environment modules to load (samtools required)"
        threads:         "Number of samtools threads"
        memory:          "Memory in GB"
        timeout:         "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      bam_list=(~{sep=" " bams})

      merge_stream() {
        if [ "${#bam_list[@]}" -gt 1 ]; then
          samtools merge -f -c -p -u -@ ~{threads} - "${bam_list[@]}"
        else
          samtools view -h -u "${bam_list[0]}"
        fi
      }

      # @PG stripping: bwa embeds tab-separated fields in CL: causing htsjdk to produce a
      # spurious second ID: field -> SAMFormatException in Redux.
      # Flag fix: per-chromosome fixmate may clear 0x1 while leaving 0x40/0x80 set on
      # discordant-mate reads; htsjdk treats this as a validation error.
      if ~{sanitize_header}; then
        merge_stream \
          | samtools view -h -@ ~{threads} \
          | awk 'BEGIN{OFS="\t"} /^@PG/{next} /^@/{print;next} {
              f=int($2)
              if (and(f,1)==0) { f=f-and(f,64)-and(f,128); $2=f }
              print
            }' \
          | samtools view -@ ~{threads} -O BAM -o merged.bam
      else
        merge_stream | samtools view -@ ~{threads} -O BAM -o merged.bam
      fi

      samtools index -@ ~{threads} merged.bam
    >>>

    output {
        File bam = "merged.bam"
        File bai = "merged.bam.bai"
    }

    runtime {
        cpu:     threads
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}

# Extract a WG tarball (amber/, cobalt/, purple/, pave/, sage/) into a local directory,
# return its absolute path, and work out which tumour sample the outputs belong to.
task extract_wgts {
    input {
        File    tarball
        String? expected_tumor_sample_id
        Int     memory  = 8
        Int     timeout = 1
    }

    parameter_meta {
        tarball:                  "Gzipped tar archive of oncoanalyser WGTS outputs produced by a prior WG run; extracted to a local directory whose absolute path is returned as output_dir"
        expected_tumor_sample_id: "Optional cross-check. If given and it disagrees with the sample ID the tarball is actually named for, the task fails immediately instead of letting SAGE_APPEND discover it hours later"
        memory:                   "Memory in GB"
        timeout:                  "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      mkdir -p wgts_extracted
      tar -xzf ~{tarball} -C wgts_extracted/
      echo "$(pwd)/wgts_extracted" > output_dir.txt

      # Work out the primary tumour sample ID from the PURPLE filenames rather than trusting a
      # hand-typed input. SAGE_APPEND reads <tumor_sample_id>.purple.somatic.vcf.gz, so a wrong
      # ID fails with "invalid path for config: input_vcf = ..." only AFTER REDUX and COBALT
      # have run. The WG step names these files from the tumour BAM's @RG SM tag, which is
      # rarely the same string as group_id -- an easy and expensive mistake to make by hand.
      shopt -s nullglob
      ids=()
      for f in wgts_extracted/purple/*.purple.purity.tsv; do
        ids+=("$(basename "$f" .purple.purity.tsv)")
      done

      if [ "${#ids[@]}" -ne 1 ]; then
        echo "ERROR: expected exactly one *.purple.purity.tsv under purple/, found ${#ids[@]}" >&2
        printf '  %s\n' "${ids[@]+${ids[@]}}" >&2
        exit 1
      fi

      tumor_sample_id="${ids[0]}"
      echo "${tumor_sample_id}" > tumor_sample_id.txt

      # SAGE_APPEND needs this specific file; check now, not an hour in.
      vcf="wgts_extracted/purple/${tumor_sample_id}.purple.somatic.vcf.gz"
      if [ ! -s "${vcf}" ]; then
        echo "ERROR: tarball has no ${vcf}, which SAGE_APPEND requires" >&2
        exit 1
      fi

      expected="~{expected_tumor_sample_id}"
      if [ -n "${expected}" ] && [ "${expected}" != "${tumor_sample_id}" ]; then
        echo "ERROR: tumor_sample_id '${expected}' does not match the WG outputs, which are" >&2
        echo "       named for '${tumor_sample_id}'. Omit tumor_sample_id to use the value" >&2
        echo "       derived from the tarball." >&2
        exit 1
      fi

      echo "primary tumour sample ID: ${tumor_sample_id}"
    >>>

    output {
        String output_dir       = read_string("output_dir.txt")
        String tumor_sample_id  = read_string("tumor_sample_id.txt")
    }

    runtime {
        cpu:     2
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
    }
}

task collect_results {
    input {
        String       outputFileNamePrefix
        File         subject_summary
        File         subject_tarball
        Array[File]  control_summaries
        Array[File]  control_tarballs
        Int memory  = 16
        Int timeout = 10
    }

    parameter_meta {
        outputFileNamePrefix: "Prefix of output file name"
        subject_summary:      "WISP summary TSV for the subject longitudinal sample"
        subject_tarball:      "WISP tarball for the subject longitudinal sample"
        control_summaries:    "WISP summary TSVs for each control sample"
        control_tarballs:     "WISP tarballs for each control sample"
        memory:               "Memory in GB"
        timeout:              "Wall-clock timeout in hours"
    }

    command <<<
        set -euo pipefail

        ## concatenate the summaries
        cat ~{subject_summary} > ~{outputFileNamePrefix}.wisp_summary.tsv
        for f in ~{sep=" " control_summaries}
        do
            tail -n1 $f >> ~{outputFileNamePrefix}.wisp_summary.tsv
        done

        ## retar the tarballs
        mkdir -p wisp
        for tgz in ~{subject_tarball} ~{sep=" " control_tarballs}
        do
            tar -xvf $tgz -C wisp/ --strip-components=1
        done
        tar -cvzf ~{outputFileNamePrefix}.wisp.tar.gz wisp/*
    >>>

    output {
        File wisp_tarballs = "~{outputFileNamePrefix}.wisp.tar.gz"
        File wisp_summary  = "~{outputFileNamePrefix}.wisp_summary.tsv"
    }

    runtime {
        cpu:     2
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
    }
}

task run_wgts {
    input {
        String  group_id
        String  subject_id
        File    tumor_bam
        File    tumor_bai
        String  tumor_sample_id
        File?   normal_bam
        File?   normal_bai
        String? normal_sample_id
        String? normal_platform
        Boolean run_redux
        Boolean nextflow_stub
        Boolean allow_mixed_platforms
        String  sequencing_platform
        String  outdir
        String  ref_data_dir
        String  images_dir
        String  pipeline_dir
        String  nextflow_bin
        String  nextflow_home
        String  oicr_config
        String  qsub_wrapper
        String  modules
        Int memory  = 32
        Int timeout = 24
    }

    parameter_meta {
        group_id:            "Sample group identifier; used as output subdirectory and samplesheet group_id"
        subject_id:          "Subject/patient identifier"
        tumor_bam:           "Primary tumour BAM (merged if there were several inputs)"
        tumor_bai:           "Index for tumor_bam"
        tumor_sample_id:     "Tumour sample ID used in the samplesheet and as the output file prefix"
        normal_bam:          "Matched normal BAM; omit for tumour-only WGTS"
        normal_bai:          "Index for normal_bam; required when normal_bam is provided"
        normal_sample_id:    "Normal sample ID; required when normal_bam is provided"
        normal_platform:     "Sequencing platform of the normal, used only to refuse a mixed-platform run. oncoanalyser takes one --sequencing_platform per run, and this run processes both samples"
        run_redux:           "When true REDUX processes the alignments; when false the inputs are declared bam_redux with generate_redux_tsvs_only so REDUX only regenerates its TSVs"
        allow_mixed_platforms: "When false, refuse to launch if the samples in this single pipeline run span more than one sequencing platform"
        nextflow_stub:       "Run oncoanalyser with -stub --create_stub_placeholders: placeholder outputs, no real compute"
        sequencing_platform: "Value for --sequencing_platform: illumina, sbx or ultima"
        outdir:              "Output directory; the pipeline writes to outdir/group_id/"
        ref_data_dir:        "Staged HMF reference data directory, used for --igenomes_base, --hmf_genomes_base and --ref_data_hmf_data_path"
        images_dir:          "Singularity image cache directory (NXF_SINGULARITY_CACHEDIR)"
        pipeline_dir:        "oncoanalyser checkout containing main.nf"
        nextflow_bin:        "Path to the Nextflow executable; must be 25.10.0 or newer for this pipeline"
        nextflow_home:       "NXF_HOME holding the pre-cached nf-schema plugin, so the run works with NXF_OFFLINE=true"
        oicr_config:         "OICR overlay config passed with -c (SGE settings, container overrides, unused reference paths nulled out)"
        qsub_wrapper:        "qsub shim copied onto PATH; strips h_rss/mem_free from Nextflow's generated .command.run so SGE's cgroup does not kill JVM tasks"
        modules:             "Environment modules to load"
        memory:              "Memory in GB for the Cromwell/SGE task (hosts the Nextflow driver JVM only)"
        timeout:             "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      mkdir -p ~{outdir}
      abs_outdir=$(readlink -f ~{outdir})

      WORKDIR=$(pwd)

      # Optional inputs render as the empty string when absent, so capture them in shell
      # variables and branch on those rather than testing the interpolation directly.
      normal_bam="~{normal_bam}"
      normal_bai="~{normal_bai}"
      normal_sample_id="~{normal_sample_id}"
      normal_platform="~{normal_platform}"

      # oncoanalyser takes ONE --sequencing_platform per run, and this run processes both
      # samples below. Refuse rather than silently apply the wrong error model to one of them.
      if [ -n "${normal_platform}" ] && [ "${normal_platform}" != "~{sequencing_platform}" ]; then
        {
          echo "mixed sequencing platforms within a SINGLE oncoanalyser run:"
          echo "  this run uses --sequencing_platform ~{sequencing_platform}"
          echo "  but the normal sample is ${normal_platform}"
          echo "oncoanalyser has no per-sample platform setting, so one of the two samples"
          echo "would be analysed with the wrong error model."
        } >&2
        if ~{allow_mixed_platforms}; then
          echo "WARNING: proceeding anyway because allow_mixed_platforms=true." >&2
          echo "         Results for the ${normal_platform} sample are not trustworthy." >&2
        else
          {
            echo "REFUSING TO LAUNCH. Options:"
            echo "  - omit normal_alignments so this run is single-platform (loses LOH)"
            echo "  - use samples sequenced on one platform"
            echo "  - set allow_mixed_platforms=true for a deliberate experiment"
          } >&2
          exit 1
        fi
      fi

      ln -s "~{tumor_bam}" "${WORKDIR}/~{tumor_sample_id}.bam"
      ln -s "~{tumor_bai}" "${WORKDIR}/~{tumor_sample_id}.bam.bai"
      if [ -n "${normal_bam}" ]; then
        ln -s "${normal_bam}" "${WORKDIR}/${normal_sample_id}.bam"
        ln -s "${normal_bai}" "${WORKDIR}/${normal_sample_id}.bam.bai"
      fi

      # run_redux=false means "these are already REDUX alignments". REDUX still runs, but
      # with -bqr_jitter_msi_only, regenerating only the TSVs and symlinking the input as
      # the redux BAM. --skip_msi_jitter was removed in 3.x and must not be passed.
      if ~{run_redux}; then
        filetype="bam"
        info=""
      else
        filetype="bam_redux"
        info="generate_redux_tsvs_only"
      fi

      # Index rows are omitted: oncoanalyser finds <alignment>.bai by adjacency, and the
      # symlinks above put it there.
      {
        echo "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath"
        if [ -n "${normal_bam}" ]; then
          echo "~{group_id},~{subject_id},${normal_sample_id},normal,dna,${filetype},${info},${WORKDIR}/${normal_sample_id}.bam"
        fi
        echo "~{group_id},~{subject_id},~{tumor_sample_id},tumor,dna,${filetype},${info},${WORKDIR}/~{tumor_sample_id}.bam"
      } > samplesheet.csv

      # -stub makes every process emit placeholder outputs, so the wrapper can be exercised
      # end to end in minutes. --create_stub_placeholders also fabricates the genome files
      # the pipeline would otherwise insist on finding.
      stub_args=()
      if ~{nextflow_stub}; then
        stub_args=(-stub --create_stub_placeholders)
      fi

      export NXF_OFFLINE=true
      export NXF_OPTS="-Xms512m -Xmx8g"
      export NXF_SINGULARITY_CACHEDIR=~{images_dir}
      export NXF_HOME=~{nextflow_home}
      # A loaded oncoanalyser 2.x module points these at its read-only tree.
      unset NXF_DIST NXF_LAUNCHER NXF_PLUGINS_DIR || true

      bin_dir="$(pwd)/bin"
      mkdir -p "${bin_dir}"
      cp ~{qsub_wrapper} "${bin_dir}/qsub"
      chmod +x "${bin_dir}/qsub"
      export PATH="${bin_dir}:$PATH"

      ~{nextflow_bin} run ~{pipeline_dir}/main.nf \
          --mode wgts \
          --sequencing_platform ~{sequencing_platform} \
          --input samplesheet.csv \
          --outdir ~{outdir} \
          --genome GRCh38_hmf \
          --processes_manual redux,amber,cobalt,sage,pave,purple \
          --igenomes_base ~{ref_data_dir} \
          --hmf_genomes_base ~{ref_data_dir} \
          --ref_data_hmf_data_path ~{ref_data_dir} \
          -profile singularity \
          -c ~{oicr_config} \
          "${stub_args[@]}" \
          -resume

      echo "${abs_outdir}/~{group_id}" > output_dir.txt

      # alignments/ is excluded on purpose: it holds full REDUX BAMs. A chained PE run reads
      # the primary outputs from disk, not from this tarball.
      #
      # Collect only the directories that exist, so a run that legitimately produces fewer
      # of them (e.g. tumour-only, where nothing germline is written) does not fail at the
      # very last step after hours of compute.
      tar_dirs=()
      for d in amber cobalt purple pave sage; do
        if [ -d "${abs_outdir}/~{group_id}/${d}" ]; then
          tar_dirs+=("${d}/")
        else
          echo "note: ${d}/ not present in the output, omitting from the tarball" >&2
        fi
      done
      if [ "${#tar_dirs[@]}" -eq 0 ]; then
        echo "ERROR: none of amber/ cobalt/ purple/ pave/ sage/ were produced" >&2
        exit 1
      fi
      tar -czf ~{outdir}.wgts.tar.gz \
          -C "${abs_outdir}/~{group_id}" \
          "${tar_dirs[@]}"

      tar -czf ~{outdir}.pipeline_info.tar.gz \
          -C "${abs_outdir}" \
          pipeline_info/
    >>>

    output {
        String output_dir             = read_string("output_dir.txt")
        File   wgts_tarball           = "~{outdir}.wgts.tar.gz"
        File   pipeline_info_tarball  = "~{outdir}.pipeline_info.tar.gz"
    }

    runtime {
        cpu:     2
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}

task run_purity_estimate {
    input {
        String  group_id
        String  subject_id
        String  tumor_sample_id
        File?   normal_bam
        File?   normal_bai
        String? normal_sample_id
        String? normal_platform
        File    longitudinal_bam
        File    longitudinal_bai
        String  longitudinal_sample_id
        Boolean run_redux
        Boolean nextflow_stub
        Boolean allow_mixed_platforms
        String  sequencing_platform
        String  wgts_outdir
        String  outdir
        String  ref_data_dir
        String  images_dir
        String  pipeline_dir
        String  nextflow_bin
        String  nextflow_home
        String  oicr_config
        String  qsub_wrapper
        String  modules
        Int memory  = 32
        Int timeout = 10
    }

    parameter_meta {
        group_id:               "Sample group identifier; used as output subdirectory and samplesheet group_id"
        subject_id:             "Subject/patient identifier; prepended to the wisp summary row"
        tumor_sample_id:        "Primary tumour sample ID; must match the purple/amber output filenames in wgts_outdir"
        normal_bam:             "Primary normal BAM. Required to use LOH: oncoanalyser rejects an amber_dir input unless a primary normal alignment is also given. Omit for COBALT-only estimation"
        normal_bai:             "Index for normal_bam; required when normal_bam is provided"
        normal_sample_id:       "Normal sample ID; required when normal_bam is provided"
        normal_platform:        "Sequencing platform of the normal, used only to refuse a mixed-platform run. oncoanalyser takes one --sequencing_platform per run, and REDUX here processes both the longitudinal sample and this normal"
        longitudinal_bam:       "Longitudinal ctDNA BAM (merged if there were several inputs)"
        longitudinal_bai:       "Index for longitudinal_bam"
        longitudinal_sample_id: "Longitudinal sample ID used in the samplesheet"
        run_redux:              "When true REDUX processes the alignments; when false they are declared bam_redux with generate_redux_tsvs_only so REDUX only regenerates its TSVs"
        allow_mixed_platforms:  "When false, refuse to launch if the samples in this single pipeline run span more than one sequencing platform"
        nextflow_stub:          "Run oncoanalyser with -stub --create_stub_placeholders: placeholder outputs, no real compute"
        sequencing_platform:    "Value for --sequencing_platform: illumina, sbx or ultima"
        wgts_outdir:            "Path to the WG output directory (purple/ and amber/ subdirectories expected)"
        outdir:                 "Output directory; the pipeline writes to outdir/group_id/"
        ref_data_dir:           "Staged HMF reference data directory, used for --igenomes_base, --hmf_genomes_base and --ref_data_hmf_data_path"
        images_dir:             "Singularity image cache directory (NXF_SINGULARITY_CACHEDIR)"
        pipeline_dir:           "oncoanalyser checkout containing main.nf"
        nextflow_bin:           "Path to the Nextflow executable; must be 25.10.0 or newer for this pipeline"
        nextflow_home:          "NXF_HOME holding the pre-cached nf-schema plugin, so the run works with NXF_OFFLINE=true"
        oicr_config:            "OICR overlay config passed with -c (SGE settings, container overrides, unused reference paths nulled out)"
        qsub_wrapper:           "qsub shim copied onto PATH; strips h_rss/mem_free from Nextflow's generated .command.run so SGE's cgroup does not kill JVM tasks"
        modules:                "Environment modules to load"
        memory:                 "Memory in GB for the Cromwell/SGE task (hosts the Nextflow driver JVM only)"
        timeout:                "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      WORKDIR=$(pwd)

      # Optional inputs render as the empty string when absent, so capture them in shell
      # variables and branch on those rather than testing the interpolation directly.
      normal_bam="~{normal_bam}"
      normal_bai="~{normal_bai}"
      normal_sample_id="~{normal_sample_id}"
      normal_platform="~{normal_platform}"

      # oncoanalyser takes ONE --sequencing_platform per run, and this run processes both
      # samples below. Refuse rather than silently apply the wrong error model to one of them.
      if [ -n "${normal_platform}" ] && [ "${normal_platform}" != "~{sequencing_platform}" ]; then
        {
          echo "mixed sequencing platforms within a SINGLE oncoanalyser run:"
          echo "  this run uses --sequencing_platform ~{sequencing_platform}"
          echo "  but the normal sample is ${normal_platform}"
          echo "oncoanalyser has no per-sample platform setting, so one of the two samples"
          echo "would be analysed with the wrong error model."
        } >&2
        if ~{allow_mixed_platforms}; then
          echo "WARNING: proceeding anyway because allow_mixed_platforms=true." >&2
          echo "         Results for the ${normal_platform} sample are not trustworthy." >&2
        else
          {
            echo "REFUSING TO LAUNCH. Options:"
            echo "  - omit normal_alignments so this run is single-platform (loses LOH)"
            echo "  - use samples sequenced on one platform"
            echo "  - set allow_mixed_platforms=true for a deliberate experiment"
          } >&2
          exit 1
        fi
      fi

      ln -s "~{longitudinal_bam}" "${WORKDIR}/~{longitudinal_sample_id}.bam"
      ln -s "~{longitudinal_bai}" "${WORKDIR}/~{longitudinal_sample_id}.bam.bai"
      if [ -n "${normal_bam}" ]; then
        ln -s "${normal_bam}" "${WORKDIR}/${normal_sample_id}.bam"
        ln -s "${normal_bai}" "${WORKDIR}/${normal_sample_id}.bam.bai"
      fi

      if ~{run_redux}; then
        filetype="bam"
        redux_info=""
      else
        filetype="bam_redux"
        redux_info="generate_redux_tsvs_only"
      fi

      # The longitudinal sample carries longitudinal_sample, plus the redux hint when
      # applicable; multiple info fields are separated by ';'.
      long_info="longitudinal_sample"
      if [ -n "${redux_info}" ]; then
        long_info="longitudinal_sample;${redux_info}"
      fi

      # AMBER only earns its keep when a primary normal is available for LOH, and
      # oncoanalyser errors on an amber_dir input without one.
      if [ -n "${normal_bam}" ]; then
        processes="redux,amber,cobalt,sage_append,wisp"
      else
        processes="redux,cobalt,sage_append,wisp"
      fi

      # amber_dir / cobalt_dir / sage_append_dir must never appear on the longitudinal
      # sample: oncoanalyser treats that as a fatal input clash.
      {
        echo "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath"
        echo "~{group_id},~{subject_id},~{tumor_sample_id},tumor,dna,purple_dir,,~{wgts_outdir}/purple/"
        if [ -n "${normal_bam}" ]; then
          echo "~{group_id},~{subject_id},~{tumor_sample_id},tumor,dna,amber_dir,,~{wgts_outdir}/amber/"
          echo "~{group_id},~{subject_id},${normal_sample_id},normal,dna,${filetype},${redux_info},${WORKDIR}/${normal_sample_id}.bam"
        fi
        echo "~{group_id},~{subject_id},~{longitudinal_sample_id},tumor,dna,${filetype},${long_info},${WORKDIR}/~{longitudinal_sample_id}.bam"
      } > samplesheet_purity.csv

      # -stub makes every process emit placeholder outputs, so the wrapper can be exercised
      # end to end in minutes. --create_stub_placeholders also fabricates the genome files
      # the pipeline would otherwise insist on finding.
      stub_args=()
      if ~{nextflow_stub}; then
        stub_args=(-stub --create_stub_placeholders)
      fi

      export NXF_OFFLINE=true
      export NXF_OPTS="-Xms512m -Xmx8g"
      export NXF_SINGULARITY_CACHEDIR=~{images_dir}
      export NXF_HOME=~{nextflow_home}
      unset NXF_DIST NXF_LAUNCHER NXF_PLUGINS_DIR || true

      bin_dir="$(pwd)/bin"
      mkdir -p "${bin_dir}"
      cp ~{qsub_wrapper} "${bin_dir}/qsub"
      chmod +x "${bin_dir}/qsub"
      export PATH="${bin_dir}:$PATH"

      ~{nextflow_bin} run ~{pipeline_dir}/main.nf \
          --mode purity_estimate \
          --purity_estimate_mode wgts \
          --sequencing_platform ~{sequencing_platform} \
          --input samplesheet_purity.csv \
          --outdir ~{outdir} \
          --genome GRCh38_hmf \
          --processes_manual ${processes} \
          --igenomes_base ~{ref_data_dir} \
          --hmf_genomes_base ~{ref_data_dir} \
          --ref_data_hmf_data_path ~{ref_data_dir} \
          -profile singularity \
          -c ~{oicr_config} \
          "${stub_args[@]}" \
          -resume

      abs_outdir=$(readlink -f ~{outdir})
      wisp_dir="~{outdir}/~{group_id}/wisp"
      tar -czf ~{outdir}.wisp.tar.gz -C "$(dirname "$wisp_dir")" wisp/

      # prepend subject_id as first column in summary
      head -n1 "$wisp_dir"/*.wisp.summary.tsv | sed 's/^/sample_id\t/' >  ~{outdir}.wisp_summary.tsv
      tail -n1 "$wisp_dir"/*.wisp.summary.tsv | sed 's/^/~{subject_id}\t/' >> ~{outdir}.wisp_summary.tsv

      tar -czf ~{outdir}.pipeline_info.tar.gz \
          -C "${abs_outdir}" \
          pipeline_info/
    >>>

    output {
        File wisp_tarball            = "~{outdir}.wisp.tar.gz"
        File wisp_summary            = "~{outdir}.wisp_summary.tsv"
        File pipeline_info_tarball   = "~{outdir}.pipeline_info.tar.gz"
    }

    runtime {
        cpu:     2
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: modules
    }
}
