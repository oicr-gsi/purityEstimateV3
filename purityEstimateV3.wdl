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
#
#   * --skip_msi_jitter is gone. run_redux=false now means filetype bam_redux plus info
#     generate_redux_tsvs_only, with redux STILL in --processes_manual, so REDUX regenerates
#     the TSVs instead of us faking them with touch.
#   * --sequencing_platform must be passed for SBX / Ultima.
#   * CRAM is converted to BAM against the reference it was ENCODED with. For Ultima that is
#     hs38DH, which differs from oncoanalyser's masked GRCh38 on chr9, chr13, chr16, chr21 and
#     chrX, so decoding against the wrong one silently corrupts bases.
#   * SAGE_APPEND is directly selectable, so 'orange' is no longer needed to trigger it.
#
# Platform is resolved per sample, not per run: fixmate applies only to paired reads, so
# running it on a single-end sample would drop every record. Each head job is passed the
# platform of the sample it processes.
#
# Upstream limitation behind allow_mixed_platforms: oncoanalyser has one --sequencing_platform
# per run, so a PE run that also processes a normal of a different platform cannot be expressed
# correctly.
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
        Boolean       doFixmate = true        # see parameter_meta
        Boolean       run_control = false     # also estimate purity for each control (PE/WG_PE only)
        Array[Pair[String, String]]? controls # PE mode: (control_id, bam_path) pairs; BAI at bam_path+".bai"
        String        outputFileNamePrefix
        Boolean       include_germline_outputs = false  # see parameter_meta
        Boolean       allow_mixed_platforms = false  # production guardrail; see parameter_meta
        Boolean       use_primary_filters = true  # PE mode; see parameter_meta
        Boolean       nextflow_stub = false   # run oncoanalyser with -stub; see parameter_meta
        String        scheduler = ""          # empty = detect; see parameter_meta
        String?       slurm_partition         # required when scheduler = "slurm"
        String?       slurm_account           # see parameter_meta
        Array[String]? singularity_binds      # null = keep the binds the module's overlay sets
        Array[String]? nextflow_config        # null = the module's config for `scheduler`
        String        modules = "java/17 singularity/3.9.4 samtools/1.16.1 oncoanalyser/3.0.0-rc.3 oncoanalyser-data/3.0.0"
    }

    parameter_meta {
        mode:                   "Run mode: WG=WGTS only; PE=purity estimation against a pre-existing WG tarball; WG_PE=run both sequentially"
        group_id:               "Sample group identifier used as the output subdirectory name and samplesheet group_id"
        subject_id:             "Subject/patient identifier"
        tumor_alignments:       "Primary tumour BAM/CRAM(s) with indices; required for WG and WG_PE mode. Multiple entries (e.g. per flowcell) are merged"
        normal_alignments:      "Matched normal BAM/CRAM(s) with indices. MANDATORY for WG and WG_PE, with no override: without it the primary is called tumour-only, germline variants are not subtracted from the somatic call set, and the MRD result is a false positive. Not used by the PE step, which does not pass the normal to WISP; supplying it in PE mode wastes staging effort"
        longitudinal_alignments: "ctDNA BAM/CRAM(s) with indices; required for PE and WG_PE mode. Multiple entries are merged"
        wgts_tarball:           "Tarball produced by a prior WG run (amber/, cobalt/, purple/, pave/, sage/); required for PE mode"
        tumor_sample_id:        "Primary tumour sample ID. Normally leave unset: in WG modes it is read from the tumour @RG SM tag, and in PE mode it is derived from the WG tarball's purple/ filenames. If given, it overrides the @RG SM tag in WG modes and is cross-checked against the tarball in PE mode"
        sequencing_platform:    "Sequencing platform passed to --sequencing_platform: illumina, sbx or ultima. Leave unset to detect it from the @RG PL tag"
        run_redux:              "When true, REDUX processes the alignments normally. When false, inputs are treated as already REDUX-processed and REDUX only regenerates its TSVs (-bqr_jitter_msi_only). Does not affect control BAMs"
        doFixmate:              "Whether to add mate CIGAR (MC) tags before REDUX. Only consulted when run_redux is true and the sample is Illumina, since single-end reads have no mates to fix. Leave true for alignments that lack MC tags; set false when the aligner already wrote them, as bwa-mem2 does, to skip the per-chromosome fixmate scatter entirely. The inputs are still merged, because oncoanalyser's samplesheet accepts one file per sample and filetype, FASTQ excepted; REDUX itself would take a list"
        run_control:            "When true, also run purity estimation for each control BAM in the controls array (PE/WG_PE mode only)"
        controls:               "PE mode: array of (control_id, bam_path) pairs; BAI assumed at bam_path+'.bai'. Controls are BAM only and always run with run_redux=false; used only when run_control=true"
        outputFileNamePrefix:   "Output directory prefix; the pipeline writes to outputFileNamePrefix/group_id/"
        include_germline_outputs: "Whether to keep germline variant calls in the WG archive. Defaults false: germline calls play no part in MRD, since WISP reads the PURPLE somatic VCF and upstream purity_estimate.nf disables germline calling itself. Note they are still GENERATED -- oncoanalyser hardcodes germline calling on and it cannot be disabled from configuration -- so this drops sage/germline, pave/germline and the PURPLE germline files from the archive rather than skipping the work"
        allow_mixed_platforms:  "Guardrail. oncoanalyser applies ONE --sequencing_platform per pipeline run, so two samples of different platforms in the same run means one of them gets the wrong error model. Left false (the production default) such a combination fails before Nextflow starts. Set true only for deliberate experiments: the run proceeds with a loud warning"
        use_primary_filters:    "PE and WG_PE mode: whether the plasma stage works from the prefiltered primary VCF rather than the full call set. The WG stage always writes both, so this only chooses between them. A WG tarball produced before pre-filtering existed contains no prefiltered VCF; the run then falls back to the full set with a warning rather than failing"
        scheduler:              "Which batch scheduler Nextflow submits to, sge or slurm. Leave empty and validate_inputs decides from the submit command present on the cluster, which is what makes one set of inputs portable between sites; the value it resolved is reported in that task's output and used by the head jobs. Set it only to override that. It selects whether the submit wrapper is placed on PATH: the wrapper exists to rewrite the h_rss/mem_free directives the SGE executor embeds directly in .command.run, where configuration cannot reach them. The Slurm executor emits --mem and --cpus-per-task from the memory and cpus directives, so it needs no wrapper"
        slurm_partition:        "Partition Nextflow submits its own jobs to, required when scheduler is slurm. The overlay the oncoanalyser module ships names a queue for its own scheduler, which does not exist elsewhere"
        slurm_account:          "Accounting group for the jobs Nextflow submits, when the site requires one. Also clears the resource request the module's overlay writes in the other scheduler's syntax, which sbatch would reject, so leave it null only where no account is needed"
        singularity_binds:      "Filesystem paths bound into every container, replacing the bind the module's overlay sets. Leave null where the containers can already reach the reference data and the working directory, which is the case when the run shares a filesystem with the site the module was built for"
        nextflow_config:        "Config overlays passed to nextflow, each with its own -c, in order. Leave null on SGE to use the overlay the oncoanalyser module ships. A site whose scheduler differs must supply its own, because that overlay sets the executor; split it so the settings that hold everywhere (genome paths, container overrides, per-process resources) stay in one file and only the executor and filesystem binds are per-site"
        nextflow_stub:          "When true, oncoanalyser runs with -stub --create_stub_placeholders: every process writes placeholder outputs instead of doing real work. This exercises the whole wrapper (samplesheet, samplesheet validation, output layout, tarring, Vidarr outputs) in minutes. Real input alignments are still required, but they can be tiny, because the pipeline never reads them. Not a Cromwell dry run"
        modules:                "Environment modules to load. The oncoanalyser module supplies the pipeline checkout, the container image cache, NEXTFLOW_HOME, the scheduler submit wrapper, the site config overlay and the Nextflow launcher; the oncoanalyser-data module supplies the reference bundle. The resource paths below read variables exported by both, so the module versions here and those paths must stay in step"
    }

    # Resource paths. String (not File) so Cromwell reads them off the shared filesystem
    # instead of localizing gigabytes per task.
    #
    # The first five are the LITERAL TEXT of variables exported by the oncoanalyser module,
    # not their values: a workflow-level declaration is evaluated before any task runs, so
    # the variable cannot be expanded here. They are interpolated into the command blocks and
    # the shell expands them there, after `modules` has been loaded. Under `set -euo pipefail`
    # an unset one fails the task loudly, which is the behaviour we want.
    String images_dir     = "$IMAGES_DIR"
    String pipeline_dir   = "$ONCOANALYSER_FOLDER"
    String nextflow_home  = "$NEXTFLOW_HOME"
    # The module ships an overlay for its own scheduler. The head jobs generate the overrides
    # a different scheduler needs from slurm_partition, slurm_account and singularity_binds,
    # and apply them after everything here, so this list carries only the settings that are
    # independent of the scheduler.
    Array[String] site_configs = select_first([nextflow_config, ["$ONCOANALYSER_OICR_CONFIG"]])
    String submit_wrapper = "$QSUB_WRAPPER"

    # Nextflow comes from PATH; the pipeline module supplies it and pins the version. The
    # nf-schema plugin requires a newer Nextflow than the pipeline manifest declares.
    String nextflow_bin   = "nextflow"

    # Reference data comes from the SEPARATE oncoanalyser-data module, which is versioned on
    # its own cycle. Use the variables it exports rather than deriving paths from its root, so
    # a re-layout inside the module does not break us. cram_to_bam needs the CRAM reference but
    # not the pipeline, so it loads only the data module (see its `modules` default) -- a data
    # module sets variables and adds no libraries, so it is cheap and safe to load anywhere.
    String ref_data_dir   = "$REFERENCE_FILES_DIR"
    String cram_reference = "$VENDOR_GENOME_HS38DH"

    # GRCh38 chromosomes processed individually to keep fixmate memory small.
    Array[String] chromosomes = [
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
        "chr21","chr22","chrX","chrY","chrM"
    ]

    # Fail fast and legibly on a mode/input mismatch. Without this the first symptom is an
    # opaque select_first error from a downstream call.
    call validate_inputs {
        input:
            mode                = mode,
            has_tumor           = defined(tumor_alignments),
            has_normal          = defined(normal_alignments),
            has_longitudinal    = defined(longitudinal_alignments),
            has_wgts_tarball    = defined(wgts_tarball),
            has_tumor_sample_id = defined(tumor_sample_id),
            run_control         = run_control,
            has_controls        = defined(controls),
            scheduler           = scheduler,
            slurm_partition     = select_first([slurm_partition, ""]),
            slurm_account       = select_first([slurm_account, ""]),
            singularity_binds   = select_first([singularity_binds, []]),
            nextflow_config     = select_first([nextflow_config, []])
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
    # Getting this wrong is silent and destructive: fixmate keeps only paired records, so on
    # a single-end sample every record is dropped and the merged BAM ends up header-only,
    # which fails obscurely much later. An explicit sequencing_platform input overrides the
    # detection, for data with a missing or wrong PL tag.
    if (defined(tumor_alignments)) {
        String tumor_platform  = select_first([sequencing_platform, tumor_info.platform])
        # Only probe when the answer can change the decision: doFixmate=true always
        # fixmates, and a non-Illumina sample never does.
        if (run_redux && !doFixmate && tumor_platform == "illumina") {
            call probe_mc_tags as probe_tumor {
                input: bam = select_first([tumor_alignments])[0].aln
            }
        }
        # Correct only the dangerous direction. Skipping fixmate is honoured when the tags
        # are really there; when they are not, fixmate runs anyway rather than letting
        # REDUX mis-mark duplicates.
        Boolean tumor_fixmate = run_redux && tumor_platform == "illumina"
                                 && (doFixmate || !select_first([probe_tumor.has_mc_tags, true]))
    }
    if (defined(normal_alignments)) {
        String normal_platform = select_first([sequencing_platform, normal_info.platform])
        # Only probe when the answer can change the decision: doFixmate=true always
        # fixmates, and a non-Illumina sample never does.
        if (run_redux && !doFixmate && normal_platform == "illumina") {
            call probe_mc_tags as probe_normal {
                input: bam = select_first([normal_alignments])[0].aln
            }
        }
        # Correct only the dangerous direction. Skipping fixmate is honoured when the tags
        # are really there; when they are not, fixmate runs anyway rather than letting
        # REDUX mis-mark duplicates.
        Boolean normal_fixmate = run_redux && normal_platform == "illumina"
                                 && (doFixmate || !select_first([probe_normal.has_mc_tags, true]))
    }
    if (defined(longitudinal_alignments)) {
        String longitudinal_platform = select_first([sequencing_platform, longitudinal_info.platform])
        # Only probe when the answer can change the decision: doFixmate=true always
        # fixmates, and a non-Illumina sample never does.
        if (run_redux && !doFixmate && longitudinal_platform == "illumina") {
            call probe_mc_tags as probe_longitudinal {
                input: bam = select_first([longitudinal_alignments])[0].aln
            }
        }
        # Correct only the dangerous direction. Skipping fixmate is honoured when the tags
        # are really there; when they are not, fixmate runs anyway rather than letting
        # REDUX mis-mark duplicates.
        Boolean longitudinal_fixmate = run_redux && longitudinal_platform == "illumina"
                                 && (doFixmate || !select_first([probe_longitudinal.has_mc_tags, true]))
    }

    # The tumour sample ID is resolved per mode, NOT here. A workflow-level declaration is
    # evaluated eagerly by Cromwell even when only used inside a conditional, so
    # select_first([tumor_sample_id, tumor_info.sample_id]) at this scope blows up in PE mode
    # where both are legitimately absent. WG modes derive it in the WG block below (@RG SM or
    # an explicit override); PE derives it from the tarball in extract_wgts.

    # ---------------------------------------------------------------------------------
    # Stage each sample: CRAM -> BAM, then optional fixmate, then merge when needed.
    #
    # @PG sanitisation is tied to run_redux rather than to fixmate. htsjdk rejects a header
    # whose @PG CL: field contains tabs, which an aligner can write, so anything going on to
    # REDUX needs it -- whether or not fixmate ran. Alignments that bypass REDUX keep their
    # @PG provenance.
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
                    repair_flags = true,
                    strip_bad_pg = true
            }
        }

        # Nothing to merge for a single already-fixmate-free alignment; use it as is.
        if (!select_first([tumor_fixmate]) && (run_redux || length(tumor_staged_bam) > 1)) {
            call merge_bams as merge_tumor_plain {
                input:
                    bams = tumor_staged_bam,
                    bais = tumor_staged_bai,
                    strip_bad_pg   = run_redux
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
                    repair_flags = true,
                    strip_bad_pg = true
            }
        }

        if (!select_first([normal_fixmate]) && (run_redux || length(normal_staged_bam) > 1)) {
            call merge_bams as merge_normal_plain {
                input:
                    bams = normal_staged_bam,
                    bais = normal_staged_bai,
                    strip_bad_pg   = run_redux
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
                    repair_flags = true,
                    strip_bad_pg = true
            }
        }

        if (!select_first([longitudinal_fixmate]) && (run_redux || length(longitudinal_staged_bam) > 1)) {
            call merge_bams as merge_longitudinal_plain {
                input:
                    bams = longitudinal_staged_bam,
                    bais = longitudinal_staged_bai,
                    strip_bad_pg   = run_redux
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
                scheduler           = validate_inputs.scheduler_used,
                slurm_partition     = select_first([slurm_partition, ""]),
                slurm_account       = select_first([slurm_account, ""]),
                singularity_binds   = select_first([singularity_binds, []]),
                nextflow_config     = site_configs,
                submit_wrapper        = submit_wrapper,
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
        String extracted_purple_dir = extract_wgts.output_dir + "/purple"
    }

    # The primary tumour sample id, read off the outputs in both directions: PE from the
    # tarball's purple/ filenames, WG and WG_PE from the tumour @RG SM tag.
    String primary_sample_id = select_first([extract_wgts.tumor_sample_id, wg_tumor_sample_id])

    # WG modes: filter the primary's somatic VCF, then archive. Both filter sets are applied
    # here, and the prefiltered VCF goes into the archive beside the untouched original --
    # the WG deliverable is what a clinical report is written from, and what decides whether
    # a plasma sample is worth taking. Packing is a separate task purely so it runs AFTER
    # filtering and can therefore carry both VCFs.
    if (mode == "WG" || mode == "WG_PE") {
        call pre_filtering {
            input:
                purple_dir           = select_first([run_wgts.output_dir]) + "/purple",
                tumor_sample_id      = primary_sample_id,
                outputFileNamePrefix = outputFileNamePrefix
        }

        call pack_wgts {
            input:
                wgts_dir                 = select_first([run_wgts.output_dir]),
                purple_dir               = pre_filtering.purple_dir_out,
                outputFileNamePrefix     = outputFileNamePrefix,
                include_germline_outputs = include_germline_outputs
        }
    }

    if (mode == "PE" || mode == "WG_PE") {
        # WG_PE reads the primary from what pre_filtering wrote; PE from the extracted
        # tarball.
        String pe_purple_dir = select_first([pre_filtering.purple_dir_out, extracted_purple_dir])

        call run_purity_estimate as subject_purity {
            input:
                group_id               = group_id,
                subject_id             = subject_id,
                tumor_sample_id        = primary_sample_id,
                longitudinal_bam       = select_first([longitudinal_bam]),
                longitudinal_bai       = select_first([longitudinal_bai]),
                longitudinal_sample_id = select_first([longitudinal_info.sample_id]),
                run_redux              = run_redux,
                nextflow_stub          = nextflow_stub,
                sequencing_platform    = select_first([longitudinal_platform]),
                primary_purple_dir     = pe_purple_dir,
                use_primary_filters    = use_primary_filters,
                outdir                 = outputFileNamePrefix,
                ref_data_dir           = ref_data_dir,
                images_dir             = images_dir,
                pipeline_dir           = pipeline_dir,
                nextflow_bin           = nextflow_bin,
                nextflow_home          = nextflow_home,
                scheduler              = validate_inputs.scheduler_used,
                slurm_partition     = select_first([slurm_partition, ""]),
                slurm_account       = select_first([slurm_account, ""]),
                singularity_binds   = select_first([singularity_binds, []]),
                nextflow_config        = site_configs,
                submit_wrapper           = submit_wrapper,
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
                        tumor_sample_id        = primary_sample_id,
                        longitudinal_bam       = control.right,
                        longitudinal_bai       = control.right + ".bai",
                        longitudinal_sample_id = control_info.sample_id,
                        run_redux              = false,
                        nextflow_stub          = nextflow_stub,
                        sequencing_platform    = select_first([sequencing_platform, control_info.platform]),
                        primary_purple_dir     = pe_purple_dir,
                use_primary_filters    = use_primary_filters,
                        outdir                 = outputFileNamePrefix,
                        ref_data_dir           = ref_data_dir,
                        images_dir             = images_dir,
                        pipeline_dir           = pipeline_dir,
                        nextflow_bin           = nextflow_bin,
                        nextflow_home          = nextflow_home,
                        scheduler              = validate_inputs.scheduler_used,
                        slurm_partition     = select_first([slurm_partition, ""]),
                        slurm_account       = select_first([slurm_account, ""]),
                        singularity_binds   = select_first([singularity_binds, []]),
                        nextflow_config        = site_configs,
                        submit_wrapper           = submit_wrapper,
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
        description: "Runs HMF oncoanalyser 3.0.0-rc.3 to estimate tumour purity in longitudinal ctDNA samples, for Illumina or Ultima Genomics data. In WG mode it runs WGTS (REDUX, AMBER, COBALT, SAGE, PAVE, PURPLE) on a primary tumour with an optional matched normal and produces a tarball of the results. In PE mode it runs WISP against a pre-existing WG tarball to report the ctDNA fraction of a longitudinal sample. WG_PE does both in sequence. BAM and CRAM are both accepted.\n\n![purityEstimateV3 workflow](docs/purityEstimateV3.svg)\n\nIn the chart the two head-job boxes are where Cromwell stops and Nextflow starts: `run_wgts` and `run_purity_estimate` are each a SINGLE Cromwell task that runs `nextflow run`, and every process inside them is submitted to the cluster by Nextflow itself. A run directory therefore holds far fewer `call-` directories than there are tools. Diagram source is Graphviz, in docs/.\n\n### Valid input combinations\n\n| mode | tumor_alignments | normal_alignments | longitudinal_alignments | wgts_tarball |\n|---|---|---|---|---|\n| WG | required | **required** | - | - |\n| PE | - | not used | required | required |\n| WG_PE | required | **required** | required | - |\n\n`normal_alignments` is used only by the WG step, for tumour/normal somatic calling, and it is REQUIRED there. Do not supply it in PE mode: the PE step does not pass the normal to WISP, so it would be staged (CRAM conversion, fixmate, merge) at real cost and then discarded. \n\nNormal alignments is required because without a matched normal, SAGE has no reference against which to subtract germline variants, so the primary somatic call set is dominated by germline sites. Those are present in the patient own cfDNA at heterozygous and homozygous frequencies, and WISP measures them at high VAF and reports the result as tumour fraction.\n\n### Inputs with mixed platform (Illumina or Ultima)\n\noncoanalyser applies a single --sequencing_platform to a whole pipeline run and never checks it against the BAM headers, so two samples of different platforms in one run means one of them is analysed with the wrong error model, silently.\n\nNote that a run here means one oncoanalyser (Nextflow) invocation, not one WDL job. WG_PE launches two runs, so it can legitimately span platforms: the WG run uses the primary's platform and the PE run uses the longitudinal sample's.\n\nPlatform is read per sample from the @RG PL tag; The wdl input `sequencing_platform` overrides it for data with a missing or wrong tag. Note:\n\n* fixmate is applied only to Illumina samples. Ultima reads are single-end, and fixmate would drop every record and leave a header-only BAM.\n* the WG run requires its tumour and normal to agree, and refuses to launch otherwise. `allow_mixed_platforms` overrides this, at the cost of one sample being analysed with the wrong error model.\n\n### Note on deliverables of wdl\n\n**Germline calls are generated but not delivered.** oncoanalyser calls germline variants whenever a matched normal is present, and this cannot be switched off from configuration. They are therefore still produced, but excluded from the WG results because MRD assay does not use them. Set `include_germline_outputs` to true to keep `sage/germline/`, `pave/germline/` and the PURPLE germline files. \n\n**The WG archive carries two somatic VCFs.** `<sample>.purple.somatic.vcf.gz` is PURPLE's full call set, untouched. `<sample>.purple.somatic.prefiltered.vcf.gz` is the same call set reduced to the sites that can carry MRD signal, by the primary filters (mappability, repeat count, SNV only, tier, nearby indel, subclonal) and by the germline filter that WISP documents but cannot apply, since oncoanalyser gives it no reference genotype. `primary_site_report` gives the per-filter breakdown. In PE mode `use_primary_filters` chooses which of the two the plasma stage works from; note that WISP applies the primary filters itself either way, and records the reason per site, so the prefiltered VCF is a deliverable rather than a correction.\n\n**LOH is not available in any configuration this workflow can currently produce.** Purity therefore comes from SNVs and COBALT copy number only. "
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
            primary_site_report: {
                description: "Plain-text report of the primary tumour's variant list: how many candidate sites survive each filter in turn, primary filters and the germline filter alike. The final count is the number of sites available for MRD assessment, which is what says whether a plasma sample is worth taking. Produced in WG and WG_PE mode.",
                vidarr_label: "primarySiteReport"
            },
            pipeline_info: {
                description: "Tarball of the Nextflow pipeline_info/ directory (execution report, timeline, trace, DAG, params JSON, software versions); always produced. In WG_PE mode the PE run's copy is used.",
                vidarr_label: "pipelineInfo"
            }
        }
    }

    # Whichever stage ran last produced the pipeline_info; in WG_PE that is the PE stage.
    # Every mode runs at least one Nextflow stage, so one of the two is always defined. It is
    # computed here rather than in the output block because an output whose value is a
    # compound expression cannot be wrapped by the Vidarr output preprocessor; outputs must
    # be plain references.
    File pipeline_info_selected = select_first([subject_purity.pipeline_info_tarball,
                                               run_wgts.pipeline_info_tarball])

    output {
        File? wg_tarball          = pack_wgts.wgts_tarball   # WG and WG_PE only
        File? wisp_tarballs       = collect_results.wisp_tarballs
        File? wisp_summary        = collect_results.wisp_summary
        File? primary_site_report = pre_filtering.site_report  # WG and WG_PE only
        File  pipeline_info       = pipeline_info_selected
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
        String  scheduler
        String  slurm_partition = ""
        String  slurm_account   = ""
        Array[String] singularity_binds = []
        Array[String] nextflow_config = []
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
        scheduler:           "Requested batch scheduler"
        slurm_partition:     "Partition for the jobs nextflow submits; required for slurm"
        slurm_account:       "Accounting group for those jobs; checked only for consistency here"
        singularity_binds:   "Container bind paths; checked only for consistency here"
        nextflow_config:     "Extra nextflow config files, checked for readability"
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

      # A tumour-only primary does not merely degrade the MRD result, it INVERTS it. With no
      # matched normal SAGE cannot subtract germline variants, so the somatic call set is
      # dominated by germline sites; those sit in the patient's own cfDNA at heterozygous and
      # homozygous frequencies, and WISP reports that as tumour fraction. Unconditional hard
      # error, no override: a tumour-only primary has no legitimate use here.
      if [ "${mode}" != "PE" ] && ~{if has_normal then "false" else "true"}; then
        errors+=("mode ${mode} requires normal_alignments: without a matched normal the primary is called tumour-only, germline variants are not subtracted, and WISP reports a FALSE MRD-POSITIVE. There is deliberately no override")
      fi

      # These paths are String, not File, so Cromwell neither localizes nor checks them; a
      # wrong one would otherwise surface only when nextflow parses -c, by which point the
      # alignment work is done. An entry holding an env var is resolved in the head job, where
      # the module is loaded, so here it can only be counted.
      config_count=0
      while IFS= read -r cfg; do
        [ -n "${cfg}" ] || continue
        config_count=$((config_count + 1))
        case "${cfg}" in *'$'*) continue ;; esac
        [ -f "${cfg}" ] && [ -r "${cfg}" ] || errors+=("nextflow_config not readable: ${cfg}")
      done < "~{write_lines(nextflow_config)}"

      # Which scheduler the head jobs tell nextflow to submit to. Decided here, once, from
      # the submit command the cluster provides, and passed on to those tasks; a caller that
      # had to name it per site would get it wrong the first time the inputs moved. sbatch is
      # preferred where both exist, and the choice is always logged so it can be checked.
      sched="~{scheduler}"
      if [ -z "${sched}" ]; then
        if   command -v sbatch >/dev/null 2>&1; then sched=slurm
        elif command -v qsub   >/dev/null 2>&1; then sched=sge
        else
          errors+=("cannot tell which scheduler nextflow should submit to: neither sbatch nor qsub is on PATH. Set the scheduler input")
        fi
        [ -z "${sched}" ] || echo "detected scheduler ${sched}" >&2
      fi
      echo "${sched}" > scheduler.txt

      # Settings that only a slurm run reads. Supplied elsewhere they are ignored rather than
      # refused, so one set of inputs can carry them and still run at a site that does not
      # need them.
      slurm_only=()
      if [ -n "~{slurm_partition}" ]; then slurm_only+=("slurm_partition"); fi
      if [ -n "~{slurm_account}" ];   then slurm_only+=("slurm_account"); fi
      binds=(~{sep=" " singularity_binds})
      if [ "${#binds[@]}" -gt 0 ];    then slurm_only+=("singularity_binds"); fi

      case "${sched}" in
        sge)
          if [ "${#slurm_only[@]}" -gt 0 ]; then
            echo "note: $(IFS=,; echo "${slurm_only[*]}") ignored; those apply only to slurm" >&2
          fi
          ;;
        slurm)
          # The overlay the module ships names a queue for its own scheduler. The head jobs
          # generate the rest of the slurm settings, but the partition has to be supplied.
          [ -n "~{slurm_partition}" ] || \
            errors+=("scheduler slurm requires slurm_partition: the overlay the oncoanalyser module ships names a queue for its own scheduler, which does not exist here")
          ;;
        "") ;;
        *) errors+=("scheduler must be sge or slurm, got '${sched}'") ;;
      esac

      if [ "${#errors[@]}" -gt 0 ]; then
        echo "ERROR: inputs do not satisfy mode ${mode}:" >&2
        for e in "${errors[@]}"; do
          echo "  - ${e}" >&2
        done
        exit 1
      fi

      if [ "${mode}" = "PE" ] && ~{if has_normal then "true" else "false"}; then
        {
          echo "note: normal_alignments is not used in PE mode. The PE run processes only the"
          echo "      longitudinal sample, so the normal would be staged (CRAM conversion,"
          echo "      fixmate, merge) at real cost and then discarded. Omit it."
        } >&2
      fi

      echo "inputs OK for mode ${mode}"
    >>>

    output {
        String checked   = read_string(stdout())
        String scheduler_used = read_string("scheduler.txt")
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
        String modules = "samtools/1.16.1 oncoanalyser-data/3.0.0"
        Int threads = 8
        Int memory  = 16
        Int timeout = 24
    }

    parameter_meta {
        aln:            "Input CRAM file"
        idx:            "CRAM index (.crai) for aln"
        cram_reference: "Path to the FASTA the CRAM was encoded against, hs38DH for Ultima vendor CRAMs. Normally the literal $VENDOR_GENOME_HS38DH, expanded by the shell after the data module loads. A String so it is read from shared storage without localization"
        modules:        "Environment modules to load. samtools is required; oncoanalyser-data supplies $VENDOR_GENOME_HS38DH for cram_reference"
        threads:        "Number of samtools threads"
        memory:         "Memory in GB"
        timeout:        "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # Use no more threads than the scheduler actually allocated. A backend that does not
      # map the cpu runtime attribute leaves samtools oversubscribed on one core, which with
      # its block queues is far slower than simply running single-threaded -- and silently so.
      # The note makes a misconfigured backend visible in the log instead of as a mystery.
      # Prefer the scheduler's own statement of the allocation over nproc. Slurm only
      # constrains CPU affinity when configured to, so with a CFS quota instead nproc reports
      # the whole node and the cap would silently do nothing.
      # SLURM_CPUS_PER_TASK exists only when --cpus-per-task was passed -- which is precisely
      # what a backend missing that flag does not do -- so SLURM_CPUS_ON_NODE is the one to
      # rely on: it is always set, and equals the allocation for a one-task job. NSLOTS is
      # SGE's equivalent.
      avail=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${NSLOTS:-}}}
      [ -n "${avail}" ] || avail=$(nproc 2>/dev/null || echo 1)
      threads=~{threads}
      if [ "${avail}" -lt "${threads}" ]; then
        echo "NOTE: ~{threads} threads requested but only ${avail} cpu(s) allocated; using ${avail}." >&2
        echo "      If this is not deliberate, the backend is not mapping the cpu attribute." >&2
        threads=${avail}
      fi

      # Cromwell localizes the CRAM and its index into different directories, so symlink
      # both here with matching basenames for samtools to find the index.
      ln -s ~{aln} ./input.cram
      ln -s ~{idx} ./input.cram.crai

      [ -s "~{cram_reference}" ] || {
        echo "ERROR: CRAM reference not readable: ~{cram_reference}" >&2
        exit 1
      }

      samtools view -@ "${threads}" -T ~{cram_reference} -b -o converted.bam ./input.cram
      samtools index -@ "${threads}" converted.bam
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

# Does this alignment already carry mate CIGAR tags?
#
# REDUX needs MC tags to mark duplicates correctly, and fixmate is how they get there. Some
# aligners write them already, in which case the per-chromosome fixmate scatter is pure cost.
# Deciding that from the data rather than from an input means doFixmate=false cannot silently
# hand REDUX an alignment it will mis-process.
#
# Reads records rather than the header, so a bounded slice is used: MC is per-read, and the
# question is whether the writer emits it at all, which the first records answer.
task probe_mc_tags {
    input {
        File   bam
        Int    records = 100000
        String modules = "samtools/1.16.1"
        Int memory  = 4
        Int timeout = 1
    }

    parameter_meta {
        bam:     "Alignment to inspect. Always a BAM, because the probe is only reached for platforms that get fixmate, and those inputs are BAMs; so no reference is needed to read records"
        records: "How many leading records to inspect"
        modules: "Environment modules to load (samtools required)"
        memory:  "Memory in GB"
        timeout: "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # head closes the pipe, and grep -c exits 1 when it matches nothing, so neither can be
      # allowed to fail the task.
      set +o pipefail
      n=$(samtools view "~{bam}" | head -~{records} | grep -c 'MC:Z:' || true)
      set -o pipefail
      if [ "${n}" -gt 0 ]; then
        echo "true" > has_mc_tags.txt
        echo "mate CIGAR tags present (${n} in the first ~{records} records)" >&2
      else
        echo "false" > has_mc_tags.txt
        echo "no mate CIGAR tags in the first ~{records} records" >&2
      fi
    >>>

    output {
        Boolean has_mc_tags = read_boolean("has_mc_tags.txt")
    }

    runtime {
        cpu:     1
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
      # Use no more threads than the scheduler actually allocated. A backend that does not
      # map the cpu runtime attribute leaves samtools oversubscribed on one core, which with
      # its block queues is far slower than simply running single-threaded -- and silently so.
      # The note makes a misconfigured backend visible in the log instead of as a mystery.
      # Prefer the scheduler's own statement of the allocation over nproc. Slurm only
      # constrains CPU affinity when configured to, so with a CFS quota instead nproc reports
      # the whole node and the cap would silently do nothing.
      # SLURM_CPUS_PER_TASK exists only when --cpus-per-task was passed -- which is precisely
      # what a backend missing that flag does not do -- so SLURM_CPUS_ON_NODE is the one to
      # rely on: it is always set, and equals the allocation for a one-task job. NSLOTS is
      # SGE's equivalent.
      avail=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${NSLOTS:-}}}
      [ -n "${avail}" ] || avail=$(nproc 2>/dev/null || echo 1)
      threads=~{threads}
      if [ "${avail}" -lt "${threads}" ]; then
        echo "NOTE: ~{threads} threads requested but only ${avail} cpu(s) allocated; using ${avail}." >&2
        echo "      If this is not deliberate, the backend is not mapping the cpu attribute." >&2
        threads=${avail}
      fi

      # Symlink BAM and BAI preserving the original filenames so samtools can
      # locate the index automatically (it derives the index path from the BAM
      # path; renaming would break that lookup).
      ln -s ~{bam} .
      ln -s ~{bai} .
      bam_name=$(basename ~{bam})
      # -e 'rnext == rname' rather than awk on $7=="=": RNEXT is written as "=" when the
      # mate is on the same reference, but a writer may spell the reference name out instead,
      # and the samtools expression means what is intended either way.
      samtools view -h -@ 2 -e 'rnext == rname' "${bam_name}" ~{chr} \
        | samtools sort -n -u -@ 2 -m 1G - \
        | samtools fixmate -m -u -@ 2 - - \
        | samtools sort -@ "${threads}" -m 2G -o ~{chr}.fixedmate.bam -
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
# from one lane BAM. Extracts only reads where RNEXT != "=", a small fraction of the total,
# so the intermediate file stays small.
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
      # Use no more threads than the scheduler actually allocated. A backend that does not
      # map the cpu runtime attribute leaves samtools oversubscribed on one core, which with
      # its block queues is far slower than simply running single-threaded -- and silently so.
      # The note makes a misconfigured backend visible in the log instead of as a mystery.
      # Prefer the scheduler's own statement of the allocation over nproc. Slurm only
      # constrains CPU affinity when configured to, so with a CFS quota instead nproc reports
      # the whole node and the cap would silently do nothing.
      # SLURM_CPUS_PER_TASK exists only when --cpus-per-task was passed -- which is precisely
      # what a backend missing that flag does not do -- so SLURM_CPUS_ON_NODE is the one to
      # rely on: it is always set, and equals the allocation for a one-task job. NSLOTS is
      # SGE's equivalent.
      avail=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${NSLOTS:-}}}
      [ -n "${avail}" ] || avail=$(nproc 2>/dev/null || echo 1)
      threads=~{threads}
      if [ "${avail}" -lt "${threads}" ]; then
        echo "NOTE: ~{threads} threads requested but only ${avail} cpu(s) allocated; using ${avail}." >&2
        echo "      If this is not deliberate, the backend is not mapping the cpu attribute." >&2
        threads=${avail}
      fi

      # -f 1:  paired
      # -F 12: neither read nor mate unmapped
      # -e:    mate on a different reference; the complement of the fixmate_chr expression
      # No index needed: full-file sequential scan.
      samtools view -h -@ 2 -f 1 -F 12 -e 'rnext != rname' ~{bam} \
        | samtools sort -n -u -@ 2 -m 1G - \
        | samtools fixmate -m -u -@ 2 - - \
        | samtools sort -@ "${threads}" -m 2G -o discordant.fixedmate.bam -
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

# Merge coordinate-sorted BAMs, and check the header is fit for what consumes it.
#
#   repair_flags     per-chromosome fixmate can clear 0x1 while leaving 0x40/0x80 set on
#                    discordant-mate reads, which htsjdk treats as a validation error. Only
#                    the fixmate path needs this.
#   strip_bad_pg     set when the output feeds REDUX, i.e. htsjdk. An @PG line is only
#                    stripped if it is ACTUALLY malformed, so a well-formed header keeps its
#                    provenance; see the detection in the command.
task merge_bams {
    input {
        Array[File] bams
        Array[File] bais
        Boolean     repair_flags = false
        Boolean     strip_bad_pg = false
        String      modules = "samtools/1.16.1"
        Int threads = 8
        Int memory  = 16
        Int timeout = 24
    }

    parameter_meta {
        bams:           "Coordinate-sorted BAMs to merge"
        bais:           "Index files corresponding to each entry in bams (same order)"
        repair_flags:   "Repair reads whose paired flag was cleared while the first/second-in-pair flags were left set. Needed only after per-chromosome fixmate"
        strip_bad_pg:   "Check the @PG header and strip it if malformed. Set when the output feeds REDUX; a well-formed header is left alone"
        modules:        "Environment modules to load (samtools required)"
        threads:        "Number of samtools threads"
        memory:         "Memory in GB"
        timeout:        "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # Use no more threads than the scheduler actually allocated. A backend that does not
      # map the cpu runtime attribute leaves samtools oversubscribed on one core, which with
      # its block queues is far slower than simply running single-threaded -- and silently so.
      # The note makes a misconfigured backend visible in the log instead of as a mystery.
      # Prefer the scheduler's own statement of the allocation over nproc. Slurm only
      # constrains CPU affinity when configured to, so with a CFS quota instead nproc reports
      # the whole node and the cap would silently do nothing.
      # SLURM_CPUS_PER_TASK exists only when --cpus-per-task was passed -- which is precisely
      # what a backend missing that flag does not do -- so SLURM_CPUS_ON_NODE is the one to
      # rely on: it is always set, and equals the allocation for a one-task job. NSLOTS is
      # SGE's equivalent.
      avail=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${NSLOTS:-}}}
      [ -n "${avail}" ] || avail=$(nproc 2>/dev/null || echo 1)
      threads=~{threads}
      if [ "${avail}" -lt "${threads}" ]; then
        echo "NOTE: ~{threads} threads requested but only ${avail} cpu(s) allocated; using ${avail}." >&2
        echo "      If this is not deliberate, the backend is not mapping the cpu attribute." >&2
        threads=${avail}
      fi

      bam_list=(~{sep=" " bams})

      # Strip @PG only when it is actually malformed. An aligner can embed tabs inside CL:,
      # and because the header line is itself tab-delimited those become extra fields --
      # including a second ID: -- which is what htsjdk rejects. Counting ID: fields detects
      # exactly that, so a well-formed header keeps its provenance.
      strip_pg=false
      if ~{strip_bad_pg}; then
        bad_pg=$(samtools view -H "${bam_list[0]}" \
                 | awk -F'\t' '/^@PG/ {n=0; for (i=1;i<=NF;i++) if ($i ~ /^ID:/) n++; if (n>1) c++}
                               END {print c+0}')
        if [ "${bad_pg}" -gt 0 ]; then
          echo "@PG header has an embedded tab (${bad_pg} line(s) with more than one ID:);" >&2
          echo "stripping @PG so htsjdk can read the header" >&2
          strip_pg=true
        else
          echo "@PG header is well formed; keeping it" >&2
        fi
      fi

      merge_stream() {
        if [ "${#bam_list[@]}" -gt 1 ]; then
          samtools merge -f -c -p -u -@ "${threads}" - "${bam_list[@]}"
        else
          samtools view -h -u "${bam_list[0]}"
        fi
      }

      if ${strip_pg} || ~{repair_flags}; then
        merge_stream \
          | samtools view -h -@ "${threads}" \
          | awk -v strip_pg="${strip_pg}" -v repair="~{repair_flags}" 'BEGIN{OFS="\t"}
              /^@PG/ { if (strip_pg == "true") next; print; next }
              /^@/   { print; next }
              {
                if (repair == "true") {
                  f = int($2)
                  if (and(f,1) == 0) { f = f - and(f,64) - and(f,128); $2 = f }
                }
                print
              }' \
          | samtools view -@ "${threads}" -O BAM -o merged.bam
      else
        merge_stream | samtools view -@ "${threads}" -O BAM -o merged.bam
      fi

      samtools index -@ "${threads}" merged.bam
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

      # Derive the primary tumour sample ID from the PURPLE filenames rather than trusting a
      # hand-typed input. SAGE_APPEND reads <tumor_sample_id>.purple.somatic.vcf.gz, so a wrong
      # ID would fail only after the expensive steps have run. The WG step names these files
      # from the tumour @RG SM tag, which is rarely the same string as group_id.
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

# Filter the primary's somatic VCF and report what each filter cost.
#
# TWO FILTER SETS, both applied here because both are primary-side and both need the same
# VCF and the same sample-column resolution:
#
#   PRIMARY FILTERS   the criteria WISP uses to decide which sites can carry MRD signal:
#                     mappability, repeat count, SNV only, tier, nearby indel, subclonal.
#
#                     Two of WISP's ten conditions are missing ON PURPOSE, because both are
#                     measured on the ctDNA sample: the high-confidence check
#                     (RC_QUAL/RC_CNT >= 18) and average edge distance (AED[1] >= 0.06).
#                     WISP evaluates them on the cfDNA-annotated VCF, which only exists after
#                     SAGE_APPEND, and applies both itself.
#
#                     AED is a trap worth naming. It is a FORMAT field, so it HAS a value on
#                     the primary (Number=2, "[alt,total]"), and FORMAT/AED[<tumour>:1] >= 0.06
#                     runs happily while removing almost nothing -- because the tumour's edge
#                     distance is not the quantity the condition means. In WISP's own per-site
#                     output the equivalent column is zero at every site with no plasma alt
#                     read, which is what settles it. A field being present is not the same as
#                     the condition being applicable.
#   GERMLINE FILTER   sites with support in the matched normal. WISP documents this one but
#                     never applies it, because oncoanalyser does not pass the reference
#                     sample id through and WISP therefore has no genotype to test.
#
# The result is written beside the original as <sample>.purple.somatic.prefiltered.vcf.gz.
# The original is left untouched and BOTH go into the WG archive, so a reader of that archive
# sees the full call set and the MRD-usable subset.
#
# Only the PASS step drops records. Every condition after it SOFT-FLAGS, adding its name to
# the FILTER column, so the prefiltered VCF is self-describing: it says per site which
# conditions failed instead of silently omitting the site, and FILTER="PASS" selects the
# usable subset. That is also what WISP does with whatever VCF it is handed, so passing this
# file to the plasma stage gives the same site set as hard filtering would. The per-filter table is what says how
# much signal a primary can support, which is the question a clinical report has to answer
# before a plasma sample is worth taking.
#
# Note that WISP applies the primary filters itself, on the sites it is given, and records
# the reason per site. Pre-filtering here does not replace that; it produces the reduced call
# set as a deliverable. Which of the two VCFs a later PE run works from is chosen there, by
# use_primary_filters, not here.
task pre_filtering {
    input {
        String  purple_dir
        String  tumor_sample_id
        String  outputFileNamePrefix
        String  modules = "bcftools/1.9"
        Int memory  = 4
        Int timeout = 2
    }

    parameter_meta {
        purple_dir:       "PURPLE output directory of the primary, holding <tumor_sample_id>.purple.somatic.vcf.gz"
        tumor_sample_id:  "Primary tumour sample ID; also used to identify which VCF column is the tumour, and hence which is the normal"
        outputFileNamePrefix: "Prefix for the report filename, so runs of different samples do not provision the same name"
        modules:          "Environment modules to load (bcftools required)"
        memory:           "Memory in GB"
        timeout:          "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail

      SRC="~{purple_dir}"
      TUMOR="~{tumor_sample_id}"
      VCF_NAME="${TUMOR}.purple.somatic.vcf.gz"
      PREFILTERED_NAME="${TUMOR}.purple.somatic.prefiltered.vcf.gz"
      SRC_VCF="${SRC}/${VCF_NAME}"
      DEST="$(pwd)/purple_prefiltered"
      REPORT="~{outputFileNamePrefix}.primary_site_report.txt"
      : > "${REPORT}"

      # Anything that stops us filtering is reported and then ignored: an unusable VCF must
      # not take down a run that would otherwise have produced a result. The caller falls
      # back to the unfiltered directory.
      give_up() {
        echo "NOTE: pre-filtering skipped: $1" >&2
        { echo "prefiltering: SKIPPED"; echo "reason: $1"; } >> "${REPORT}"
        echo "-1" > usable_sites.txt
        echo "0"  > germline_removed.txt
        echo "${SRC}" > purple_dir_out.txt
        exit 0
      }

      [ -s "${SRC_VCF}" ] || give_up "no ${VCF_NAME} in ${SRC}"

      # Sample columns resolved BY NAME, never by position. If the columns are the other way
      # round a hardcoded index tests the tumour instead of the normal, which removes the
      # entire call set and yields a confident MRD-negative with no error anywhere.
      mapfile -t SAMPLES < <(bcftools query -l "${SRC_VCF}")
      [ "${#SAMPLES[@]}" -eq 2 ] || give_up "expected 2 samples in ${VCF_NAME}, found ${#SAMPLES[@]}"
      TUM_IDX=""; NORM_IDX=""
      for i in 0 1; do
        if [ "${SAMPLES[$i]}" = "${TUMOR}" ]; then TUM_IDX=$i; NORM_IDX=$((1 - i)); fi
      done
      [ -n "${TUM_IDX}" ] || give_up "tumour '${TUMOR}' not among VCF samples: ${SAMPLES[*]}"

      HDR=$(bcftools view -h "${SRC_VCF}")
      has() { grep -q "^##$1=<ID=$2," <<<"${HDR}"; }

      # The quality field has been renamed across hmftools versions. ARCBQ is the name the
      # WISP filter is documented against and what current SAGE emits; ABQ is the older
      # equivalent and is still what externally supplied VCFs may carry. RABQ is
      # deliberately NOT accepted: it is the raw, pre-recalibration value, so applying the
      # threshold to it would filter on the wrong quantity. Better to skip the filter and
      # say so than to apply it to a field that does not mean what the threshold assumes.
      QUAL=""
      if has FORMAT ARCBQ; then QUAL=ARCBQ; elif has FORMAT ABQ; then QUAL=ABQ; fi
      REPC_FIELDS=()
      has INFO RC_REPC && REPC_FIELDS+=(RC_REPC)
      has INFO REP_C   && REPC_FIELDS+=(REP_C)

      GERMLINE_EXPR=""
      if [ -n "${QUAL}" ]; then
        GERMLINE_EXPR="(FORMAT/AD[${NORM_IDX}:1]/FORMAT/DP[${NORM_IDX}]) > 0.01 & FORMAT/${QUAL}[${NORM_IDX}:1] > 30"
      fi

      # tag|mode|expression. mode i = flag when the expression is FALSE, e = flag when TRUE.
      # The germline filter has to be an exclude: bcftools cannot negate an indexed FORMAT
      # expression. Tags are what land in the VCF FILTER column, so they are uppercase and
      # match the condition names the WISP documentation uses.
      FILTERS=()
      has INFO MAPPABILITY  && FILTERS+=("MAPPABILITY|i|INFO/MAPPABILITY>=0.5")
      for rf in "${REPC_FIELDS[@]+"${REPC_FIELDS[@]}"}"; do
        FILTERS+=("${rf}|i|INFO/${rf}<4 || INFO/${rf}==\".\"")
      done
      FILTERS+=("NON_SNV|i|TYPE=\"snp\"")
      has INFO TIER         && FILTERS+=("LOW_CONFIDENCE|i|INFO/TIER!=\"LOW_CONFIDENCE\"")
      has INFO NEARBY_INDEL && FILTERS+=("NEARBY_INDEL|i|INFO/NEARBY_INDEL=0")
      if has INFO SUBCL && has INFO PURPLE_VCN; then
        FILTERS+=("SUBCLONAL|i|INFO/SUBCL<=0.5 || INFO/PURPLE_VCN>=0.7")
      fi
      [ -n "${GERMLINE_EXPR}" ] && FILTERS+=("GERMLINE|e|${GERMLINE_EXPR}")

      {
        echo "primary sample: ${TUMOR}"
        echo "vcf:            ${VCF_NAME}"
        echo "samples:        [0]=${SAMPLES[0]} [1]=${SAMPLES[1]}"
        echo "normal column:  ${NORM_IDX} (${SAMPLES[$NORM_IDX]})"
        echo "qual field:     ${QUAL:-NONE, germline filter not applied}"
        echo "repeat fields:  ${REPC_FIELDS[*]+"${REPC_FIELDS[*]}"}"
        echo
        echo "Two of WISP's conditions are measurements on the cfDNA sample and cannot be"
        echo "evaluated here: average edge distance (AED[1] >= 0.06) and the quality ratio"
        echo "(RC_QUAL[0]+[1]+[3])/(RC_CNT[0]+[1]+[3]) >= 18. WISP applies both itself, so the"
        echo "count below is an upper bound on the sites it will actually use."
        echo
        printf '%-16s %10s %10s\n' "filter" "remaining" "lost"
      } >> "${REPORT}"

      WORK=$(mktemp -d)
      TOTAL=$(bcftools view -H "${SRC_VCF}" | wc -l)
      printf '%-16s %10s %10s\n' "(total)" "${TOTAL}" "-" >> "${REPORT}"

      # Only this first step DROPS records. Everything after it SOFT-FLAGS: the record set is
      # fixed here and later steps only add tags to the FILTER column, so the delivered VCF
      # says per site which conditions it failed instead of silently omitting it. A consumer
      # selects FILTER="PASS" to get the usable subset, which is also what WISP does with
      # whatever VCF it is given.
      bcftools view -i 'FILTER="PASS"' -Oz -o "${WORK}/cur.vcf.gz" "${SRC_VCF}"
      bcftools index -t -f "${WORK}/cur.vcf.gz"
      PREV=$(bcftools view -H "${WORK}/cur.vcf.gz" | wc -l)
      KEPT=${PREV}
      printf '%-16s %10s %10s\n' "PASS" "${PREV}" "$(( TOTAL - PREV ))" >> "${REPORT}"

      GERMLINE_LOST=""
      for entry in "${FILTERS[@]}"; do
        tag="${entry%%|*}"; rest="${entry#*|}"; mode="${rest%%|*}"; expr="${rest#*|}"
        # A field can exist in the header and still be unusable by this bcftools build, so
        # a failing expression is reported and skipped rather than killing the run.
        if ! bcftools filter -s "${tag}" -m + "-${mode}" "${expr}" \
             -Oz -o "${WORK}/next.vcf.gz" "${WORK}/cur.vcf.gz" 2>"${WORK}/err"; then
          printf '%-16s %10s %10s\n' "${tag}" "SKIPPED" "$(head -1 "${WORK}/err" | cut -c1-40)" >> "${REPORT}"
          continue
        fi
        mv "${WORK}/next.vcf.gz" "${WORK}/cur.vcf.gz"
        bcftools index -t -f "${WORK}/cur.vcf.gz"
        # "remaining" is how many records are still PASS, not how many are left in the file.
        n=$(bcftools view -H -i 'FILTER="PASS"' "${WORK}/cur.vcf.gz" | wc -l)
        printf '%-16s %10s %10s\n' "${tag}" "${n}" "$(( PREV - n ))" >> "${REPORT}"

        if [ "${tag}" = "GERMLINE" ]; then GERMLINE_LOST=$(( PREV - n )); fi
        PREV=${n}
      done
      USABLE=${PREV}
      echo "${USABLE}" > usable_sites.txt
      echo "${GERMLINE_LOST:-0}" > germline_removed.txt

      # Copy rather than write in place: the source is another task's output or a tarball
      # extraction, and neither may be mutated.
      mkdir -p "${DEST}"
      cp -r "${SRC}/." "${DEST}/"
      cp "${WORK}/cur.vcf.gz" "${DEST}/${PREFILTERED_NAME}"
      bcftools index -t -f "${DEST}/${PREFILTERED_NAME}"
      echo "${DEST}" > purple_dir_out.txt

      { echo
        echo "sites usable for MRD: ${USABLE} of ${TOTAL}"
        echo "written as ${PREFILTERED_NAME}, alongside the unmodified ${VCF_NAME}."
        echo "That file holds all ${KEPT} PASS-in-purple records; the ones that failed a"
        echo "condition carry it in the FILTER column rather than being dropped, so select"
        echo "FILTER=\"PASS\" for the ${USABLE} usable sites."
      } >> "${REPORT}"

      rm -rf "${WORK}"
      cat "${REPORT}" >&2
    >>>

    output {
        Int    usable_sites    = read_int("usable_sites.txt")
        Int    germline_removed = read_int("germline_removed.txt")
        String purple_dir_out  = read_string("purple_dir_out.txt")
        File   site_report     = "~{outputFileNamePrefix}.primary_site_report.txt"
    }

    runtime {
        cpu:     1
        timeout: "~{timeout}"
        memory:  "~{memory} GB"
        modules: "~{modules}"
    }
}

# Archive the WG results. Separate from run_wgts so it can run AFTER pre_filtering and
# therefore carry both the original and the prefiltered somatic VCF.
task pack_wgts {
    input {
        String  wgts_dir
        String  purple_dir
        String  outputFileNamePrefix
        Boolean include_germline_outputs
        Int memory  = 8
        Int timeout = 4
    }

    parameter_meta {
        wgts_dir:                 "The WG output directory, supplying amber/, cobalt/, pave/ and sage/"
        purple_dir:               "purple/ as written by pre_filtering: the pipeline's own output plus the prefiltered somatic VCF"
        outputFileNamePrefix:     "Prefix for the archive filename"
        include_germline_outputs: "Whether to keep germline variant calls in the archive; see the workflow input of the same name"
        memory:                   "Memory in GB"
        timeout:                  "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      # Stage symlinks and dereference them into the archive, so the WG output is not copied
      # a second time. alignments/ is excluded on purpose: it holds full REDUX BAMs, and a
      # chained PE run reads the primary outputs from disk rather than from this archive.
      stage="$(pwd)/stage"
      mkdir -p "${stage}"

      # Collect only the directories that exist, so a run that legitimately produces fewer
      # of them does not fail at the very last step after hours of compute.
      tar_dirs=()
      for d in amber cobalt pave sage; do
        if [ -d "~{wgts_dir}/${d}" ]; then
          ln -s "~{wgts_dir}/${d}" "${stage}/${d}"
          tar_dirs+=("${d}/")
        else
          echo "note: ${d}/ not present in the output, omitting from the archive" >&2
        fi
      done
      if [ -d "~{purple_dir}" ]; then
        ln -s "~{purple_dir}" "${stage}/purple"
        tar_dirs+=("purple/")
      else
        echo "note: purple/ not present, omitting from the archive" >&2
      fi
      if [ "${#tar_dirs[@]}" -eq 0 ]; then
        echo "ERROR: none of amber/ cobalt/ purple/ pave/ sage/ were produced" >&2
        exit 1
      fi

      # Germline calls are not part of the MRD deliverable and WISP does not use them, so
      # they are dropped from the archive. They ARE still generated: oncoanalyser hardcodes
      # germline calling on and it cannot be switched off from configuration.
      # The pattern catches sage/germline/, pave/germline/ and the PURPLE germline files.
      exclude_args=()
      if ! ~{include_germline_outputs}; then
        exclude_args=(--exclude='*germline*')
      fi

      # -h dereferences the staged symlinks; without it the archive holds four links.
      tar -czhf ~{outputFileNamePrefix}.wgts.tar.gz \
          "${exclude_args[@]}" \
          -C "${stage}" \
          "${tar_dirs[@]}"
    >>>

    output {
        File wgts_tarball = "~{outputFileNamePrefix}.wgts.tar.gz"
    }

    runtime {
        cpu:     1
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
        String  scheduler
        String  slurm_partition = ""
        String  slurm_account   = ""
        Array[String] singularity_binds = []
        Array[String] nextflow_config
        String  submit_wrapper
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
        ref_data_dir:        "HMF reference data directory, used for --igenomes_base, --hmf_genomes_base and --ref_data_hmf_data_path; normally the literal $REFERENCE_FILES_DIR"
        images_dir:          "Singularity image cache directory (NXF_SINGULARITY_CACHEDIR); normally the literal $IMAGES_DIR, expanded by the shell after the oncoanalyser module loads"
        pipeline_dir:        "oncoanalyser checkout containing main.nf; normally the literal $ONCOANALYSER_FOLDER"
        nextflow_bin:        "Nextflow executable; defaults to `nextflow` on PATH, with the version pinned by the module's NXF_VER. Must resolve to 25.10.0 or newer"
        nextflow_home:       "NXF_HOME holding the pre-cached nf-schema plugin, so the run works with NXF_OFFLINE=true"
        scheduler:       "Batch scheduler Nextflow submits to; the submit wrapper is installed only for sge"
        slurm_partition:   "Partition for the jobs nextflow submits, when scheduler is slurm"
        slurm_account:     "Accounting group for those jobs; also clears the request the module's overlay writes in the other scheduler's syntax"
        singularity_binds: "Paths bound into every container, replacing the bind the module's overlay sets. Empty keeps that bind"
        nextflow_config:   "Config overlays, each passed with its own -c, in the order given"
        submit_wrapper:      "Wrapper placed on PATH ahead of the scheduler submit command, to adjust the resource requests Nextflow generates"
        modules:             "Environment modules to load"
        memory:              "Memory in GB for this task, which hosts the Nextflow driver only; the pipeline processes get their own allocations"
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

      # This run processes the tumour AND the normal, but oncoanalyser applies ONE
      # --sequencing_platform per run, so a tumour/normal pair split across platforms would
      # get the wrong error model applied to one of them. There is no way to drop the normal
      # from a tumour/normal analysis, so refuse unless explicitly overridden.
      if [ -n "${normal_platform}" ] && [ "${normal_platform}" != "~{sequencing_platform}" ]; then
        {
          echo "mixed sequencing platforms within a SINGLE oncoanalyser run:"
          echo "  this run uses --sequencing_platform ~{sequencing_platform}"
          echo "  but the normal sample is ${normal_platform}"
        } >&2
        if ~{allow_mixed_platforms}; then
          echo "WARNING: proceeding because allow_mixed_platforms=true; results for the" >&2
          echo "         ${normal_platform} sample are not trustworthy." >&2
        else
          {
            echo "REFUSING TO LAUNCH. Use a tumour/normal pair from one platform, or set"
            echo "allow_mixed_platforms=true for a deliberate experiment."
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
      # -XX:-UseContainerSupport: the JVM's cgroup probe aborts startup on a node where it
      # cannot enumerate the controllers. The heap is pinned here, so nothing depends on the
      # limits that probe would report.
      export NXF_OPTS="-Xms512m -Xmx8g -XX:-UseContainerSupport"
      export NXF_SINGULARITY_CACHEDIR=~{images_dir}
      # NXF_HOME must hold the pre-cached plugins; with NXF_OFFLINE=true a wrong one fails
      # obscurely inside plugin resolution. Prefer the module's value, falling back to a
      # sibling directory that does have plugins/.
      nxf_home="~{nextflow_home}"
      if [ ! -d "${nxf_home}/plugins" ]; then
        for cand in "${ONCOANALYSER_ROOT:-}"/nextflow_home*; do
          if [ -d "${cand}/plugins" ]; then
            echo "WARNING: NEXTFLOW_HOME (${nxf_home}) has no plugins/; using ${cand}" >&2
            echo "         Fix NEXTFLOW_HOME in the oncoanalyser modulefile." >&2
            nxf_home="${cand}"
            break
          fi
        done
      fi
      [ -d "${nxf_home}/plugins" ] || {
        echo "ERROR: no NXF_HOME with a plugins/ directory. Tried ${nxf_home} and" >&2
        echo "       ${ONCOANALYSER_ROOT:-<ONCOANALYSER_ROOT unset>}/nextflow_home*" >&2
        exit 1
      }
      # NXF_HOME must be WRITABLE. The nextflow launcher does mkdir -p $NXF_HOME/tmp and the
      # local secrets provider wants $NXF_HOME/secrets, but a module tree is read-only, so
      # pointing NXF_HOME straight at it fails with "Read-only file system". Give nextflow a
      # writable directory in the task workdir. plugins/ itself must be a real directory:
      # nextflow calls createDirectories on it, which rejects a symlink even when it resolves
      # to a directory. Its entries are only read, so they can be symlinks into the module.
      nxf_home_rw="$(pwd)/nxf_home"
      mkdir -p "${nxf_home_rw}/plugins"
      for plugin in "${nxf_home}"/plugins/*; do
        if [ -e "${plugin}" ]; then
          ln -sfn "${plugin}" "${nxf_home_rw}/plugins/"
        fi
      done
      export NXF_HOME="${nxf_home_rw}"
      # A loaded oncoanalyser 2.x module points these at its read-only tree.
      unset NXF_DIST NXF_LAUNCHER NXF_PLUGINS_DIR || true

      # The SGE executor embeds h_rss/mem_free directly in .command.run, where configuration
      # cannot reach them, so a shim ahead of qsub on PATH rewrites the request. The Slurm
      # executor emits --mem and --cpus-per-task from the memory and cpus directives, so it
      # needs no wrapper.
      if [ "~{scheduler}" = "sge" ]; then
        bin_dir="$(pwd)/bin"
        mkdir -p "${bin_dir}"
        cp ~{submit_wrapper} "${bin_dir}/qsub"
        chmod +x "${bin_dir}/qsub"
        export PATH="${bin_dir}:$PATH"
      fi

      # One -c per overlay, in the order given: later files win, which is how a site config
      # overrides the shared one.
      configs=(~{sep=" " nextflow_config})
      config_args=()
      for cfg in "${configs[@]}"; do
        config_args+=(-c "${cfg}")
      done

      # The overlay the module ships selects its own scheduler, its queue and its container
      # binds. Generate the overrides here instead of requiring a file on the filesystem, so
      # the workflow carries everything it needs, and apply them last so they win. Only the
      # values that vary by site are substituted; the shape is the same everywhere.
      if [ "~{scheduler}" = "slurm" ]; then
        bind_list=(~{sep=" " singularity_binds})
        {
          echo "executor { name = 'slurm' }"
          echo "process {"
          echo "    queue = '~{slurm_partition}'"
          # Slurm enforces the memory request as RSS. A JVM sized at the module default of
          # ~95% of it leaves nothing for the helper processes the tools fork, and the cgroup
          # kills the job. Keep a quarter of the request outside the heap. Processes that do
          # not read this setting are unaffected.
          echo "    ext.xmx_mod = 0.75"
          # A cgroup kill reports no exit status, which the pipeline's own list does not
          # match, so the job would not be retried at all -- and the per-process memory
          # directives scale with task.attempt, which is otherwise unreachable. Add that one
          # case to the list rather than retrying everything, so a tool error still fails at
          # once. Selectors in the pipeline's own config stay more specific and still win.
          echo "    errorStrategy = { task.exitStatus == Integer.MAX_VALUE || task.exitStatus in ((130..145) + 104 + 175) ? 'retry' : 'finish' }"
          echo "    maxRetries = 1"
          # Always emitted, with or without an account: it also replaces the request the
          # shipped overlay writes in the other scheduler's syntax, which sbatch rejects.
          if [ -n "~{slurm_account}" ]; then
            echo "    clusterOptions = '--account=~{slurm_account}'"
          else
            echo "    clusterOptions = ''"
          fi
          echo "}"
          if [ "${#bind_list[@]}" -gt 0 ]; then
            binds=$(IFS=,; echo "${bind_list[*]}")
            echo "singularity { runOptions = '-B ${binds}' }"
          fi
        } > scheduler.config
        config_args+=(-c "$(pwd)/scheduler.config")
      fi

      # -ansi-log false because stdout is a Cromwell log file, not a terminal: the ANSI live
      # display rewrites lines and truncates process names to a nominal width, e.g.
      # "NFC...rityEstimateV3_test_01_WG)", which makes the log useless for grepping. Plain
      # mode prints one untruncated line per process event. NOTE the comment has to sit HERE,
      # above the invocation: a # inside the backslash continuation would comment out the
      # rest of the command, and womtool cannot catch that because it does not parse bash.
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
          "${config_args[@]}" \
          "${stub_args[@]}" \
          -ansi-log false \
          -resume

      echo "${abs_outdir}/~{group_id}" > output_dir.txt

      # The WG results are archived by pack_wgts, not here, so that the archive can be built
      # after pre_filtering and carry the prefiltered somatic VCF alongside the original.
      tar -czf ~{outdir}.pipeline_info.tar.gz \
          -C "${abs_outdir}" \
          pipeline_info/
    >>>

    output {
        String output_dir             = read_string("output_dir.txt")
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
        File    longitudinal_bam
        File    longitudinal_bai
        String  longitudinal_sample_id
        Boolean run_redux
        Boolean nextflow_stub
        Boolean use_primary_filters
        String  sequencing_platform
        String  primary_purple_dir
        String  outdir
        String  ref_data_dir
        String  images_dir
        String  pipeline_dir
        String  nextflow_bin
        String  nextflow_home
        String  scheduler
        String  slurm_partition = ""
        String  slurm_account   = ""
        Array[String] singularity_binds = []
        Array[String] nextflow_config
        String  submit_wrapper
        String  modules
        Int memory  = 32
        Int timeout = 10
    }

    parameter_meta {
        group_id:               "Sample group identifier; used as output subdirectory and samplesheet group_id"
        subject_id:             "Subject/patient identifier; prepended to the wisp summary row"
        tumor_sample_id:        "Primary tumour sample ID; must match the purple output filenames in wgts_outdir"
        longitudinal_bam:       "Longitudinal ctDNA BAM (merged if there were several inputs)"
        longitudinal_bai:       "Index for longitudinal_bam"
        longitudinal_sample_id: "Longitudinal sample ID used in the samplesheet"
        run_redux:              "When true REDUX processes the alignments; when false they are declared bam_redux with generate_redux_tsvs_only so REDUX only regenerates its TSVs"
        nextflow_stub:          "Run oncoanalyser with -stub --create_stub_placeholders: placeholder outputs, no real compute"
        sequencing_platform:    "Value for --sequencing_platform: illumina, sbx or ultima. Taken from the longitudinal sample, the only sample this run processes"
        primary_purple_dir:     "PURPLE directory of the primary, holding both the full somatic VCF and the prefiltered one. Which of the two SAGE_APPEND sees is decided by use_primary_filters"
        use_primary_filters:    "When true, stage the prefiltered somatic VCF under the canonical filename so SAGE_APPEND works from the reduced site list. oncoanalyser resolves that VCF by exact name, which is why this is a staged copy rather than a different path"
        outdir:                 "Output directory; the pipeline writes to outdir/group_id/"
        ref_data_dir:           "HMF reference data directory, used for --igenomes_base, --hmf_genomes_base and --ref_data_hmf_data_path; normally the literal $REFERENCE_FILES_DIR"
        images_dir:             "Singularity image cache directory (NXF_SINGULARITY_CACHEDIR); normally the literal $IMAGES_DIR, expanded by the shell after the oncoanalyser module loads"
        pipeline_dir:           "oncoanalyser checkout containing main.nf; normally the literal $ONCOANALYSER_FOLDER"
        nextflow_bin:           "Nextflow executable; defaults to `nextflow` on PATH, with the version pinned by the module's NXF_VER. Must resolve to 25.10.0 or newer"
        nextflow_home:          "NXF_HOME holding the pre-cached nf-schema plugin, so the run works with NXF_OFFLINE=true"
        scheduler:          "Batch scheduler Nextflow submits to; the submit wrapper is installed only for sge"
        slurm_partition:      "Partition for the jobs nextflow submits, when scheduler is slurm"
        slurm_account:        "Accounting group for those jobs; also clears the request the module's overlay writes in the other scheduler's syntax"
        singularity_binds:    "Paths bound into every container, replacing the bind the module's overlay sets. Empty keeps that bind"
        nextflow_config:      "Config overlays, each passed with its own -c, in the order given"
        submit_wrapper:         "Wrapper placed on PATH ahead of the scheduler submit command, to adjust the resource requests Nextflow generates"
        modules:                "Environment modules to load"
        memory:                 "Memory in GB for this task, which hosts the Nextflow driver only; the pipeline processes get their own allocations"
        timeout:                "Wall-clock timeout in hours"
    }

    command <<<
      set -euo pipefail
      WORKDIR=$(pwd)

      ln -s "~{longitudinal_bam}" "${WORKDIR}/~{longitudinal_sample_id}.bam"
      ln -s "~{longitudinal_bai}" "${WORKDIR}/~{longitudinal_sample_id}.bam.bai"

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

      # No AMBER: it is only of use for the LOH contribution, which cannot currently be
      # obtained (see the samplesheet comment below), and running it costs a second REDUX for
      # the normal plus AMBER itself for nothing.
      processes="redux,cobalt,sage_append,wisp"

      # The primary normal is deliberately NOT part of this run. It would only serve WISP's
      # AMBER_LOH method, and that method cannot engage as things stand: wisp_analysis derives
      # primary_normal_aln from REDUX_DIR_NORMAL, so the normal has to be supplied as a
      # redux_dir, not as an alignment. Supplying it as a BAM passes samplesheet validation
      # (Utils.groovy accepts hasNormalDnaBam OR hasNormalDnaReduxInput) and then silently
      # yields no LOH -- while costing a second REDUX for the normal plus AMBER. The SNV
      # result is unaffected either way, because SAGE_APPEND is given only the longitudinal
      # REDUX output.
      #
      # amber_dir / cobalt_dir / sage_append_dir must never appear on the longitudinal
      # sample: oncoanalyser treats that as a fatal input clash.
      # Which of the primary's two somatic VCFs SAGE_APPEND works from. oncoanalyser resolves
      # it by EXACT filename -- file(purple_dir).resolve("<id>.purple.somatic.vcf.gz") -- so
      # using the prefiltered one means staging a copy of the directory in which that name
      # holds the prefiltered content. No filtering happens here; the WG stage wrote both.
      purple_dir="~{primary_purple_dir}"
      if ~{use_primary_filters}; then
        prefiltered="${purple_dir}/~{tumor_sample_id}.purple.somatic.prefiltered.vcf.gz"
        if [ -s "${prefiltered}" ]; then
          staged="${WORKDIR}/purple_staged"
          mkdir -p "${staged}"
          cp -r "${purple_dir}/." "${staged}/"
          cp "${prefiltered}" "${staged}/~{tumor_sample_id}.purple.somatic.vcf.gz"
          if [ -f "${prefiltered}.tbi" ]; then
            cp "${prefiltered}.tbi" "${staged}/~{tumor_sample_id}.purple.somatic.vcf.gz.tbi"
          fi
          purple_dir="${staged}"
          echo "using the prefiltered primary VCF" >&2
        else
          # A WG tarball made before pre-filtering existed has no prefiltered VCF. Warn and
          # carry on with the full call set rather than failing a run that is still valid.
          echo "WARNING: no ~{tumor_sample_id}.purple.somatic.prefiltered.vcf.gz in" >&2
          echo "         ${purple_dir}; using the full call set instead." >&2
        fi
      fi

      {
        echo "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath"
        echo "~{group_id},~{subject_id},~{tumor_sample_id},tumor,dna,purple_dir,,${purple_dir}"
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
      # -XX:-UseContainerSupport: the JVM's cgroup probe aborts startup on a node where it
      # cannot enumerate the controllers. The heap is pinned here, so nothing depends on the
      # limits that probe would report.
      export NXF_OPTS="-Xms512m -Xmx8g -XX:-UseContainerSupport"
      export NXF_SINGULARITY_CACHEDIR=~{images_dir}
      # NXF_HOME must hold the pre-cached plugins; with NXF_OFFLINE=true a wrong one fails
      # obscurely inside plugin resolution. Prefer the module's value, falling back to a
      # sibling directory that does have plugins/.
      nxf_home="~{nextflow_home}"
      if [ ! -d "${nxf_home}/plugins" ]; then
        for cand in "${ONCOANALYSER_ROOT:-}"/nextflow_home*; do
          if [ -d "${cand}/plugins" ]; then
            echo "WARNING: NEXTFLOW_HOME (${nxf_home}) has no plugins/; using ${cand}" >&2
            echo "         Fix NEXTFLOW_HOME in the oncoanalyser modulefile." >&2
            nxf_home="${cand}"
            break
          fi
        done
      fi
      [ -d "${nxf_home}/plugins" ] || {
        echo "ERROR: no NXF_HOME with a plugins/ directory. Tried ${nxf_home} and" >&2
        echo "       ${ONCOANALYSER_ROOT:-<ONCOANALYSER_ROOT unset>}/nextflow_home*" >&2
        exit 1
      }
      # NXF_HOME must be WRITABLE. The nextflow launcher does mkdir -p $NXF_HOME/tmp and the
      # local secrets provider wants $NXF_HOME/secrets, but a module tree is read-only, so
      # pointing NXF_HOME straight at it fails with "Read-only file system". Give nextflow a
      # writable directory in the task workdir. plugins/ itself must be a real directory:
      # nextflow calls createDirectories on it, which rejects a symlink even when it resolves
      # to a directory. Its entries are only read, so they can be symlinks into the module.
      nxf_home_rw="$(pwd)/nxf_home"
      mkdir -p "${nxf_home_rw}/plugins"
      for plugin in "${nxf_home}"/plugins/*; do
        if [ -e "${plugin}" ]; then
          ln -sfn "${plugin}" "${nxf_home_rw}/plugins/"
        fi
      done
      export NXF_HOME="${nxf_home_rw}"
      unset NXF_DIST NXF_LAUNCHER NXF_PLUGINS_DIR || true

      # The SGE executor embeds h_rss/mem_free directly in .command.run, where configuration
      # cannot reach them, so a shim ahead of qsub on PATH rewrites the request. The Slurm
      # executor emits --mem and --cpus-per-task from the memory and cpus directives, so it
      # needs no wrapper.
      if [ "~{scheduler}" = "sge" ]; then
        bin_dir="$(pwd)/bin"
        mkdir -p "${bin_dir}"
        cp ~{submit_wrapper} "${bin_dir}/qsub"
        chmod +x "${bin_dir}/qsub"
        export PATH="${bin_dir}:$PATH"
      fi

      # One -c per overlay, in the order given: later files win, which is how a site config
      # overrides the shared one.
      configs=(~{sep=" " nextflow_config})
      config_args=()
      for cfg in "${configs[@]}"; do
        config_args+=(-c "${cfg}")
      done

      # The overlay the module ships selects its own scheduler, its queue and its container
      # binds. Generate the overrides here instead of requiring a file on the filesystem, so
      # the workflow carries everything it needs, and apply them last so they win. Only the
      # values that vary by site are substituted; the shape is the same everywhere.
      if [ "~{scheduler}" = "slurm" ]; then
        bind_list=(~{sep=" " singularity_binds})
        {
          echo "executor { name = 'slurm' }"
          echo "process {"
          echo "    queue = '~{slurm_partition}'"
          # Slurm enforces the memory request as RSS. A JVM sized at the module default of
          # ~95% of it leaves nothing for the helper processes the tools fork, and the cgroup
          # kills the job. Keep a quarter of the request outside the heap. Processes that do
          # not read this setting are unaffected.
          echo "    ext.xmx_mod = 0.75"
          # A cgroup kill reports no exit status, which the pipeline's own list does not
          # match, so the job would not be retried at all -- and the per-process memory
          # directives scale with task.attempt, which is otherwise unreachable. Add that one
          # case to the list rather than retrying everything, so a tool error still fails at
          # once. Selectors in the pipeline's own config stay more specific and still win.
          echo "    errorStrategy = { task.exitStatus == Integer.MAX_VALUE || task.exitStatus in ((130..145) + 104 + 175) ? 'retry' : 'finish' }"
          echo "    maxRetries = 1"
          # Always emitted, with or without an account: it also replaces the request the
          # shipped overlay writes in the other scheduler's syntax, which sbatch rejects.
          if [ -n "~{slurm_account}" ]; then
            echo "    clusterOptions = '--account=~{slurm_account}'"
          else
            echo "    clusterOptions = ''"
          fi
          echo "}"
          if [ "${#bind_list[@]}" -gt 0 ]; then
            binds=$(IFS=,; echo "${bind_list[*]}")
            echo "singularity { runOptions = '-B ${binds}' }"
          fi
        } > scheduler.config
        config_args+=(-c "$(pwd)/scheduler.config")
      fi

      # -ansi-log false because stdout is a Cromwell log file, not a terminal: the ANSI live
      # display rewrites lines and truncates process names to a nominal width, e.g.
      # "NFC...rityEstimateV3_test_01_WG)", which makes the log useless for grepping. Plain
      # mode prints one untruncated line per process event. NOTE the comment has to sit HERE,
      # above the invocation: a # inside the backslash continuation would comment out the
      # rest of the command, and womtool cannot catch that because it does not parse bash.
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
          "${config_args[@]}" \
          "${stub_args[@]}" \
          -ansi-log false \
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
