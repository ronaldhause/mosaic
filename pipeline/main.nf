#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// mosaic: msi-h/mss classification from paired tumor/normal exome bams.
// hause et al., nature medicine 2016. see pipeline/README or repo root.

def help_message() {
    log.info """
    usage:
      nextflow run pipeline/main.nf \\
        --samplesheet samples.csv \\
        --reference hg19.fa \\
        --microsatellite_bed loci.bed \\
        -profile local

    required:
      --samplesheet           csv with columns: sample_name,tumor_bam,normal_bam,tumor_type
      --reference             hg19 fasta (indexed)
      --microsatellite_bed    bed of microsatellite loci

    optional:
      --outdir                    output dir (default: ./results)
      --num_loci                  total loci in bed (default: 516876)
      --defb_locus                defb105a/b coordinates (default: 8:7679723-7679741)
      --peak_avg_threshold_high   msi-h threshold (default: 0.0055)
      --peak_avg_threshold_low    defb-conditional threshold (default: 0.0029)

    profiles: local, slurm, docker, singularity
    """.stripIndent()
}

if (params.help) {
    help_message()
    exit 0
}

// --- process definitions ---

// run_msings: per sample-pair. implements the paper's high-sensitivity
// paired approach — msings analyzer on each bam separately, then python
// compares the per-locus length distributions (see bin/compare_paired_msings.py).
//
// msings versions vary on column names and exact output formatting; the
// awk reshape below pulls locus + length-count blob into a 2-col tsv that
// compare_paired_msings.py auto-detects. if msi is installed under a
// different entrypoint (e.g. msings/run_msings), set params.msings_cmd.
process run_msings {
    tag "$sample_name"
    label 'process_medium'

    input:
    tuple val(sample_name), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(tumor_type)
    path reference
    path reference_fai
    path microsatellite_bed

    output:
    tuple val(sample_name), val(tumor_type), path("${sample_name}.msings.tsv")

    script:
    def msings_cmd = params.msings_cmd ?: 'msi analyzer'
    """
    # per-sample analyzer outputs. produce one row per locus with a trailing
    # column of 'offset:depth:supporting_calls' pairs (msings default).
    for kind in tumor normal; do
        if [ "\$kind" = tumor ]; then
            bam=${tumor_bam}
        else
            bam=${normal_bam}
        fi
        ${msings_cmd} \\
            \$bam \\
            -l ${microsatellite_bed} \\
            -r ${reference} \\
            -o \${kind}.analyzer.txt
    done

    # reshape to 'locus<TAB>offset:count;...' then compute peak_diff.
    python3 ${projectDir}/bin/reshape_msings_analyzer.py \\
        --input tumor.analyzer.txt  --output tumor.tall.tsv
    python3 ${projectDir}/bin/reshape_msings_analyzer.py \\
        --input normal.analyzer.txt --output normal.tall.tsv

    python3 ${projectDir}/bin/compare_paired_msings.py \\
        --tumor tumor.tall.tsv \\
        --normal normal.tall.tsv \\
        --sample_name ${sample_name} \\
        --min_rel_abundance ${params.min_rel_abundance} \\
        --output ${sample_name}.msings.tsv
    """
}

// compute_features: per sample. wraps bin/compute_features.R.
process compute_features {
    tag "$sample_name"
    label 'process_low'

    input:
    tuple val(sample_name), val(tumor_type), path(msings_tsv)

    output:
    path "${sample_name}.features.csv"

    script:
    """
    Rscript ${projectDir}/bin/compute_features.R \\
        --input ${msings_tsv} \\
        --sample_name ${sample_name} \\
        --tumor_type ${tumor_type} \\
        --num_loci ${params.num_loci} \\
        --defb_locus "${params.defb_locus}" \\
        --output ${sample_name}.features.csv
    """
}

// classify_samples: gather step. runs once over all per-sample features.
process classify_samples {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path feature_files

    output:
    path "mosaic_results.csv"

    script:
    """
    mkdir -p features
    # stage inputs into a single dir for the pattern glob
    for f in ${feature_files}; do cp -L "\$f" features/; done

    Rscript ${projectDir}/bin/classify_msi.R \\
        --input_dir features \\
        --pattern "*.features.csv" \\
        --threshold_high ${params.peak_avg_threshold_high} \\
        --threshold_low ${params.peak_avg_threshold_low} \\
        --output mosaic_results.csv
    """
}

// resolves bai index, checking both <bam>.bai and <bam>.bai (s/.bam$/.bai/) conventions.
def resolve_bai(bam_path) {
    def bam = file(bam_path)
    def bai_a = file("${bam}.bai")
    def bai_b = file(bam.toString().replaceAll(/\.bam$/, '.bai'))
    if (bai_a.exists()) return bai_a
    if (bai_b.exists()) return bai_b
    // fall back to .bai sibling; nextflow will error at stage time if missing
    return bai_a
}

workflow {
    // required params check
    if (!params.samplesheet)         { exit 1, "missing --samplesheet" }
    if (!params.reference)           { exit 1, "missing --reference" }
    if (!params.microsatellite_bed)  { exit 1, "missing --microsatellite_bed" }

    ref     = file(params.reference)
    ref_fai = file("${params.reference}.fai")
    ms_bed  = file(params.microsatellite_bed)

    samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_name,
                file(row.tumor_bam),
                resolve_bai(row.tumor_bam),
                file(row.normal_bam),
                resolve_bai(row.normal_bam),
                row.tumor_type ?: 'NA'
            )
        }

    msings_out   = run_msings(samples, ref, ref_fai, ms_bed)
    features_out = compute_features(msings_out)
    classify_samples(features_out.collect())
}
