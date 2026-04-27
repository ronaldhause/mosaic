#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// build_reference: regenerate the microsatellite bed used by main.nf
// from an hg19 fasta and a capture-bait bed, following the paper's
// methods (misa for tract discovery -> pad +/- 5 bp -> merge compound/
// complex with <= 10 bp gap -> retain loci within 50 bp of capture baits
// -> optional annovar annotation).
//
// outputs: microsatellites.bed, locus_annotation.tsv, build_summary.tsv

def help_message() {
    log.info """
    usage:
      nextflow run pipeline/build_reference.nf \\
        --reference hg19.fa \\
        --capture_bait_bed baits.bed \\
        --outdir ref_build \\
        -profile local,docker

    required:
      --reference         hg19 fasta (a matching .fai must sit next to it)
      --capture_bait_bed  capture kit bait bed (e.g. nimblegen seqcap_ez_exome_v3)

    optional:
      --outdir            output dir [default: ./ref_build]
      --misa_motifs       misa motif spec [default: 1-5 2-5 3-5 4-5 5-5]
      --misa_interrupts   misa interrupts param [default: 0]
      --pad_bp            +/- bp padding around each tract [default: 5]
      --merge_gap_bp      max gap (bp) to merge adjacent tracts [default: 10]
      --bait_slop_bp      bait-proximity slop applied before intersect [default: 50]
      --annovar_dir       path to an annovar install (has annotate_variation.pl) [default: null]
      --annovar_db        annovar db dir (e.g. humandb) [default: null]
      --annovar_build     annovar genome build [default: hg19]

    profiles: local, slurm, docker, singularity
    """.stripIndent()
}

if (params.help) {
    help_message()
    exit 0
}


process validate_inputs {
    label 'process_low'

    input:
    path reference
    path fai
    path bait_bed

    output:
    path "validate.ok"

    script:
    """
    python3 ${projectDir}/bin/validate_build_inputs.py \\
        --reference ${reference} \\
        --fai ${fai} \\
        --capture_bait_bed ${bait_bed}
    touch validate.ok
    """
}


// misa runs a perl script against the full fasta, writing <input>.misa next
// to the input. config is an ini file in cwd; we generate it inline.
process run_misa {
    label 'process_medium'

    input:
    path reference
    val misa_motifs
    val misa_interrupts
    path validation_ok

    output:
    path "misa_output.misa"

    script:
    """
cat > misa.ini <<EOF
definition(unit_size,min_repeats):                   ${misa_motifs}
interruptions(max_difference_between_2_SSRs):        ${misa_interrupts}
GFF:                                                 false
EOF

    # misa writes <input>.misa alongside the input. stage a local copy so
    # the output lands in the work dir.
    ln -sf ${reference} ref.fa
    perl \$(command -v misa.pl || echo /opt/misa/misa.pl) ref.fa
    mv ref.fa.misa misa_output.misa
    """
}


process misa_to_bed {
    label 'process_low'

    input:
    path misa_output
    path fai
    val pad_bp

    output:
    path "loci_raw.tsv"

    script:
    """
    python3 ${projectDir}/bin/misa_to_bed.py \\
        --misa ${misa_output} \\
        --fai ${fai} \\
        --pad ${pad_bp} \\
        --output loci_raw.tsv
    """
}


process merge_compound_complex {
    label 'process_low'

    input:
    path loci_raw
    val merge_gap_bp

    output:
    tuple path("loci_merged.tsv"), path("loci_merged.bed")

    script:
    """
    python3 ${projectDir}/bin/merge_compound_complex.py \\
        --input ${loci_raw} \\
        --max_gap ${merge_gap_bp} \\
        --output_tsv loci_merged.tsv \\
        --output_bed loci_merged.bed
    """
}


process filter_by_baits {
    label 'process_low'

    input:
    path loci_merged_bed
    path bait_bed
    path fai
    val bait_slop_bp

    output:
    path "microsatellites.bed"

    script:
    """
    # genome file for bedtools slop (chrom sizes from fai)
    cut -f1,2 ${fai} > chrom.sizes

    # normalise bait bed to 3 cols and slop by bait_slop_bp on each side,
    # clamped to chrom bounds.
    awk 'BEGIN{OFS="\\t"} /^(#|track|browser)/ {next} {print \$1,\$2,\$3}' ${bait_bed} \\
        | sort -k1,1 -k2,2n > baits.sorted.bed

    bedtools slop -b ${bait_slop_bp} -i baits.sorted.bed -g chrom.sizes \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i - > baits.slopped.bed

    sort -k1,1 -k2,2n ${loci_merged_bed} > loci.sorted.bed
    bedtools intersect -u -a loci.sorted.bed -b baits.slopped.bed \\
        > microsatellites.bed
    """
}


process annotate_annovar {
    label 'process_low'

    input:
    path ms_bed
    val annovar_dir
    val annovar_db
    val annovar_build

    output:
    path "annovar_annotation.tsv"

    when:
    annovar_dir && annovar_db

    script:
    """
    # build an avinput from the bed: chrom, start, end, 0, 0
    awk 'BEGIN{OFS="\\t"} {print \$1,\$2+1,\$3,"0","0"}' ${ms_bed} > loci.avinput

    ${annovar_dir}/annotate_variation.pl \\
        --geneanno --dbtype refGene \\
        --buildver ${annovar_build} \\
        loci.avinput \\
        ${annovar_db}

    # variant_function output: func, gene, chrom, start, end, ref, alt
    awk 'BEGIN{OFS="\\t"; print "chrom","start","end","genomic_class","gene"} \\
         {print \$3,\$4-1,\$5,\$1,\$2}' \\
        loci.avinput.variant_function > annovar_annotation.tsv
    """
}


process emit_annotation {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path loci_merged_tsv
    path ms_bed
    path annovar_tsv

    output:
    path "locus_annotation.tsv"

    script:
    def annovar_arg = annovar_tsv.name != 'NO_FILE' ? "--annovar_tsv ${annovar_tsv}" : ''
    """
    python3 ${projectDir}/bin/emit_locus_annotation.py \\
        --merged_tsv ${loci_merged_tsv} \\
        --ms_bed ${ms_bed} \\
        ${annovar_arg} \\
        --output locus_annotation.tsv
    """
}


process emit_summary {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path loci_raw
    path loci_merged_tsv
    path ms_bed
    val defb_locus

    output:
    path "build_summary.tsv"
    path "microsatellites.bed"

    script:
    """
    cp ${ms_bed} microsatellites.bed

    python3 ${projectDir}/bin/emit_build_summary.py \\
        --loci_raw ${loci_raw} \\
        --merged_tsv ${loci_merged_tsv} \\
        --ms_bed ${ms_bed} \\
        --defb_locus "${defb_locus}" \\
        --output build_summary.tsv
    """
}


workflow {
    if (!params.reference)          { exit 1, "missing --reference" }
    if (!params.capture_bait_bed)   { exit 1, "missing --capture_bait_bed" }

    ref      = file(params.reference)
    fai      = file("${params.reference}.fai")
    baits    = file(params.capture_bait_bed)
    outdir   = params.outdir

    ok = validate_inputs(ref, fai, baits)

    misa_out   = run_misa(ref, params.misa_motifs, params.misa_interrupts, ok)
    raw_loci   = misa_to_bed(misa_out, fai, params.pad_bp)
    merged     = merge_compound_complex(raw_loci, params.merge_gap_bp)

    merged_tsv = merged.map { tsv, bed -> tsv }
    merged_bed = merged.map { tsv, bed -> bed }

    ms_bed = filter_by_baits(merged_bed, baits, fai, params.bait_slop_bp)

    // optional annovar — emit a sentinel "NO_FILE" when disabled so the
    // downstream process knows to skip annotation.
    if (params.annovar_dir && params.annovar_db) {
        annovar_tsv = annotate_annovar(
            ms_bed, params.annovar_dir, params.annovar_db, params.annovar_build
        )
    } else {
        annovar_tsv = Channel.of(file("${projectDir}/assets/NO_FILE"))
    }

    emit_annotation(merged_tsv, ms_bed, annovar_tsv)
    emit_summary(raw_loci, merged_tsv, ms_bed, params.defb_locus)
}
