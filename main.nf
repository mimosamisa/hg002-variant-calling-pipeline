params.reads     = "$projectDir/data/raw/HG002.fastq"
params.ref       = "$projectDir/data/reference/GRCh38.primary_assembly.genome.fa"
params.ref_fai   = "$projectDir/data/reference/GRCh38.primary_assembly.genome.fa.fai"
params.truth_vcf = "$projectDir/data/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.truth_tbi = "$projectDir/data/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
params.truth_bed = "$projectDir/data/giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

process ALIGN {
    tag "Align HG002"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    path reads
    path ref

    output:
    path "HG002.sam"

    script:
    """
    minimap2 -ax map-hifi -t ${task.cpus} ${ref} ${reads} > HG002.sam
    """
}

process SORT_INDEX {
    tag "Sort and Index BAM"
    publishDir "results/alignment", mode: 'copy'
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    input:
    path sam

    output:
    tuple path("HG002.sorted.bam"), path("HG002.sorted.bam.bai")

    script:
    """
    samtools sort -@ ${task.cpus} -o HG002.sorted.bam ${sam}
    samtools index HG002.sorted.bam
    """
}

process CLAIR3 {
    tag "Call variants with Clair3"
    publishDir "results/clair3", mode: 'copy'
    container "hkubal/clair3:latest"

    input:
    tuple path(bam), path(bai)
    path ref
    path ref_fai

    output:
    path "clair3_output/merge_output.vcf.gz"

    script:
    """
    run_clair3.sh \
      --bam_fn=${bam} \
      --ref_fn=${ref} \
      --threads=${task.cpus} \
      --platform=hifi \
      --model_path=/opt/models/hifi \
      --output=clair3_output \
      --include_all_ctgs
    """
}

process DEEPVARIANT {
    tag "Call variants with DeepVariant"
    publishDir "results/deepvariant", mode: 'copy'
    container "google/deepvariant:1.6.0"

    input:
    tuple path(bam), path(bai)
    path ref
    path ref_fai

    output:
    tuple path("HG002_deepvariant.vcf.gz"), path("HG002_deepvariant.g.vcf.gz")

    script:
    """
    run_deepvariant \
      --model_type=PACBIO \
      --ref=${ref} \
      --reads=${bam} \
      --num_shards=${task.cpus} \
      --output_vcf=HG002_deepvariant.vcf.gz \
      --output_gvcf=HG002_deepvariant.g.vcf.gz
    """
}

process HAPY_CLAIR3 {
    tag "Benchmark Clair3 with hap.py"
    publishDir "results/benchmark", mode: 'copy'
    container "quay.io/biocontainers/hap.py:0.3.14--py27h5c5a3ab_0"

    input:
    path vcf
    path ref
    path ref_fai
    path truth_vcf
    path truth_tbi
    path truth_bed

    output:
    path "benchmark_clair3.*"

    script:
    """
    hap.py \
        ${truth_vcf} \
        ${vcf} \
        -f ${truth_bed} \
        -r ${ref} \
        -o benchmark_clair3 \
        --engine=xcmp \
        --threads ${task.cpus}
    """
}

process HAPY_DV {
    tag "Benchmark DeepVariant with hap.py"
    publishDir "results/benchmark", mode: 'copy'
    container "quay.io/biocontainers/hap.py:0.3.14--py27h5c5a3ab_0"

    input:
    path vcf
    path ref
    path ref_fai
    path truth_vcf
    path truth_tbi
    path truth_bed

    output:
    path "benchmark_deepvariant.*"

    script:
    """
    hap.py \
        ${truth_vcf} \
        ${vcf} \
        -f ${truth_bed} \
        -r ${ref} \
        -o benchmark_deepvariant \
        --engine=xcmp \
        --threads ${task.cpus}
    """
}

workflow {
    reads_ch     = Channel.fromPath(params.reads)
    ref_ch       = Channel.fromPath(params.ref)
    ref_fai_ch   = Channel.fromPath(params.ref_fai)
    truth_vcf_ch = Channel.fromPath(params.truth_vcf)
    truth_tbi_ch = Channel.fromPath(params.truth_tbi)
    truth_bed_ch = Channel.fromPath(params.truth_bed)

    sam_ch     = ALIGN(reads_ch, ref_ch)
    bam_ch     = SORT_INDEX(sam_ch)
    clair3_vcf = CLAIR3(bam_ch, ref_ch, ref_fai_ch)
    dv_vcf     = DEEPVARIANT(bam_ch, ref_ch, ref_fai_ch)

    HAPY_CLAIR3(clair3_vcf, ref_ch, ref_fai_ch, truth_vcf_ch, truth_tbi_ch, truth_bed_ch)
    HAPY_DV(dv_vcf.map { vcf, gvcf -> vcf }, ref_ch, ref_fai_ch, truth_vcf_ch, truth_tbi_ch, truth_bed_ch)
}

