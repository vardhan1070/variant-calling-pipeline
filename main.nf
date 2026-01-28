nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc.nf'
include { TRIMMING } from './modules/trimming.nf'
include { INDEX_REF } from './modules/index_ref.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { VARIANT_CALLING } from './modules/variant_calling.nf'

workflow {

    reads_ch = Channel.fromFilePairs(
        "data/samples/*_{1,2}.fastq",
        flat: true
    )

    ref_ch = Channel.fromPath("data/reference/*.fasta")

    varscan_jar = Channel.fromPath("tools/VarScan.v2.4.6.jar")

    indexed_ref = INDEX_REF(ref_ch)

    FASTQC(reads_ch)

    trimmed_reads = TRIMMING(reads_ch)

    bam_files = ALIGNMENT(trimmed_reads, indexed_ref)

    VARIANT_CALLING(bam_files, indexed_ref, varscan_jar)
}
