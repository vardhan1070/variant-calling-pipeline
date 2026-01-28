process ALIGNMENT {

    tag "$sample_id"

    publishDir "results/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    tuple path(ref), path(amb), path(ann), path(bwt), path(pac), path(sa)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    bwa mem ${ref} ${r1} ${r2} | \
    samtools sort -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}
