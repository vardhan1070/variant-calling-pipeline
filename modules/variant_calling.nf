process VARIANT_CALLING {

    tag "$sample_id"

    publishDir "results/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    tuple path(ref), path(amb), path(ann), path(bwt), path(pac), path(sa)
    path varscan

    output:
    path "${sample_id}.vcf"

    script:
    """
    samtools mpileup -f ${ref} ${bam} > ${sample_id}.mpileup

    java -jar ${varscan} \
        mpileup2vcf ${sample_id}.mpileup \
        --output-vcf 1 > ${sample_id}.vcf
    """
}
