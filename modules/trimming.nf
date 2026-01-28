process TRIMMING {

    tag "$sample_id"

    publishDir "results/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("R1.trimmed.fastq"),
          path("R2.trimmed.fastq")

    script:
    """
    fastp \
        -i ${r1} \
        -I ${r2} \
        -o R1.trimmed.fastq \
        -O R2.trimmed.fastq
    """
}
