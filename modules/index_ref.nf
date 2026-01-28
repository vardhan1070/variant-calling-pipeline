process INDEX_REF {

    tag "index_ref"

    publishDir "results/reference", mode: 'copy'

    input:
    path ref

    output:
    tuple path(ref),
          path("${ref}.amb"),
          path("${ref}.ann"),
          path("${ref}.bwt"),
          path("${ref}.pac"),
          path("${ref}.sa")

    script:
    """
    bwa index ${ref}
    """
}
