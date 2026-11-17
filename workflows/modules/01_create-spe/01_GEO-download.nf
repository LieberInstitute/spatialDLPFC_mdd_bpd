process FETCH_GEO {
    container 'mbv-r:test'

    input:
    path spe_script
    val n_samples
    val timeout

    output:
    path 'GSE333974', emit: geo_bundle
    path 'GEO_sample_ids.txt', emit: geo_samples

    script:
    """
    export TMPDIR=\$PWD # Force R's tempdir() to use the current Nextflow task work directory
    Rscript --vanilla "$spe_script" $n_samples $timeout
    """
}
