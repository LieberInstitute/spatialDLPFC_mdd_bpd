process FETCH_GEO {
    label 'mbv_r'
    container r_image

    input:
    path spe_script
    val r_image
    val n_samples
    val timeout
    path raw_archive

    output:
    path 'GSE333974', emit: geo_bundle
    path 'GEO_sample_ids.txt', emit: geo_samples

    script:
    """
    export TMPDIR=\$PWD # Force R's tempdir() to use the current Nextflow task work directory
    Rscript --vanilla "$spe_script" "$n_samples" "$timeout" "$raw_archive"
    """
}
