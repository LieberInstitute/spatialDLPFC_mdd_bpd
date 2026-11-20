process REVIEW_EDGES_OUTLIER_AREAS {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    val r_image
    path spe_dir
    val n_samples
    path qc_coldata

    output:
    path 'check_edges-outlier-areas.pdf', emit: review_pdf

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$spe_dir" "$n_samples" "$qc_coldata"
    test -s check_edges-outlier-areas.pdf
    """
}
