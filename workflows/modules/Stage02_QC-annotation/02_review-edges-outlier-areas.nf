process REVIEW_EDGES_OUTLIER_AREAS {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    val r_image
    tuple path(assay_h5), path(se_rds)
    path qc_df

    output:
    path 'check_edges-outlier-areas.pdf', emit: review_pdf

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$assay_h5" "$qc_df"
    test -s check_edges-outlier-areas.pdf
    """
}
