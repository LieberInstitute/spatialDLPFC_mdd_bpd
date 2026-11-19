process MAD_OUTLIER_DETECTION {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    path edge_functions
    val r_image
    path spe_dir
    val n_samples

    output:
    path 'colData_edges-outlier-areas.csv', emit: coldata

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$spe_dir" "$n_samples" "$edge_functions"
    test -s colData_edges-outlier-areas.csv
    """
}
