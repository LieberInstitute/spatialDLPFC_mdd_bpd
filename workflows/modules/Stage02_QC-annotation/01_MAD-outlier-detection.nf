process MAD_OUTLIER_DETECTION {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    path edge_functions
    val r_image
    tuple path(assay_h5), path(se_rds)

    output:
    path 'colData_QC.csv', emit: coldata

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$assay_h5" "$edge_functions"
    test -s colData_QC.csv
    """
}
