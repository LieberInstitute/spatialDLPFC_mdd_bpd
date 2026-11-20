process SPOTSWEEPER {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    val r_image
    tuple path(assay_h5), path(se_rds)
    path qc_df, stageAs: 'colData_MAD.csv' // stage as unique name so symlink isn't overwritten
    val ignore_edge_samples

    output:
    path 'colData_QC.csv', emit: coldata

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$assay_h5" "$qc_df" "$ignore_edge_samples"
    test -s colData_QC.csv
    """
}
