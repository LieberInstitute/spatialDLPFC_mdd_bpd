process SPOTSWEEPER {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    val r_image
    path spe_dir
    val n_samples
    path qc_df, stageAs: 'colData_MAD.csv'
    val ignore_edge_samples

    output:
    path 'colData_QC.csv', emit: coldata

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$spe_dir" "$n_samples" "$qc_df" "$ignore_edge_samples"
    test -s colData_QC.csv
    """
}
