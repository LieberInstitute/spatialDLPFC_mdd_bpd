process ACTUALIZE_QC_FILTERS {
    label 'mbv_r'
    container r_image

    // quickResave must update this task's RDS, not an upstream symlink.
    stageInMode 'copy'

    input:
    path qc_script
    val r_image
    val n_samples
    tuple path(assay_h5), path(se_rds)
    path qc_df, stageAs: 'colData_generated_QC.csv' // stage as unique name so symlink isn't overwritten
    val keep_generated_qc

    output:
    tuple path("spe_n${n_samples}_assays.h5"), path("spe_n${n_samples}_se.rds"), emit: spe_files
    path 'colData_QC.csv', emit: coldata
    path 'QC_filter_comparison.csv', emit: comparison
    path 'QC_filter_differences.csv', emit: differences

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$assay_h5" "$qc_df" "$keep_generated_qc"
    test -s "spe_n${n_samples}_assays.h5"
    test -s "spe_n${n_samples}_se.rds"
    test -s colData_QC.csv
    """
}
