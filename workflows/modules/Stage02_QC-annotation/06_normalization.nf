process NORMALIZATION {
    label 'mbv_r'
    container r_image

    input:
    path norm_script
    val r_image
    val n_samples
    tuple path(assay_h5), path(se_rds)
    path qc_df

    output:
    tuple path("spe_n${n_samples}_postQC_norm_assays.h5"), path("spe_n${n_samples}_postQC_norm_se.rds"), emit: spe_files

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$norm_script" "$assay_h5" "$qc_df"
    test -s "spe_n${n_samples}_postQC_norm_assays.h5"
    test -s "spe_n${n_samples}_postQC_norm_se.rds"
    """
}
