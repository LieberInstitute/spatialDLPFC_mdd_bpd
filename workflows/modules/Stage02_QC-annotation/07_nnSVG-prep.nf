process NNSVG_PREP {
    label 'mbv_r'
    container r_image

    input:
    path prep_script
    val r_image
    tuple path(assay_h5), path(se_rds)

    output:
    path 'per-sample_spe/*', type: 'dir', emit: sample_spe
    path 'nnSVG_run_specifics.csv', emit: run_specifics

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$prep_script" "$assay_h5"
    """
}
