process CREATE_COLOR_PALETTES {
    label 'mbv_r'
    container r_image

    input:
    path palette_script
    val r_image

    output:
    path 'colorPalettes.rds', emit: palettes

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$palette_script" colorPalettes.rds
    test -s colorPalettes.rds
    """
}
