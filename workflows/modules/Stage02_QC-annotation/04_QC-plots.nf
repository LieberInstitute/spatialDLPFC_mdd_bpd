process QC_PLOTS {
    label 'mbv_r'
    container r_image

    input:
    path qc_script
    val r_image
    path qc_df
    path color_palettes

    output:
    path 'qc_summary_plots.pdf', emit: summary_pdf

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$qc_script" "$qc_df" "$color_palettes"
    test -s qc_summary_plots.pdf
    """
}
