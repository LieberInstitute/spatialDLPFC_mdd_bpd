process GENE_REF {
    label 'mbv_r'
    container r_image
    
    stageInMode 'copy'
    
    input:
    path spe_script
    val r_image
    val n_samples
    val timeout
    tuple path(assay_h5), path(se_rds)
    val reference_gtf
    path staged_gtf

    output:
    tuple path("spe_n${n_samples}_assays.h5"), path("spe_n${n_samples}_se.rds"), emit: spe_files

    script:
    """
    export TMPDIR=\$PWD # Force R's tempdir() to use the current Nextflow task work directory
    Rscript --vanilla "$spe_script" "$timeout" "$assay_h5" "${staged_gtf ? staged_gtf.name : reference_gtf}"
    test -s "spe_n${n_samples}_assays.h5"
    test -s "spe_n${n_samples}_se.rds"
    """
}
