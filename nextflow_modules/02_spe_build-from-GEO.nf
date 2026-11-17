process SPE_GEO_BUILD {
    container 'mbv-r:test'

    input:
    path spe_script
    val n_samples
    path geo_bundle
    path sample_set

    output:
    tuple path("spe_n${n_samples}_assays.h5"), path("spe_n${n_samples}_se.rds"), emit: spe_files

    script:
    """
    export TMPDIR=\$PWD # Force R's tempdir() to use the current Nextflow task work directory
    Rscript --vanilla "$spe_script" $geo_bundle $sample_set
    test -s "spe_n${n_samples}_assays.h5"
    test -s "spe_n${n_samples}_se.rds"
    """
}
