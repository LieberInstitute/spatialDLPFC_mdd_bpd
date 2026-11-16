process SPE_CONTAINER_CHECK {
    container 'mbv-r:test'

    input:
    path spe_script
    val sample_set

    output:
    path 'GEO_donor_data.csv', emit: geo_samples
    tuple val(sample_set), path("02_build_spe/spe_n${sample_set}_assays.h5"), path("02_build_spe/spe_n${sample_set}_se.rds"), emit: spe_files

    script:
    """
    export TMPDIR=\$PWD # Force R's tempdir() to use the current Nextflow task work directory
    set -o pipefail
    Rscript --vanilla "$spe_script" $sample_set
    test -s "02_build_spe/spe_n${sample_set}_assays.h5"
    test -s "02_build_spe/spe_n${sample_set}_se.rds"
    """
}
