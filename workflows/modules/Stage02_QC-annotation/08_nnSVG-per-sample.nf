process NNSVG_PER_SAMPLE {
    tag sample_id
    label 'mbv_r'
    container r_image
    // Finish the other samples; workflow.failOnIgnore keeps the run unsuccessful.
    errorStrategy 'ignore'

    input:
    path svg_script
    val r_image
    tuple val(sample_id), path(sample_spe), val(seed)

    output:
    // Because samples run in parallel and can finish in any order, saving as tuple lets downstream processes match results by sample ID without parsing filenames.
    tuple val(sample_id), path("${sample_id}_nnSVG-results.csv"), emit: svg_results

    script:
    """
    export TMPDIR=\$PWD
    Rscript --vanilla "$svg_script" "$sample_spe" "$sample_id" "${task.cpus}" "$seed"
    test -s "${sample_id}_nnSVG-results.csv"
    """
}
