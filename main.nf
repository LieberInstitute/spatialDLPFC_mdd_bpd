include { SPE_CONTAINER_CHECK } from './nextflow_modules/spe_test.nf'

// pipeline parameters
params {
    sample_set: Integer
    //sample_input: String
}


workflow {
    main:
    spe_script = file(
        "${projectDir}/code/02_build_spe/01_spe-from-GEO.R",
        checkIfExists: true
    )

    SPE_CONTAINER_CHECK(spe_script, params.sample_set)

    publish:
    first_output = SPE_CONTAINER_CHECK.out.geo_samples
    second_output = SPE_CONTAINER_CHECK.out.spe_files
}

output {
    first_output {
        path '02_build_spe'
    }
    second_output {
        path ''
    }
}
