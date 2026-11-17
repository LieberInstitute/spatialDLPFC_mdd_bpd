include { FETCH_GEO } from './modules/01_create-spe/01_GEO-download.nf'
include { SPE_GEO_BUILD } from './modules/01_create-spe/02_spe_build-from-GEO.nf'
include { GENE_REF } from './modules/01_create-spe/03_spe_update-rowData.nf'

// pipeline parameters
params {
    // global workflow parameters
    container_runtime: String = ''
    r_docker_image: String = ''
    r_apptainer_image: String = ''
    // full sample set (over-written if test profile used)
    n_samples: Integer = 119
    // module-specific parameters
    stage_subdir: String = '01_create-spe'
    timeout: Integer = 300
}


workflow {
    main:
    // check container runtime declared with profile and set image
    if (params.container_runtime == 'docker') {
        if (!params.r_docker_image) {
            error 'Set r_docker_image in your Nextflow params .yaml file.'
        }
        r_image = params.r_docker_image
    } else if (params.container_runtime == 'apptainer') {
        if (!params.r_apptainer_image.startsWith('/')) {
            error 'Set r_apptainer_image to an absolute SIF path in your Nextflow params .yaml file.'
        }
        r_image = "file://${params.r_apptainer_image}"
    } else {
        error 'Select the docker or apptainer profile.'
    }

    // set working directory for R scripts
    rDir = file("${launchDir}/code/")

    // run spe creation workflow modules
    FETCH_GEO("${rDir}/02_build_spe/01_GSE333974.R", r_image, params.n_samples, params.timeout)
    SPE_GEO_BUILD("${rDir}/02_build_spe/02_spe-from-GEO.R", r_image, params.n_samples,
        FETCH_GEO.out.geo_bundle, FETCH_GEO.out.geo_samples)
    GENE_REF("${rDir}/02_build_spe/03_update-rowData.R", r_image, params.n_samples, params.timeout,
        SPE_GEO_BUILD.out.spe_files)

    publish:
    geo_path = FETCH_GEO.out.geo_bundle
    geo_samples = FETCH_GEO.out.geo_samples
    spe_files = GENE_REF.out.spe_files
}

output {
    geo_path {
        path params.stage_subdir
    }
    geo_samples {
        path params.stage_subdir
    }
    spe_files {
        path params.stage_subdir
    }
}
