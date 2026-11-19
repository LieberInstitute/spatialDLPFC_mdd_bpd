include { MAD_OUTLIER_DETECTION } from './modules/02_QC-annotation/01_MAD-outlier-detection.nf'


// pipeline parameters
params {
    // global workflow parameters
    container_runtime: String = ''
    r_docker_image: String = ''
    r_apptainer_image: String = ''
    // Accepted when reusing the stage 1 container params file.
    reference_gtf: String = ''
    // HDF5 spe
    n_samples: Integer = 119
    spe_dir: Path
    // module-specific parameters
    stage_subdir: String = '02_QC-annotation'
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

    spe_dir = file(params.spe_dir, checkIfExists: true, glob: false)
    MAD_OUTLIER_DETECTION(
        file("${rDir}/03_QC/01_find_edges-outlier-areas.R", checkIfExists: true),
        file("${rDir}/03_QC/raster_functions.R", checkIfExists: true),
        r_image, spe_dir, params.n_samples
    )

    publish:
    mad_outliers = MAD_OUTLIER_DETECTION.out.coldata
}

output {
    mad_outliers {
        path params.stage_subdir
    }
}
