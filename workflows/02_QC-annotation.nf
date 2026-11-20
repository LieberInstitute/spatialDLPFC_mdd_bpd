include { MAD_OUTLIER_DETECTION } from './modules/02_QC-annotation/01_MAD-outlier-detection.nf'
include { REVIEW_EDGES_OUTLIER_AREAS } from './modules/02_QC-annotation/02_review-edges-outlier-areas.nf'


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
    //plots_dir: String = ''
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

    // Mirror output runs beneath plots; an explicit plots_dir overrides this.
    //output_root = workflow.outputDir
    //default_output_root = file("${launchDir}/nextflow_outputs")
    //plot_suffix = output_root.startsWith(default_output_root)
    //    ? default_output_root.relativize(output_root).toString()
    //    : output_root.fileName.toString()
    //plots_dir = params.plots_dir
    //    ? file(params.plots_dir)
    //    : file("${launchDir}/plots").resolve(plot_suffix)
    
    
    // run QC and annotation modules
    MAD_OUTLIER_DETECTION(
        file("${rDir}/03_QC/01_find_edges-outlier-areas.R", checkIfExists: true),
        file("${rDir}/03_QC/raster_functions.R", checkIfExists: true),
        r_image, spe_dir, params.n_samples
    )

    REVIEW_EDGES_OUTLIER_AREAS(
        file("${rDir}/03_QC/02_review_edges-outlier-areas.R", checkIfExists: true),
        r_image, spe_dir, params.n_samples, MAD_OUTLIER_DETECTION.out.coldata,
        //plots_dir.resolve(params.stage_subdir)
    )

    publish:
    mad_outliers = MAD_OUTLIER_DETECTION.out.coldata
    review_pdf = REVIEW_EDGES_OUTLIER_AREAS.out.review_pdf
}

output {
    mad_outliers {
        path { "processed-data/${params.stage_subdir}" }
    }
    review_pdf {
        path { "plots/${params.stage_subdir}" }
    }
}
