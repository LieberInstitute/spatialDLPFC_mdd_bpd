include { CREATE_COLOR_PALETTES } from './modules/Stage02_QC-annotation/00_create-color-palettes.nf'
include { MAD_OUTLIER_DETECTION } from './modules/Stage02_QC-annotation/01_MAD-outlier-detection.nf'
include { REVIEW_EDGES_OUTLIER_AREAS } from './modules/Stage02_QC-annotation/02_review-edges-outlier-areas.nf'
include { SPOTSWEEPER } from './modules/Stage02_QC-annotation/03_spotsweeper.nf'
include { QC_PLOTS } from './modules/Stage02_QC-annotation/04_QC-plots.nf'


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
    spe_dir: Path = 'processed-data/Stage01_create-spe'
    // module-specific parameters
    stage_subdir: String = 'Stage02_QC-annotation'
    ignore_edge_samples: List = ['V13B23-302_A1', 'V13B23-302_C1']
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
    
    
    // run QC and annotation modules
    MAD_OUTLIER_DETECTION(
        file("${rDir}/03_QC/01_find_edges-outlier-areas.R", checkIfExists: true),
        file("${rDir}/03_QC/raster_functions.R", checkIfExists: true),
        r_image, spe_dir, params.n_samples
    )

    REVIEW_EDGES_OUTLIER_AREAS(
        file("${rDir}/03_QC/02_review_edges-outlier-areas.R", checkIfExists: true),
        r_image, spe_dir, params.n_samples, MAD_OUTLIER_DETECTION.out.coldata
    )

    SPOTSWEEPER(
        file("${rDir}/03_QC/03_spotsweeper.R", checkIfExists: true),
        r_image, spe_dir, params.n_samples, MAD_OUTLIER_DETECTION.out.coldata,
        params.ignore_edge_samples.join(',')
    )

    CREATE_COLOR_PALETTES(
        file("${launchDir}/plots/create_colorPalettes.R", checkIfExists: true),
        r_image
    )

    QC_PLOTS(
        file("${rDir}/03_QC/04_QC-plots.R", checkIfExists: true),
        r_image, SPOTSWEEPER.out.coldata, CREATE_COLOR_PALETTES.out.palettes
    )

    publish:
    color_palettes = CREATE_COLOR_PALETTES.out.palettes
    review_pdf = REVIEW_EDGES_OUTLIER_AREAS.out.review_pdf
    qc_coldata = SPOTSWEEPER.out.coldata
    qc_summary_pdf = QC_PLOTS.out.summary_pdf
}

output {
    color_palettes {
        path 'plots'
    }
    review_pdf {
        path { "plots/${params.stage_subdir}" }
    }
    qc_coldata {
        path { "processed-data/${params.stage_subdir}" }
    }
    qc_summary_pdf {
        path { "plots/${params.stage_subdir}" }
    }
}
