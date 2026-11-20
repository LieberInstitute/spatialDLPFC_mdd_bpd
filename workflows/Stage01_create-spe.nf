include { DOWNLOAD_GEO_ARCHIVE } from './modules/Stage01_create-spe/00_GEO-download-archive.nf'
include { FETCH_GEO } from './modules/Stage01_create-spe/01_GEO-download.nf'
include { SPE_GEO_BUILD } from './modules/Stage01_create-spe/02_spe_build-from-GEO.nf'
include { GENE_REF } from './modules/Stage01_create-spe/03_spe_update-rowData.nf'

// pipeline parameters
params {
    // global workflow parameters
    container_runtime: String = ''
    r_docker_image: String = ''
    r_apptainer_image: String = ''
    // full sample set (over-written if test profile used)
    n_samples: Integer = 119
    // module-specific parameters
    stage_subdir: String = 'Stage01_create-spe'
    timeout: Integer = 300
    geo_archive_resolve_ip: String = ''
    reference_gtf: String = ''
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

    // GEOquery performance is unstable so try curl tar fetch
    DOWNLOAD_GEO_ARCHIVE(r_image, params.timeout, params.geo_archive_resolve_ip ?: 'default')

    // run spe creation workflow modules
    FETCH_GEO("${rDir}/02_build_spe/01_GSE333974.R", r_image, 
        params.n_samples, params.timeout, DOWNLOAD_GEO_ARCHIVE.out.raw_archive)
    
    SPE_GEO_BUILD("${rDir}/02_build_spe/02_spe-from-GEO.R", r_image, 
        params.n_samples,
        FETCH_GEO.out.geo_bundle, FETCH_GEO.out.geo_samples)

    // check for GTF file
    gtf_file = params.reference_gtf ? file(params.reference_gtf, glob: false) : null
    gtf_input = gtf_file != null && gtf_file.isFile() ? gtf_file : []

    GENE_REF("${rDir}/02_build_spe/03_update-rowData.R", r_image, 
        params.n_samples, params.timeout,
        SPE_GEO_BUILD.out.spe_files, params.reference_gtf, gtf_input)

    publish:
    geo_path = FETCH_GEO.out.geo_bundle
    geo_samples = FETCH_GEO.out.geo_samples
    spe_files = GENE_REF.out.spe_files
}

output {
    geo_path {
        path { "processed-data/${params.stage_subdir}" }
    }
    geo_samples {
        path { "processed-data/${params.stage_subdir}" }
    }
    spe_files {
        path { "processed-data/${params.stage_subdir}" }
    }
}
