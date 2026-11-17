include { FETCH_GEO } from './modules/01_create-spe/01_GEO-download.nf'
include { SPE_GEO_BUILD } from './modules/01_create-spe/02_spe_build-from-GEO.nf'
include { GENE_REF } from './modules/01_create-spe/03_spe_update-rowData.nf'

// pipeline parameters
params {
    n_samples: Integer
    stage_subdir: String = '01_create-spe'
    timeout: Integer = 300
}


workflow {
    main:
    rDir = file("${launchDir}/code/")
    FETCH_GEO("${rDir}/02_build_spe/01_GSE333974.R", params.n_samples, params.timeout)
    SPE_GEO_BUILD("${rDir}/02_build_spe/02_spe-from-GEO.R", params.n_samples,
        FETCH_GEO.out.geo_bundle, FETCH_GEO.out.geo_samples)
    GENE_REF("${rDir}/02_build_spe/03_update-rowData.R", params.n_samples, params.timeout,
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
