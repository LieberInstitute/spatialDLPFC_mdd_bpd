include { FETCH_GEO } from './nextflow_modules/01_GEO-download.nf'
include { SPE_GEO_BUILD } from './nextflow_modules/02_spe_build-from-GEO.nf'
include { GENE_REF } from './nextflow_modules/03_spe_update-rowData.nf'

// pipeline parameters
params {
    n_samples: Integer
    //sample_input: String
}


workflow {
    main:
    rDir = file("${projectDir}/code/")
    FETCH_GEO("${rDir}/02_build_spe/01_GSE333974.R", params.n_samples)
    SPE_GEO_BUILD("${rDir}/02_build_spe/02_spe-from-GEO.R", params.n_samples,
        FETCH_GEO.out.geo_bundle, FETCH_GEO.out.geo_samples)
    GENE_REF("${rDir}/02_build_spe/03_update-rowData.R", params.n_samples,
        SPE_GEO_BUILD.out.spe_files)

    publish:
    geo_path = FETCH_GEO.out.geo_bundle
    geo_samples = FETCH_GEO.out.geo_samples
    spe_files = GENE_REF.out.spe_files
}

output {
    geo_path {
        path '01_create-spe'
    }
    geo_samples {
        path '01_create-spe'
    }
    spe_files {
        path '01_create-spe'
    }
}
