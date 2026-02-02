## modifying DeconvoBuddies code because it keeps not giving me all the genes for my full size sn sce object but I already have the means calculated in the dotplot version so I'll just import that instead
get_mean_ratio_mod <- function(sce,
    cellType_col,
    #assay_name = "logcounts",
    min_prop.detected = .2,
    gene_ensembl = NULL,
    gene_name = NULL) {
    # RCMD fix
    cellType.target <- NULL
    cellType <- NULL
    ratio <- NULL
    rank_ratio <- NULL
    anno_ratio <- NULL

    #specify logcounts.mean
    assay_name = "logcounts.mean"
    
    ## check inputs are valid
    stopifnot(cellType_col %in% colnames(colData(sce)))
    stopifnot(assay_name %in% names(SummarizedExperiment::assays(sce)))

    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

#    ct_table <- table(sce[[cellType_col]])
#
#    if (any(ct_table < 10)) warning("One or more cell types has < 10 cells, this may result in unstable marker genes results. Check details of get_mean_ratio() for more info")

    sce_assay <- SummarizedExperiment::assays(sce)[[assay_name]]

    ## Get mean expression for each gene for each cellType
#    cell_means <- map(cell_types, ~ as.data.frame(MatrixGenerics::rowMeans(sce_assay[, sce[[cellType_col]] == .x])))
    cell_means <- purrr::map(cell_types, ~ as.data.frame(sce_assay[, sce[[cellType_col]] == .x]))

    cell_means <- do.call("rbind", cell_means)
    colnames(cell_means) <- "mean"
    ## Define columns
    cell_means$cellType <- rep(cell_types, each = nrow(sce))
    cell_means$gene <- rep(rownames(sce), length(cell_types))
    # print(head(cell_means))

    prop_assay = assays(sce)[["logcounts.prop.detected"]]
    cell_props <- purrr::map(cell_types, ~ as.data.frame(prop_assay[, sce[[cellType_col]] == .x]))
    cell_props <- do.call("rbind", cell_props)
    colnames(cell_props) <- "prop.detected"
    cell_means = cbind(cell_means, cell_props)

    ## Filter and calculate ratio for each celltype
    ratio_tables <- purrr::map(cell_types, ~ .get_ratio_table_mod(
        .x,
        sce,
        sce_assay,
        cellType_col,
        cell_means,
	min_prop.detected
    ))

    ratio_tables <- do.call("rbind", ratio_tables) |>
        mutate(anno_ratio = paste0(cellType.target, "/", cellType, ": ", base::round(ratio, 3))) |>
        dplyr::rename(
            cellType.2nd = cellType,
            mean.2nd = mean,
            MeanRatio = ratio,
            MeanRatio.rank = rank_ratio,
            MeanRatio.anno = anno_ratio
        )

    ## Add gene ensemble and gene_name if specified
    if (!is.null(gene_ensembl)) {
        if (gene_ensembl %in% colnames(SummarizedExperiment::rowData(sce))) {
            ratio_tables$gene_ensembl <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ][[gene_ensembl]]
        } else {
            warning("'", gene_ensembl, "' not in col rowData, gene_ensembl not included in output")
        }
    }

    if (!is.null(gene_name)) {
        if (gene_name %in% colnames(SummarizedExperiment::rowData(sce))) {
            ratio_tables$gene_name <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ][[gene_name]]
        } else {
            warning("'", gene_name, "' not in col rowData, gene_name not included in output")
        }
    }

    return(ratio_tables)
}


.get_ratio_table_mod <- function(x, sce, sce_assay, cellType_col, cell_means, min_prop.detected) {
    # RCMD Fix
    mean.target <- NULL
    gene <- NULL
    ratio <- NULL
    cellType.target <- NULL
    cellType <- NULL
    
    
#    # filter target median != 0
#    median_index <- MatrixGenerics::rowMedians(sce_assay[, sce[[cellType_col]] == x]) != 0
#    # message("Median == 0: ", sum(!median_index))
#    # filter for target means
    target_mean <- cell_means[cell_means$cellType == x, ]
#    target_mean <- target_mean[median_index, ]

    # filter by prop detected instead of median
    prop_index = target_mean[["prop.detected"]]>= min_prop.detected
    target_mean = target_mean[prop_index,]

    colnames(target_mean) <- c("mean.target", "cellType.target", "gene", "prop.detected")

    nontarget_mean <- cell_means[cell_means$cellType != x, 1:3]

    ratio_table <- dplyr::left_join(target_mean, nontarget_mean, by = "gene") |>
        mutate(ratio = mean.target / mean) |>
        dplyr::group_by(gene) |>
        arrange(ratio) |>
        dplyr::slice(1) |>
        dplyr::select(gene, cellType.target, mean.target, prop.detected, cellType, mean, ratio) |>
        arrange(-ratio) |>
        dplyr::ungroup() |>
        mutate(rank_ratio = dplyr::row_number())

    return(ratio_table)
}
