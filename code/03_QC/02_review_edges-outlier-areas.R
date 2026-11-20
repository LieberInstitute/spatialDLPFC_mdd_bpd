args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Usage: Rscript 02_review_edges-outlier-areas.R <spe_dir> <n_samples> <qc_df>")
}
spe_dir <- args[[1L]]
n_samples <- as.integer(args[[2L]])
spe_prefix <- paste0("spe_n", n_samples, "_")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(HDF5Array)
    library(ggspavis)
    library(gridExtra)
    library(ggplot2)
    library(ggrastr)
})


# LOAD SPE AND QC ANNOTATIONS #

spe <- loadHDF5SummarizedExperiment(dir=spe_dir, prefix=spe_prefix)
stopifnot(length(unique(spe$sample_id)) == n_samples)
spe$lg10.umi = log10(spe$sum_umi+1)

qc.areas = read.csv(args[[3L]], row.names=1)
stopifnot(identical(colnames(spe), rownames(qc.areas)))


# TRANSFER QC ANNOTATIONS TO SPE AND FORMAT FOR PLOTTING #

spe$edge_outlier = qc.areas$edge_outlier
spe$outlier_area = qc.areas$outlier_area
spe$check_area = as.character(
    factor(paste(spe$edge_outlier, spe$outlier_area),
    levels = c("TRUE TRUE", "TRUE FALSE", "FALSE TRUE", "FALSE FALSE"),
    labels = c("edge artifact", "edge artifact", "low UMI area", "none"))
)
spe$check_area = ifelse(spe$sum_umi<100 & spe$check_area=="none", "<100 UMI", spe$check_area)
spe$check_area <- factor(spe$check_area, levels=c("edge artifact", "low UMI area", "<100 UMI", "none"))


# PLOT ONLY SAMPLES WITH EDGE OR CONTIGUOUS OUTLIER AREA; SAVE PDF (EMPTY IF NONE) #

sids = unique(colData(spe)[spe$check_area %in% c("edge artifact", "low UMI area"), "sample_id"])

slist <- lapply(sids, function(x) spe[,spe$sample_id==x & spe$in_tissue])
names(slist) <- sids
if (length(sids) == 0L) {
    pdf("check_edges-outlier-areas.pdf", height=3, width=7)
    grid::grid.text("No samples have edge artifacts or low UMI areas.")
    dev.off()
} else {
    # Use a fileless device while constructing grobs to avoid Rplots.pdf.
    pdf(NULL, height=3, width=7)
    plist.sort <- unlist(lapply(slist, function(x) {
        sample_images <- imgData(x)
        use_image <- nrow(sample_images) > 0L && all(vapply(sample_images$data, function(image) {
            is(image, "LoadedSpatialImage") ||
                (is(image, "StoredSpatialImage") && file.exists(imgSource(image)))
        }, logical(1)))
        if (!use_image) warning("Images unavailable for ", unique(x$sample_id), "; plotting spots only.")
        plots <- list(
            suppressWarnings(plotVisium(x, annotate = "lg10.umi", highlight = NULL,
                image = use_image, spots = TRUE, point_size = .5, pal = "viridis")+
                scale_fill_viridis_c(direction = -1)),
            plotVisium(x, annotate = "check_area", highlight = NULL, 
                image = use_image, spots = TRUE, point_size = .5,
                pal=c("edge artifact"="navy", "low UMI area"="#00a000", "<100 UMI"="red", "none"="grey"))
        )

        for (p in seq_along(plots)) {
            for (layer in seq_along(plots[[p]]$layers)) {
                if (inherits(plots[[p]]$layers[[layer]]$geom, "GeomPoint")) {
                    plots[[p]]$layers[[layer]]$aes_params$stroke <- 0
                }
            }
        }
        lapply(plots, function(x) rasterize(x, layers = "Point", dpi = 250))
    }), recursive=FALSE)
    review_pages <- marrangeGrob(grobs=plist.sort, nrow=1, ncol=2, top=quote(sids[[g]]))
    dev.off() # end fileless device
    ggsave(filename="check_edges-outlier-areas.pdf", plot=review_pages,
        height=3, width=7)
}
