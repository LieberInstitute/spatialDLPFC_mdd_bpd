library(data.table)
library(ggplot2)
library(ggpubr)

# Read csv's and generate cell centroid coordinate 
indir <- "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/15_RNAscope_Validation/ObjectData"
files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)

dt_list <- lapply(files, function(f) {
  dt <- fread(f)
  fname <- basename(f)
  brain <- sub("^Br([0-9]+)_.*", "Br\\1", fname)
  dt[, file := fname]
  dt[, brain := brain]
  dt[, x_center := (XMin + XMax) / 2]
  
# Use x-axis for cortical depths since L1 is at xrel=0
  dr <- range(dt$x_center, na.rm = TRUE)
  dt[, x_rel := (x_center - dr[1]) / (dr[2] - dr[1])]
  dt
})

combined_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Melt to long format 
plot_dt <- melt(
  combined_dt,
  id.vars = c("brain", "file", "Object Id", "x_rel"),
  measure.vars = c("PCP4_620 Classification", "TAC1_690 Classification"),
  variable.name = "celltype",
  value.name = "positive"
)

# Clean marker names from HALO output
plot_dt[, celltype := fcase(
  celltype == "PCP4_620 Classification", "PCP4+",
  celltype == "TAC1_690 Classification", "TAC1+"
)]

# Generate mutually exclusive classification of PCP4+ vs. TAC1+
cell_class_dt <- plot_dt[, .(
  pos_sum = sum(positive),
  marker = as.character(celltype[positive == 1][1])
), by = .(brain, file, `Object Id`, x_rel)]

cell_class_dt[, category := "Negative"]
cell_class_dt[pos_sum == 1, category := marker]
cell_class_dt[pos_sum >= 2, category := "Multi-Positive"]

# Check cell counts per marker category
print(cell_class_dt[, .N, by = category])

# Colors each classification type
marker_cols <- c(
  "PCP4+"          = "orange",
  "TAC1+"          = "blue",
  "Multi-Positive" = "black",
  "Negative"       = "grey90"
)

# Keep only single-positive cells
single_pos_dt <- cell_class_dt[category %in% c("PCP4+", "TAC1+")]

p_density <- ggplot(single_pos_dt, 
                    aes(y = x_rel, fill = category, color = category)) +
  geom_density(alpha = 0.4, linewidth = 0.8, adjust = 1.5) +
  scale_y_reverse(expand = c(0, 0), limits = c(1, 0)) +
  scale_fill_manual(values = marker_cols) +
  scale_color_manual(values = marker_cols) +
  labs(
    y = "Relative Depth (0 = L1, 1 = WM)",
    x = "Abundance Density",
    fill = "Classification",
    color = "Classification"
  ) +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 2,
    plot.title = element_text(hjust = 0.5)
  )

print(p_density)

# Save plot to cluster and GitHub
outdir <- "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/plots/15_RNAscope_Validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(outdir, "PCP4_TAC1_density.png"),
  plot = p_density,
  width = 6,
  height = 10,
  dpi = 300
)

