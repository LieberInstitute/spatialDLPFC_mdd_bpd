library(data.table)
library(ggplot2)
library(ggpubr)

# ── 1. Read & normalize ────────────────────────────────────────────────────────
indir <- "/Users/madeline.abramson/Desktop/HALO MBv/ObjectData"
files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)

dt_list <- lapply(files, function(f) {
  dt <- fread(f)
  fname <- basename(f)
  brain <- sub("^Br([0-9]+)_.*", "Br\\1", fname)
  dt[, file := fname]
  dt[, brain := brain]
  dt[, x_center := (XMin + XMax) / 2]
  
  # All brains use X as cortical depth, all confirmed L1 at xrel=0
  dr <- range(dt$x_center, na.rm = TRUE)
  dt[, x_rel := (x_center - dr[1]) / (dr[2] - dr[1])]
  dt
})

combined_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# ── 2. Melt to long format ────────────────────────────────────────────────────
plot_dt <- melt(
  combined_dt,
  id.vars = c("brain", "file", "Object Id", "x_rel"),
  measure.vars = c("PCP4_620 Classification", "TAC1_690 Classification"),
  variable.name = "celltype",
  value.name = "positive"
)

# Clean marker names
plot_dt[, celltype := fcase(
  celltype == "PCP4_620 Classification", "PCP4+",
  celltype == "TAC1_690 Classification", "TAC1+"
)]

# ── 3. Mutually exclusive classification (matches published method) ────────────
cell_class_dt <- plot_dt[, .(
  pos_sum = sum(positive),
  marker = as.character(celltype[positive == 1][1])
), by = .(brain, file, `Object Id`, x_rel)]

cell_class_dt[, category := "Negative"]
cell_class_dt[pos_sum == 1, category := marker]
cell_class_dt[pos_sum >= 2, category := "Multi-Positive"]

# Check cell counts per category
print(cell_class_dt[, .N, by = category])

# ── 4. Colors ─────────────────────────────────────────────────────────────────
marker_cols <- c(
  "PCP4+"          = "red",
  "TAC1+"          = "green",
  "Multi-Positive" = "black",
  "Negative"       = "grey90"
)

# ── 5. Keep only single-positive cells
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

# ── 6. Generate histogram version due to smoothing effect of geom_density()
p_hist <- ggplot(single_pos_dt,
                 aes(x = x_rel, fill = category, color = category,
                     y = after_stat(density))) +
  geom_histogram(alpha = 0.3, linewidth = 0.3,
                 position = "identity",
                 bins = 30) +
  geom_density(alpha = 0, linewidth = 0.8, adjust = 1.5) +  # smooth line on top
  coord_flip() +
  scale_x_reverse() +
  scale_fill_manual(values = marker_cols) +
  scale_color_manual(values = marker_cols) +
  labs(
    y = "Abundance Density",
    x = "Relative Depth (0 = L1, 1 = WM)",
    fill = "Classification",
    color = "Classification"
  ) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 2)

print(p_hist)