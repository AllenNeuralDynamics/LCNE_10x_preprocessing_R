#02b make publication purpose fig 
source("./func/cluster_functions.R")

input_file  <- "../scratch/LC_clustered_harmony.rds"
combined <- readRDS(input_file)
lc_markers <- c("Dbh", "Slc6a2", "Th", "Slc18a2", "Tacr3")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# -------------------------
# Theme
# -------------------------
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title  = element_text(size = base_size + 1),
      axis.text   = element_text(size = base_size),
      plot.title  = element_text(size = base_size + 2, face = "bold"),
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size - 1),
      strip.text   = element_text(size = base_size, face = "bold"),
      panel.border = element_rect(fill = NA, linewidth = 0.35),
      plot.margin  = margin(6, 6, 6, 6)
    )
}

# -------------------------
# Save helper (PNG 500dpi + PDF)
# -------------------------
save_pub <- function(p, filename_base, w, h, outdir = "../results/figures",
                     source_data = NULL) {
  dir.create(outdir, showWarnings = FALSE)
  
  ggsave(file.path(outdir, paste0(filename_base, ".png")),
         plot = p, width = w, height = h,
         units = "in", dpi = 500, bg = "white")
  ggsave(file.path(outdir, paste0(filename_base, ".svg")),
         plot = p, width = w, height = h,
         units = "in",
         device = grDevices::svg,
         bg = "white")
  
  # source data for publication (graph values only)
  if (is.null(source_data) && !is.null(p$data)) source_data <- p$data
  if (!is.null(source_data)) {
    write.csv(source_data,
              file.path(outdir, paste0(filename_base, "_source_data.csv")),
              row.names = TRUE)   # rownames = cell barcodes / feature_cluster ids
  }
}





input_file = "../scratch/LC_cluster_interest.rds"  # 5151 cells to start with 
combined_lc <- readRDS(input_file)



## some QC Plot
combined_lc$percent.mt <- as.numeric(as.character(combined_lc$percent.mt))
hist(combined_lc$percent.mt,
     # breaks = 50,
     main = "Distribution of mitochondrial percentage",
     xlab = "percent.mt",
     col = "gray")


# =========================
# 1) Cluster UMAP
# =========================
p_umap_clusters <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) + NoLegend() + coord_fixed() + labs(title = "UMAP of Clusters", x = "UMAP 1", y = "UMAP 2")
p_umap_clusters

# save_pub(p_umap_clusters, "umap_clusters", w = 4.2, h = 4.0)
save_pub(p_umap_clusters, "fig_S8d", w = 4.2, h = 4.0)



# =========================
# 2) FeaturePlots (ordered layering)
# =========================
p_features <- FeaturePlot(
  combined,
  features = lc_markers[1:4],
  reduction = "umap",
  ncol = 2,
  pt.size = 1,
  min.cutoff = "q5",
  max.cutoff = "q95",
  cols = c("lightgrey", "blue"),
  order = FALSE,
  raster = TRUE
) &
  theme_pub(base_size = 8) &
  coord_fixed() &
  labs(x = "UMAP 1", y = "UMAP 2") &
  theme(legend.key.width  = unit(0.2, "cm"), legend.key.height = unit(0.3, "cm"))

p_features

# 
# save_pub(p_features, "umap_fedatureplots_lc_markers",
#          w = 3.9 * 3, h = 3.3 * ceiling(length(lc_markers) / 3),
#          source_data = cbind(Embeddings(combined, "umap")[, 1:2], FetchData(combined, vars = lc_markers[1:4])))

save_pub(p_features, "fig_S8c",
         w = 3.9 * 3, h = 3.3 * ceiling(length(lc_markers) / 3),
         source_data = cbind(Embeddings(combined, "umap")[, 1:2], FetchData(combined, vars = lc_markers[1:4])))




################



# =========================
# 3) DotPlot
# =========================
p_dot <- DotPlot(
  combined,
  scale = FALSE,   # this is new but i think this makes more sense. 
  features = lc_markers[1:4],
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  labs(title = "Marker expression by cluster",
       y = "Seurat clusters") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_dot

n_clusters <- length(levels(factor(combined$seurat_clusters)))

# save_pub(
#   p_dot,
#   "dotplot_lc_markers",
#   w = 6.5,
#   h = max(3.5, 0.18 * n_clusters + 2.0)
# )
save_pub(
  p_dot,
  "fig_S8b",
  w = 6.5,
  h = max(3.5, 0.18 * n_clusters + 2.0)
)



#######
# 4) highlight the slected cluster!
lc_candidate_clusters <- get_candidate_clusters(combined, lc_markers[1:3])
cells_cellsincluster <- WhichCells(combined, idents = lc_candidate_clusters)  # adjust if seurat_clusters isn't Idents
## 6607

p_umap_cellsincluster <- DimPlot(
  combined,
  reduction = "umap",
  cells.highlight = cells_cellsincluster,
  cols.highlight = "firebrick3",
  cols = "grey85",
  pt.size = 0.1,
  sizes.highlight = 0.8
) +
  NoLegend() +
  coord_fixed() +
  labs(title = "Cluster 7", x = "UMAP 1", y = "UMAP 2")

p_umap_cellsincluster


