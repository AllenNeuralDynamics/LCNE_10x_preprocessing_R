#02b make publication purpose fig 
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
save_pub <- function(p, filename_base, w, h, outdir = "../figures") {
  dir.create(outdir, showWarnings = FALSE)
  
  ggsave(file.path(outdir, paste0(filename_base, ".png")),
         plot = p, width = w, height = h,
         units = "in", dpi = 500, bg = "white")
  ggsave(file.path(outdir, paste0(filename_base, ".svg")),
         plot = p, width = w, height = h,
         units = "in",
         device = grDevices::svg,
         bg = "white")
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
) + NoLegend()
p_umap_clusters


save_pub(p_umap_clusters, "umap_clusters", w = 4.2, h = 4.0)



# =========================
# 2) FeaturePlots (ordered layering)
# =========================

p_features <- FeaturePlot(
  combined,
  features = lc_markers,
  reduction = "umap",
  ncol = 3,
  pt.size = 1,
  min.cutoff = "q1",
  max.cutoff = "q99",
  cols = c("lightgrey", "firebrick3"),
  order = FALSE,
  raster = TRUE,
) &
  theme_pub(base_size = 8)
p_features


save_pub(
  p_features,
  "umap_fedatureplots_lc_markers",
  w = 3.9 * 3,
  h = 3.3 * ceiling(length(lc_markers) / 3)
)

####### testing ########
# 
# p1 <- FeaturePlot(combined, lc_markers, pt.size=0.5, order=TRUE,  raster=FALSE)
# p2 <- FeaturePlot(combined, lc_markers, pt.size=0.5, order=TRUE,  raster=TRUE)
# p3 <- FeaturePlot(combined, lc_markers, pt.size=0.5, order=FALSE, raster=FALSE)
# p4 <- FeaturePlot(combined, lc_markers, pt.size=0.5, order=FALSE, raster=TRUE)

################



# =========================
# 3) DotPlot
# =========================
p_dot <- DotPlot(
  combined,
  scale = FALSE,   # this is new but i think this makes more sense. 
  features = lc_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  labs(title = "Marker expression by cluster",
       y = "Seurat clusters") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

n_clusters <- length(levels(factor(combined$seurat_clusters)))

save_pub(
  p_dot,
  "dotplot_lc_markers",
  w = 6.5,
  h = max(3.5, 0.18 * n_clusters + 2.0)
)



