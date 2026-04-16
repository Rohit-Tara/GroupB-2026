library(Seurat)
library(Matrix)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(FNN)
library(igraph)
library(patchwork)

# load saved annotated objects
load("annotated_clusters_final_12apr26.RData")

# rebuild mapping objects
# not saved in the .RData
data <- read.csv("/home/obai/Documents/BI_Term_2/SRP/Datasets/expression_matrix_full.csv")

# Create lookup: gene_id -> gene_name
gene_map <- setNames(data$gene_name, data$gene_id)
# Create reverse lookup: gene_name -> gene_id
name_to_id <- setNames(data$gene_id, data$gene_name)

# Derive replicate metadata from geo_cell_id:
# Extract patterns "reg1" to "reg4"
# Rename to shorter format "r1"–"r4"
seu$replicate <- str_extract(seu$geo_cell_id, "reg[1-4]")
seu$replicate <- ifelse(
  is.na(seu$replicate),
  NA,                   # keep NA if no replicate pattern found
  sub("reg", "r", seu$replicate)  # convert "regX" -> "rX
)

# Initialize fetal_type metadata column
seu$fetal_type <- NA_character_
# Label cells belonging to fetal subclusters (e.g. "_fetal_10wpc_c1")
seu$fetal_type[grepl("_fetal_[0-9]+wpc_c[0-9]+", seu$geo_cell_id)] <- "fetal_subcluster"
# Label cells belonging to whole fetal samples (e.g. "_fetal_10wpc")
seu$fetal_type[grepl("_fetal_[0-9]+wpc$", seu$geo_cell_id)] <- "fetal_sample"


# Attach cell-level metadata from Seurat object to t-SNE dataframe
# Matching is done by aligning tsne_df$sample with Seurat cell names (colnames)

# Convert cluster to character to avoid factor-related issues
tsne_df$cluster <- as.character(tsne_df$cluster)

# Add replicate, day, and fetal_type annotations via index matching
tsne_df$replicate <- seu$replicate[match(tsne_df$sample, colnames(seu))]
tsne_df$day <- seu$day[match(tsne_df$sample, colnames(seu))]
tsne_df$fetal_type <- seu$fetal_type[match(tsne_df$sample, colnames(seu))]
tsne_df$cluster <- as.character(tsne_df$cluster)
tsne_df$replicate <- seu$replicate[match(tsne_df$sample, colnames(seu))]
tsne_df$day <- seu$day[match(tsne_df$sample, colnames(seu))]
tsne_df$fetal_type <- seu$fetal_type[match(tsne_df$sample, colnames(seu))]

# cluster annotation
cluster_annotation <- data.frame(
  cluster = c("0", "1", "2", "3", "4", "5", "6"),
  cluster_code = c("C1", "C2", "C3", "C4", "C5", "C6", "C7"),
  cell_type = c(
    "Late prog",
    "Dorsal NPC",
    "Dorsal neuron",
    "Mesenchymal",
    "Immature neuron",
    "Neuron subtype",
    "RSPO+ hem"
  ),
  stringsAsFactors = FALSE
)

# Build plotting dataframe with selected t-SNE coordinates and metadata
tsne_plot_stage <- tsne_df[, c(
  "sample", "tSNE1", "tSNE2",
  "cluster", "replicate", "day", "fetal_type"
)]

# Add cluster annotations by joining on cluster ID
# all.x = TRUE keeps all t-SNE rows even if annotation is missing
# sort = FALSE preserves original row order
tsne_plot_stage <- merge(
  tsne_plot_stage,
  cluster_annotation,
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

# define stage / sample symbols
tsne_plot_stage$shape_group <- NA_character_

# Assign shape groups for specific developmental days
tsne_plot_stage$shape_group[tsne_plot_stage$day == "33d"] <- "33d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "35d"] <- "35d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "37d"] <- "37d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "41d"] <- "41d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "65d"] <- "65d"

# Assign replicate-specific labels for selected days
tsne_plot_stage$shape_group[tsne_plot_stage$day == "53d" & tsne_plot_stage$replicate == "r1"] <- "r1_53d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "53d" & tsne_plot_stage$replicate == "r2"] <- "r2_53d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "58d" & tsne_plot_stage$replicate == "r3"] <- "r3_58d"
tsne_plot_stage$shape_group[tsne_plot_stage$day == "58d" & tsne_plot_stage$replicate == "r4"] <- "r4_58d"

# Assign fetal subcluster labels for remaining unlabeled entries
tsne_plot_stage$shape_group[
  is.na(tsne_plot_stage$shape_group) &
    tsne_plot_stage$fetal_type == "fetal_sample"
] <- "fetal_sample"

tsne_plot_stage$shape_group[
  is.na(tsne_plot_stage$shape_group) &
    tsne_plot_stage$fetal_type == "fetal_subcluster"
] <- "fetal_subcluster"

# check
table(tsne_plot_stage$shape_group, useNA = "ifany")

# map plotting symbols
pch_map <- c(
  "33d" = 19,
  "35d" = 17,
  "37d" = 25,
  "41d" = 18,
  "65d" = 15,
  "r1_53d" = 0,
  "r2_53d" = 5,
  "r3_58d" = 1,
  "r4_58d" = 2,
  "fetal_sample" = 4,
  "fetal_subcluster" = 8
)

tsne_plot_stage$pch <- unname(pch_map[tsne_plot_stage$shape_group])

# label positions
label_pos_stage <- aggregate(
  cbind(tSNE1, tSNE2) ~ cluster + cluster_code,
  data = tsne_plot_stage,
  FUN = median
)

# figure 3.D

# Adjust plotting parameters:
# xpd = TRUE allows drawing outside plot region (for legends)
# mar expands right margin to fit legends
par(xpd = TRUE, mar = c(5, 4, 4, 10))

# Scatter plot of t-SNE coordinates:
# Point shape (pch) encodes stage/replicate (shape_group)
# Color encodes cluster identity
plot(
  tsne_plot_stage$tSNE1,
  tsne_plot_stage$tSNE2,
  pch = tsne_plot_stage$pch,
  col = as.factor(tsne_plot_stage$cluster_code),
  xlab = "tSNE1",
  ylab = "tSNE2",
  main = "PCA and unbiased clustering using t-SNE"
)

# Add cluster labels at precomputed centroid positions
text(
  label_pos_stage$tSNE1,
  label_pos_stage$tSNE2,
  labels = label_pos_stage$cluster_code,
  cex = 1.1,
  font = 2,
  pos = 3
)

# Legend for cluster identities:
# Colors correspond to cluster codes
# Positioned outside plot area on the right
legend(
  x = par("usr")[2] + 5,
  y = par("usr")[4],
  legend = paste(cluster_annotation$cluster_code, cluster_annotation$cell_type),
  col = seq_along(cluster_annotation$cluster_code),
  pch = 16,
  cex = 0.7,
  bty = "n",
  title = "Clusters"
)

# Legend for experimental groups:
# Point shapes (pch) represent developmental stage and replicate conditions
legend(
  x = par("usr")[2] + 5,
  y = par("usr")[4] - 0.3 * diff(par("usr")[3:4]),
  legend = c(
    "33d", "35d", "37d", "41d", "65d",
    "r1 53d", "r2 53d", "r3 58d", "r4 58d",
    "fetal sample", "fetal subcluster"
  ),
  pch = c(19, 17, 25, 18, 15, 0, 5, 1, 2, 4, 8),
  cex = 0.7,
  bty = "n",
  title = "experiments"
)

# figure 3.E

# Define genes to display in Figure 3E
genes_fig3E <- c(
  "FOXG1", "OTX2", "RSPO2", "DCN",
  "ASPM", "LIN28A", "MYT1L", "NEUROD6"
)

# Convert gene symbols to Ensembl IDs for lookup in the Seurat object
gene_ids_fig3E <- name_to_id[genes_fig3E]

# Drop genes that could not be mapped to an Ensembl ID
gene_ids_fig3E <- gene_ids_fig3E[!is.na(gene_ids_fig3E)]

# Print resolved gene symbol ↔ Ensembl ID mapping for verification
print(data.frame(
  gene_name = names(gene_ids_fig3E),
  gene_id = unname(gene_ids_fig3E)
))

# List available dimensional reductions in the Seurat object
Reductions(seu)

# Compute t-SNE embedding only if it is not already present
if (!"tsne" %in% names(seu@reductions)) {
  seu <- RunTSNE(seu, dims = 1:20)
}

# Generate one feature plot per gene using the t-SNE embedding
plots <- FeaturePlot(
  seu,
  features = unname(gene_ids_fig3E),
  reduction = "tsne",
  cols = c("lightgrey", "red"),
  combine = FALSE
)

# Replace Ensembl ID plot titles with human-readable gene symbols
for (i in seq_along(plots)) {
  plots[[i]] <- plots[[i]] + ggplot2::ggtitle(names(gene_ids_fig3E)[i])
}

# Arrange the eight gene expression plots in a 2 x 4 layout
p_fig3E <- patchwork::wrap_plots(plots, ncol = 4)

# Display the combined figure
print(p_fig3E)

# only if needed!
# Troubleshooting plotting issues (e.g. blank/overlapping plots)
# Close the current graphics device (run multiple times in case several are open)
dev.off()
# Reset all graphics devices and redraw plot
while (!is.null(dev.list())) dev.off()
print(p_fig3E)
# end of troubleshooting 

# figure 3.F
# Initialize list to store per-gene plots
plots_fig3F <- list()

# Loop over selected genes (Ensembl IDs mapped to gene names)
for (i in seq_along(gene_ids_fig3F)) {
  gid <- unname(gene_ids_fig3F[i])     # Ensembl gene ID
  gname <- names(gene_ids_fig3F)[i]    # Gene symbol
  
  # Extract expression values and region group metadata from Seurat object
  df <- FetchData(seu, vars = c(gid, "region_group"))
  colnames(df) <- c("expr", "region_group")
  
  # # Remove cells without region annotation
  df <- df[!is.na(df$region_group), ]
  
  # Set consistent ordering of region groups for plotting
  df$region_group <- factor(df$region_group, levels = c("r1", "r2", "r3", "r4", "fetal"))
  
  # Create violin plot with overlaid jittered points
  p <- ggplot(df, aes(x = region_group, y = expr)) +
    geom_violin(
      fill = "grey70",    # violin fill color
      color = "black",    # outline color
      trim = FALSE,       # show full distribution tails
      width = 0.8,
      adjust = 1.5,       # smoothing factor
      scale = "width"     # normalize widths across groups
    ) +
    geom_jitter(
      width = 0.12,       # horizontal jitter
      height = 0,
      size = 0.3,
      alpha = 0.35,
      color = "black"
    ) +
    ggtitle(gname) +      # use gene symbol as title
    xlab(NULL) +
    ylab("LOG2 FPKM") +
    
    # Fix y-axis range for consistent comparison across genes
    coord_cartesian(ylim = c(-1, 11), expand = FALSE) +
    # Custom y-axis ticks and labels (sparse labeling for clarity)
    scale_y_continuous(
      breaks = c(0.5, 2.5, 5, 7.5, 9.5),
      labels = c("", "0", "", "10", "")
    ) +
    # Apply clean theme and adjust styling
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  # Store plot in list
  plots_fig3F[[i]] <- p
}

# Arrange all gene plots into a grid (3 columns)
p_fig3F <- patchwork::wrap_plots(plots_fig3F, ncol = 3)

# Enforce square aspect ratio for all subplots
p_fig3F <- p_fig3F & theme(
  aspect.ratio = 1
)

# Display final figure
print(p_fig3F)

# Extra checks for violin plot interpretation
# Select a subset of marker genes to inspect in more detail
genes_check <- c("FOXG1", "NEUROD6", "OTX2")
# Map gene symbols to Ensembl IDs and remove any unmapped genes
gene_ids_check <- name_to_id[genes_check]
gene_ids_check <- gene_ids_check[!is.na(gene_ids_check)]
# Compute average expression of selected genes by region group
avg_region <- AggregateExpression(
  seu,
  features = unname(gene_ids_check),
  group.by = "region_group"
)$RNA
# Convert to data frame and add readable gene annotations
avg_region <- as.data.frame(avg_region)
avg_region$gene_id <- rownames(avg_region)
avg_region$gene_name <- gene_map[avg_region$gene_id]
# Inspect average expression table by region
avg_region
# Recreate gene mapping for presence/absence analysis
genes_check <- c("FOXG1", "NEUROD6", "OTX2")
gene_ids_check <- name_to_id[genes_check]
gene_ids_check <- gene_ids_check[!is.na(gene_ids_check)]
# Initialize list to store percent-positive results by region group
presence_list <- list()
# For each gene, calculate the fraction of cells with detectable expression (> 0) in each region group
for (i in seq_along(gene_ids_check)) {
  gid <- unname(gene_ids_check[i])
  gname <- names(gene_ids_check)[i]
  # Extract per-cell expression and region annotation
  df <- FetchData(seu, vars = c(gid, "region_group"))
  colnames(df) <- c("expr", "region_group")
  df <- df[!is.na(df$region_group), ]
  # Compute proportion of expressing cells by region group
  out <- aggregate(expr ~ region_group, data = df, FUN = function(x) mean(x > 0))
  out$gene_name <- gname
  out$percent_positive <- out$expr * 100
  out$expr <- NULL
  presence_list[[i]] <- out
}
# Combine per-gene results into one table
presence_region <- do.call(rbind, presence_list)
# Inspect percent-positive table by region
presence_region

# check expression by cluster
genes_check <- c("FOXG1", "NEUROD6", "OTX2")
gene_ids_check <- name_to_id[genes_check]
gene_ids_check <- gene_ids_check[!is.na(gene_ids_check)]

avg_cluster <- AggregateExpression(
  seu,
  features = unname(gene_ids_check),
  group.by = "ident"
)$RNA
# Convert to data frame and add readable gene annotations
avg_cluster <- as.data.frame(avg_cluster)
avg_cluster$gene_id <- rownames(avg_cluster)
avg_cluster$gene_name <- gene_map[avg_cluster$gene_id]
# Inspect average expression table by cluster
avg_cluster
# Initialize list to store percent-positive results by cluster
presence_cluster_list <- list()
# For each gene, calculate the fraction of cells with detectable expression (> 0) in each cluster
for (i in seq_along(gene_ids_check)) {
  gid <- unname(gene_ids_check[i])
  gname <- names(gene_ids_check)[i]
  # Extract per-cell expression and cluster identity
  df <- FetchData(seu, vars = c(gid, "ident"))
  colnames(df) <- c("expr", "cluster")
  # Compute proportion of expressing cells by cluster
  out <- aggregate(expr ~ cluster, data = df, FUN = function(x) mean(x > 0))
  out$gene_name <- gname
  out$percent_positive <- out$expr * 100
  out$expr <- NULL
  
  presence_cluster_list[[i]] <- out
}
# Combine per-gene results into one table
presence_cluster <- do.call(rbind, presence_cluster_list)
# Inspect percent-positive table by cluster
presence_cluster

# back-up stage
# save primary
save(
  seu,
  tsne_df,
  cluster_annotation,
  tsne_plot_stage,
  avg_region,
  presence_region,
  avg_cluster,
  presence_cluster,
  gene_map,
  name_to_id,
  gene_ids_check,
  file = "annotated_clusters_complete_13apr26.RData"
)

# save secondary
write.csv(tsne_plot_stage, "tsne_plot_stage_13apr26.csv", row.names = FALSE)
write.csv(avg_region, "avg_region_13apr26.csv", row.names = FALSE)
write.csv(presence_region, "presence_region_13apr26.csv", row.names = FALSE)
write.csv(avg_cluster, "avg_cluster_13apr26.csv", row.names = FALSE)
write.csv(presence_cluster, "presence_cluster_13apr26.csv", row.names = FALSE)
write.csv(seu@meta.data, "seurat_metadata_13apr26.csv", row.names = TRUE)





