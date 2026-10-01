# Visualizing Tapestri with Seurat
library(Seurat)
library(ggplot2)
library(tidyverse)

import_path <-  "data/Tapestri/HC02/"

# Import data
dna_vaf <- read_csv(paste0(import_path, "mosaic/dna_variants/layers/AF_MISSING.csv"))
dna_ngt <- read_csv(paste0(import_path, "mosaic/dna_variants/layers/NGT.csv"))
protein_norm <- read_csv(paste0(import_path, "mosaic/protein_read_counts/layers/normalized_counts.csv"))


dim(dna_vaf)
dim(dna_ngt)
dim(protein_norm)

# Put into dataframes to create Seurat Assays
vaf_table <- t(as.matrix(dna_vaf[,3:ncol(dna_vaf)]))
colnames(vaf_table) <- dna_vaf %>% dplyr::pull("AF_MISSING")

ngt_table <- t(as.matrix(dna_ngt[,3:ncol(dna_ngt)]))
colnames(ngt_table) <- dna_ngt %>% dplyr::pull("NGT")

protein_table <- t(as.matrix(protein_norm[,3:ncol(protein_norm)]))
colnames(protein_table) <- protein_norm %>% dplyr::pull("normalized_counts")


# Calculate drop out rates
ado <- numeric(length = nrow(vaf_table))
names(ado) <- rownames(vaf_table)
for(i in 1:nrow(vaf_table)){
  ado[i] <- sum(vaf_table[i,] == -50)/ncol(vaf_table)
}

# Re-encode missing values as NA instead of -50 or 0 (NGT assay for hierarchical clustering)
vaf_table_missing <- vaf_table
vaf_table[vaf_table == -50] <- NA
ngt_table[ngt_table == 3] <- NA

# Create assays
dna_assay <- CreateAssayObject(data = vaf_table)
dna_missing_assay <- CreateAssayObject(data = vaf_table_missing)
ngt_assay <- CreateAssayObject(data = ngt_table)


protein_assay <- CreateAssayObject(data = protein_table)

# Construct Seurat container
seu_tapestri <- CreateSeuratObject(dna_assay, assay = "DNA")
seu_tapestri[["DNA_missing"]] <- dna_missing_assay
seu_tapestri[["NGT"]] <- ngt_assay
seu_tapestri[["NGT"]] <- ngt_assay
seu_tapestri[["ADT"]] <- protein_assay

### Get feature-level metadata ###
GQ <- read_csv(paste0(import_path, "dna_variants/layers/GQ.csv"))
DP <- read_csv(paste0(import_path, "dna_variants/layers/DP.csv"))

GQ <- GQ[,c(3:ncol(GQ))]
DP <- DP[,c(3:ncol(DP))]

GQ_tidy <- GQ %>% gather(key = "gene", value = "GQ")

DP_tidy <- DP %>% gather(key = "gene", value = "DP")

# Get variants with good GQ scores & DP
GQ_median <- GQ_tidy %>% group_by(gene) %>% summarise(GQ = median(GQ))
DP_median <- DP_tidy %>% group_by(gene) %>% summarise(DP = median(DP))

variant_filter <- intersect(GQ_median %>% filter(GQ >= 40) %>% pull(gene),
                            DP_median %>% filter(DP >= 10) %>% pull(gene))


# Include only genes from WES whitelist
whitelist_HC02 <- read.csv("HC02_whitelist.csv", row.names = 1)
whitelist_mito <- read.csv("mito_whitelist.csv", row.names = 1)
# Remove gene names to match whitelist format
variant_filter_sub <- sub(".*?:", "", variant_filter)
whitelist <- variant_filter[variant_filter_sub%in%whitelist_HC02[,1]]
whitelist_full <- variant_filter[variant_filter_sub%in%whitelist_HC02[,1]|
                                   variant_filter_sub%in%whitelist_mito[,1]]

# Filter variants mutated in less than 1% of cells
variant_filter2 <- names(which(rowSums(seu_tapestri[["NGT"]]@data, na.rm = T)>0.01*ncol(seu_tapestri)))
whitelist <- whitelist[whitelist%in%variant_filter2]

# Filter variants with drop out rate larger than 5 %
ado_filter <- names(ado[ado<0.03])
whitelist <- whitelist[whitelist%in%ado_filter]
length(whitelist)

# Construct umap based on Hamming distance within NGT space
DefaultAssay(seu_tapestri) <- "NGT"
G <- FetchData(seu_tapestri, vars = whitelist, layer = "data")
cell_names <- rownames(G)

# Binarize for mutation presence and treat missing as 0 (wt)
Gb <- as.matrix(G)
Gb[Gb == 1L | Gb == 2L] <- 1L
Gb[Gb == 0L] <- 0L
storage.mode(Gb) <- "integer"
Gb[is.na(G)] <- 0

# Add this as assay to Seurat object
seu_tapestri[["NGT0"]] <- CreateAssayObject(data = t(Gb))


# 1) UMAP directly on binary with Hamming; also return neighbor lists
k <- 30
set.seed(1)
umap_outs <- uwot::umap(
  X            = Gb,
  metric       = "hamming",
  n_neighbors  = k,
  min_dist     = 0.2,
  n_components = 2,
  ret_nn       = TRUE,   # << get neighbor indices/distances
  verbose      = TRUE
)

emb <- umap_outs$embedding
rownames(emb) <- cell_names
colnames(emb) <- c("UMAP1","UMAP2")


# Add UMAP to Seurat object
DefaultAssay(seu_tapestri) <- "NGT0"

seu_tapestri[["umap"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "umap_",
  assay = "NGT0"
)

# Cluster on this umap
seu_tapestri <- FindNeighbors(seu_tapestri, reduction = "umap", dims = c(1,2))
seu_tapestri <- FindClusters(seu_tapestri, algorithm = 2, resolution = 0.0001)
DimPlot(seu_tapestri, label = T)

# Change to missing assay to identify informative variants as FindAllMarkers cannot handle NAs
DefaultAssay(seu_tapestri) <- "DNA_missing"
sig_muts <- FindAllMarkers(seu_tapestri, only.pos = T, features = whitelist)

library(viridisLite)
DoHeatmap(seu_tapestri, features = sig_muts$gene, slot = "data", disp.min = -50, disp.max = 100)+
  scale_fill_gradient2(low = "black",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Cluster 22 are doublets (shares variants at low VAF that are normally mutually exclusive)
seu_tapestri <- subset(seu_tapestri, idents = "22", invert = T)

# Cluster 7 is not clearly defined => fuse into polyclonal
seu_tapestri <- RenameIdents(seu_tapestri, "7" = "0")
seu_tapestri$fine <- Idents(seu_tapestri)


# Extract NGT0 values for hierarchical clustering
ngt_avg <- AverageExpression(seu_tapestri, features = sig_muts$gene, assays = "NGT0", return.seurat = T)
ngt_avg_bin <- ngt_avg[["NGT0"]]@layers$data
rownames(ngt_avg_bin) <- rownames(ngt_avg)
colnames(ngt_avg_bin) <- colnames(ngt_avg)
# Binarize
ngt_avg_bin[ngt_avg_bin >0.8] <- 1L
ngt_avg_bin[ngt_avg_bin <0.8] <- 0



# Hierchical clustering
d <- dist(t(ngt_avg_bin), method = "binary")
h_clust <- hclust(d, method = "ward.D2")
cutree(h_clust, k = 23)
length(unique(Idents(seu_tapestri)))


# Use hierarchical clustering to reorder clusters in Seurat object
# Manual for now because I haven't found a better way and 25 clusters are still
# manageable
plot(h_clust)
reordered_levels <- c("0","14","10","12","24",
                      "1", "17",
                      "6", "3", "4",  
                      "5", "20","21", "8", "2", "13", "19", 
                      "23","18", "9","15", 
                      "11","16"
)
ordered_clusters <- seu_tapestri$fine
ordered_clusters <- factor(ordered_clusters, levels = reordered_levels)

seu_tapestri$ordered <- ordered_clusters
Idents(seu_tapestri) <- seu_tapestri$ordered

# Identify clone-associated variants and plot
DefaultAssay(seu_tapestri) <- "DNA_missing"
sig_muts <- FindAllMarkers(seu_tapestri, only.pos = T, features = whitelist)

DoHeatmap(seu_tapestri, features = sig_muts$gene, slot = "data",
          disp.min = -50, disp.max = 100)+
  scale_fill_gradient2(low = "black",
                       mid = "#3B4D73",
                       high = "#D7ECEE",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Cluster 17 is identical to 1, 20 to 2
# 24 likely identical to 3 (enriched for RP1 FN)
seu_tapestri <- RenameIdents(seu_tapestri, "17" = "1", "20" = "2",
                             "24" = "3")

# P2RX1 mutation calling seems mostly mono-allelic => Disregard subclones defined by absence of mutation
seu_tapestri <- RenameIdents(seu_tapestri, "5" = "2", "21" = "8", "19" = "13")

# HAS3 mutation calling seems mostly mono-allelic => Disregard subclones defined by absence of mutation
seu_tapestri <- RenameIdents(seu_tapestri, "18" = "9", "23" = "15")


# Cluster 11 contains another subclone
temp  <- subset(seu_tapestri, idents = "11")
DefaultAssay(temp) <- "NGT0"
temp <- FindNeighbors(temp, reduction = "umap", dims = c(1,2))
temp <- FindClusters(temp, algorithm = 2, resolution = 0.02)
Idents(seu_tapestri, cells = WhichCells(temp, ident = "0")) <- "11a"
Idents(seu_tapestri, cells = WhichCells(temp, ident = "1")) <- "11b"
rm(temp)

# Cluster 12 contains another subclone
temp  <- subset(seu_tapestri, idents = "12")
DefaultAssay(temp) <- "NGT0"
temp <- FindNeighbors(temp, reduction = "umap", dims = c(1,2))
temp <- FindClusters(temp, algorithm = 2, resolution = 0.02)
Idents(seu_tapestri, cells = WhichCells(temp, ident = "0")) <- "12a"
Idents(seu_tapestri, cells = WhichCells(temp, ident = "1")) <- "12b"
rm(temp)


# Final ordering
reordered_levels <- c("0",
                      "14","10","12a", "12b",
                      "1", "6", "3", "4",
                      "8", "2", "13", 
                      "15", "9",
                      "11b","16","11a"
)
ordered_clusters <- Idents(seu_tapestri)
ordered_clusters <- factor(ordered_clusters, levels = reordered_levels)

seu_tapestri$ordered <- ordered_clusters
Idents(seu_tapestri) <- seu_tapestri$ordered


# Identify clone-associated variants and plot
DefaultAssay(seu_tapestri) <- "DNA_missing"
sig_muts <- FindAllMarkers(seu_tapestri, only.pos = T, features = whitelist)
# Plot reordered mutations
DoHeatmap(seu_tapestri, features = sig_muts$gene, slot = "data",draw.lines = F,
          disp.min = -50, disp.max = 100)+
  scale_fill_gradient2(low = "black",
                       mid = "#3B4D73",
                       high = "#D7ECEE",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Combine subclusters into larger families for plotting
seu_tapestri$fine <- Idents(seu_tapestri)
Idents(seu_tapestri) <- seu_tapestri$fine
# Keep cell ordering for plotting of subclones within larger families
cell.order <- names(seu_tapestri$fine[order(seu_tapestri$fine)])

seu_tapestri <- RenameIdents(seu_tapestri, 
                             "0" = "polyclonal",
                             "14" = "A", "10" = "A", "12a" = "A", "12b" = "A",
                             "1" = "A",  "6" = "A", "4" = "A", "3" = "A",
                             "8" = "B","2" = "B", "13" = "B",
                             "15" = "C","9" = "C",
                             "11b" = "D","16" = "D", "11a" = "D"
)
seu_tapestri$coarse <- Idents(seu_tapestri)

clone_colors = c("grey", "darkblue", "darkorange", "olivedrab", "tomato1")
# Plot reordered mutations of clonal families
DoHeatmap(seu_tapestri, features = sig_muts$gene, slot = "data",draw.lines = F,
          disp.min = -50, disp.max = 100, group.colors = clone_colors, cells = cell.order)+
  scale_fill_gradient2(low = "black",
                       mid = "#3B4D73",
                       high = "#D7ECEE",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")


# Plot Protein
DefaultAssay(seu_tapestri) <- "ADT"
seu_tapestri <- ScaleData(seu_tapestri)
DoHeatmap(seu_tapestri, features = c("NKG2C", "NKG2A", "NKp30", "CD57", "CD2",
                                     "KIR2DL2-L3-S2", "KIR2DL1-S1-S3-S5",
                                     "KIR3DL1", "CD161", "KLRG1"),
          disp.min = 0, disp.max = 2, lines.width = 5,
          slot = "data")+
  scale_fill_gradientn(colors = magma(9, alpha = 1, begin = 0, end = 1, direction = 1))

## Clone sizes
freq_df <- as.data.frame(table(seu_tapestri$coarse))


# Calculate VAFs
freq_df$VAF_totalNK <- freq_df$Freq / sum(freq_df$Freq)


ggplot(freq_df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity",  size = 1, color = "white",) +
  coord_polar("y") +
  theme_void() +  # Remove background and axes
  #geom_text(aes(label = Label), position = position_stack(vjust = 0.5)) +  # Add labels
  labs(title = "Pie Chart of Category Frequencies (Ordered by Factor Levels)") +
  scale_fill_manual(values = clone_colors)

# Export fully annotated object, stash informative mutations as variable features
DefaultAssay(seu_tapestri) <- "NGT0"
VariableFeatures(seu_tapestri) <- sig_muts$gene
saveRDS(seu_tapestri, file = "HC02_tapestri.rds")
