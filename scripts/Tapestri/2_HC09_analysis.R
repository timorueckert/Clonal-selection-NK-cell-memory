# Visualizing Tapestri with Seurat
library(Seurat)
library(ggplot2)
library(tidyverse)

import_path <- "data/Tapestri/HC09/"

### Import data
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
seu_tapestri[["ADT"]] <- protein_assay

### Get feature-level metadata ###
GQ <- read_csv(paste0(import_path, "mosaic/dna_variants/layers/GQ.csv"))
DP <- read_csv(paste0(import_path, "mosaic/dna_variants/layers/DP.csv"))

GQ <- GQ[,c(3:ncol(GQ))]
DP <- DP[,c(3:ncol(DP))]

GQ_tidy <- GQ %>% gather(key = "gene", value = "GQ")

DP_tidy <- DP %>% gather(key = "gene", value = "DP")

# Get variants with good GQ scores & DP
GQ_median <- GQ_tidy %>% group_by(gene) %>% summarise(GQ = median(GQ))
DP_median <- DP_tidy %>% group_by(gene) %>% summarise(DP = median(DP))

variant_filter <- intersect(GQ_median %>% filter(GQ >= 40) %>% pull(gene),
                            DP_median %>% filter(DP >= 10) %>% pull(gene))

### Variant filtering
# Include only genes from WES whitelist
whitelist_HC09 <- read.csv("HC09_whitelist.csv", row.names = 1)
whitelist_mito <- read.csv("mito_whitelist.csv", row.names = 1)
# Remove gene names to match whitelist format
variant_filter_sub <- sub(".*?:", "", variant_filter)
whitelist <- variant_filter[variant_filter_sub%in%whitelist_HC09[,1]]
whitelist_full <- variant_filter[variant_filter_sub%in%whitelist_HC09[,1]|
                                   variant_filter_sub%in%whitelist_mito[,1]]

# Filter variants mutated in less than 1% of cells
variant_filter2 <- names(which(rowSums(seu_tapestri[["NGT"]]@data, na.rm = T)>0.01*ncol(seu_tapestri)))
whitelist <- whitelist[whitelist%in%variant_filter2]

# Filter variants with drop out rate larger than 5 %
ado_filter <- names(ado[ado<0.05])
whitelist <- whitelist[whitelist%in%ado_filter]
length(whitelist)

# MN1 and ARVCF variants pass filters but are extremely noisy => exclude
whitelist <- whitelist[!whitelist%in%c("MN1:chr22:27797702:T/C", "ARVCF:chr22:19981403:G/C")]


### Clonotype Clustering
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

# ARHGEF33 suffers from drop-outs and distorts clustering => exclude here
Gb <- Gb[,!colnames(Gb)=="ARHGEF33:chr2:38959950:C/T"]



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
seu_tapestri <- FindClusters(seu_tapestri, algorithm = 2, resolution = 0.001)
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



# Cluster 8, 25, 26 are not defined by any variants => fuse into polyclonal
seu_tapestri <- RenameIdents(seu_tapestri, "8" = "0", "25" = "0", "26" = "0")

# 2 and 4; 7 and 14; 28 and 16 are identical
seu_tapestri <- RenameIdents(seu_tapestri, "2" = "4", "14" = "7", "28" = "16")

#Fuse 29/7/12 (differences are too sparse to be clearly resolved)
seu_tapestri <- RenameIdents(seu_tapestri, "29" = "7", "12" = "7")

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
cutree(h_clust, k = 22)
length(unique(Idents(seu_tapestri)))
plot(h_clust)

#
# Use hierarchical clustering to reorder clusters in Seurat object
# Manual for now because I haven't found a better way and number of clusters is still
# manageable
reordered_levels <- c("0",
                      "22","17",
                      "4", "11", "19", "5",
                      "3", "1", "6","7", 
                      "13", "9",
                      "27","15", "23","21",  "18",
                      "10",
                      "24", "16", "20" 
)
ordered_clusters <- seu_tapestri$fine
ordered_clusters <- factor(ordered_clusters, levels = reordered_levels)

seu_tapestri$ordered <- ordered_clusters
Idents(seu_tapestri) <- seu_tapestri$ordered


# Identify clone-associated variants and plot
DefaultAssay(seu_tapestri) <- "DNA_missing"
sig_muts <- FindAllMarkers(seu_tapestri, only.pos = T, features = whitelist)
DoHeatmap(seu_tapestri, features = sig_muts$gene, slot = "data", disp.min = -50, disp.max = 100)+
  scale_fill_gradient2(low = "black",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")




# Combine subclusters into larger families
seu_tapestri$fine <- Idents(seu_tapestri)
# Keep cell ordering for plotting of subclones within larger families
cell.order <- names(seu_tapestri$fine[order(seu_tapestri$fine)])

seu_tapestri <- RenameIdents(seu_tapestri, 
                             "0" = "non-clones",
                             "22" = "A", "17" = "A", "4" ="A", "11" = "A",
                             "19" = "A",  "5" = "A", "3" = "A", "1" = "A",
                             "6" = "A", "7" = "A",  "13" = "A",
                             "9" = "A","27" = "A", "15" = "A", "23" = "A",
                             "21" = "A", "18" = "A", "10" = "A",
                             "16" = "B", "20" = "B", "24" = "B"
)
seu_tapestri$coarse <- Idents(seu_tapestri)


# Plot reordered mutations of clonal families
clone_colors <- c("darkgrey","darkblue", "darkorange")
DoHeatmap(seu_tapestri, features = sig_muts$gene, slot = "data",
          disp.min = -50, disp.max = 100, cells = cell.order, group.colors = clone_colors)+
  scale_fill_gradient2(low = "black",
                       mid = "#3B4D73",
                       high = "#D7ECEE",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Plot protein markers
DefaultAssay(seu_tapestri) <- "ADT"
seu_tapestri <- ScaleData(seu_tapestri)
DoHeatmap(seu_tapestri, features = c("NKG2C", "NKG2A", "NKp30", "CD57", "CD2",
                                     "KIR2DL2-L3-S2", "KIR2DL1-S1-S3-S5",
                                     "KIR3DL1", "CD161", "KLRG1"),
          disp.min = 0, disp.max = 2, lines.width = 5,
          slot = "data")+
  scale_fill_gradientn(colors = magma(9, alpha = 1, begin = 0, end = 1, direction = 1))


## Clone sizes
freq_df <- as.data.frame(table(Idents(seu_tapestri)))


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
saveRDS(seu_tapestri, file = "HC09_tapestri.rds")
