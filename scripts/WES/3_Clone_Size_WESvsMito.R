##### Comparative analysis of clone sizes estimated by WES / mtATAC

library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)


#### Extract clonotype sizes from scATAC data ####

#### Analysis of new donors HC09, HC12 ####

### HC12 ###
seu <- readRDS("data/ATAC/outs/H12_atac.rds")
seu$cluster <- Idents(seu)

# Find clonotypes
DefaultAssay(seu) <- "alleles"
seu <- FindClonotypes(seu, metric = "euclidean", k =4)
seu <- FindClusters(seu, resolution = 1.2, group.singletons = F, algorithm = 1)
seu$clonotype <- Idents(seu)


# Find significantly enriched mutations to identify meaningful clonotypes
sigmut  <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut <- sigmut  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)
filtered_clonotypes <- topmut$cluster[!topmut$cluster%in%c("singleton")]

# Clonotypes look well defined overall
DoHeatmap(seu, features = sigmut$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(seu, expression = clonotype %in% topmut$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()

# Identify clonotypes significantly associated with NKG2C+ adaptive population
# using Caleb's chi² approach
DimPlot(seu, label = T, group.by = "cluster")

NK_clonotypes<- data.frame(
  Clusters = seu$cluster,
  mito_cluster = paste0(seu$clonotype)
)

# Include only clonotypes which are clearly identified by mutational profile
NK_clonotypes <- NK_clonotypes %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes))

# Calculate chi² statistic
NK_clonotypes_obs <- get_chisquare_stats(NK_clonotypes); NK_clonotypes$what <- "Observed"
NK_clonotypes_perm <- get_chisquare_stats(NK_clonotypes, 1); NK_clonotypes$what <- "Permuted"

sig_clones <- NK_clonotypes_obs %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)


# Inspect significantly associated clonotypes to check which ones are in adaptive
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "27"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "30"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "37"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "34"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "40"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "39"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "43"))#yes


hc12_clones <- c("27", "30", "37", "34", "40", "39", "43")

# Filter informative clones
seu$clones_filtered <- factor(seu$clonotype,
                              levels = c(levels(seu$clonotype), "non_clones"))
seu$clones_filtered[!seu$clones_filtered%in%hc12_clones] <- "non_clones"

# subset NKG2C+ clusters as this is how we sorted for WES
FeaturePlot(seu, features = "NKG2C-Strep")
DimPlot(seu, group.by = "cluster", label = T)
seu_NKG2C <- subset(seu, cluster %in%c("3", "0"), invert = T)


# Quantify total clone frequency
clone_counts_HC12 <- as.data.frame(table(seu_NKG2C$clones_filtered))
colnames(clone_counts_HC12) <- c("clone_id", "count")

# Convert to frequency:
clone_counts_HC12 <- clone_counts_HC12 %>%
  mutate(
    frequency = count / sum(count),
    percentage = 100 * frequency
  )

# Filter for clones and add donor id
clone_counts_HC12 <- clone_counts_HC12 %>% filter(clone_id%in%hc12_clones) %>% 
  mutate(donor = "HC12")

### HC09 ###
seu <- readRDS("data/ATAC/outs/H9_atac.rds")
seu$cluster <- Idents(seu)



# Find clonotypes
DefaultAssay(seu) <- "alleles"
seu <- FindClonotypes(seu, metric = "euclidean", k = 4)
seu <- FindClusters(seu, resolution = 0.8, group.singletons = F, algorithm = 1)
seu$clonotype <- Idents(seu)


# Find significantly enriched mutations to identify meaningful clonotypes
sigmut  <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut <- sigmut  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)
filtered_clonotypes <- topmut$cluster[!topmut$cluster%in%c("singleton")]

# Clonotypes look well defined overall
DoHeatmap(seu, features = sigmut$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(seu, expression = clonotype %in% topmut$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()

# Identify clonotypes significantly associated with NKG2C+ adaptive population
# using Caleb's chi² approach
DimPlot(seu, label = T, group.by = "cluster")

NK_clonotypes<- data.frame(
  Clusters = seu$cluster,
  mito_cluster = paste0(seu$clonotype)
)

# Include only clonotypes which are clearly identified by mutational profile
NK_clonotypes <- NK_clonotypes %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes))

# Calculate chi² statistic
NK_clonotypes_obs <- get_chisquare_stats(NK_clonotypes); NK_clonotypes$what <- "Observed"
NK_clonotypes_perm <- get_chisquare_stats(NK_clonotypes, 1); NK_clonotypes$what <- "Permuted"

sig_clones <- NK_clonotypes_obs %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)


# Inspect significantly associated clonotypes to check which ones are in adaptive
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "38"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "24"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "19"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "28"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "32"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "36"))#yes
hc09_clones <- c("38", "19", "28", "32", "36")

# Filter informative clones
seu$clones_filtered <- factor(seu$clonotype,
                              levels = c(levels(seu$clonotype), "non_clones"))
seu$clones_filtered[!seu$clones_filtered%in%hc09_clones] <- "non_clones"

# subset NKG2C+ clusters as this is how we sorted for WES
FeaturePlot(seu, features = "NKG2C-Strep")
DimPlot(seu, group.by = "cluster", label = T)
seu_NKG2C <- subset(seu, cluster %in%c("1"), invert = T)


# Quantify total clone frequency
clone_counts_HC09 <- as.data.frame(table(seu_NKG2C$clones_filtered))
colnames(clone_counts_HC09) <- c("clone_id", "count")

# Convert to frequency:
clone_counts_HC09 <- clone_counts_HC09 %>%
  mutate(
    frequency = count / sum(count),
    percentage = 100 * frequency
  )


# Filter for clones and add donor id
clone_counts_HC09 <- clone_counts_HC09 %>% filter(clone_id%in%hc09_clones) %>% 
  mutate(donor = "HC09")


#### Import data from old donors published in Nat Imm and analyze with same approach ####

### CMVpos3 == HC06 ###
seu <- readRDS("data/ATAC/inputs/CMVpos3_full_alleles.rds")

seu$cluster <- Idents(seu)


# Find clonotypes
DefaultAssay(seu) <- "alleles"
seu <- FindClonotypes(seu, metric = "euclidean", k = 4)
seu <- FindClusters(seu, resolution = 1.4, group.singletons = F, algorithm = 1)
seu$clonotype <- Idents(seu)


# Find significantly enriched mutations to identify meaningful clonotypes
sigmut  <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut <- sigmut  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)
filtered_clonotypes <- topmut$cluster[!topmut$cluster%in%c("singleton")]

# Clonotypes look well defined overall
DoHeatmap(seu, features = c(sigmut$gene, "14422T>C"), slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(seu, expression = clonotype %in% topmut$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()

# Identify clonotypes significantly associated with NKG2C+ adaptive population
# using Caleb's chi² approach
DimPlot(seu, label = T, group.by = "cluster")

NK_clonotypes <- data.frame(
  Clusters = seu$cluster,
  mito_cluster = seu$clonotype
)

# Include only clonotypes which are clearly identified by mutational profile
NK_clonotypes <- NK_clonotypes %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes))

# Calculate chi² statistic
NK_clonotypes_obs <- get_chisquare_stats(NK_clonotypes); NK_clonotypes$what <- "Observed"
NK_clonotypes_perm <- get_chisquare_stats(NK_clonotypes, 1); NK_clonotypes$what <- "Permuted"

#
sig_clones <- NK_clonotypes_obs %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)

# Inspect significantly associated clonotypes to check which ones are in adaptive
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "6"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "11"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "7"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "12"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "8"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "19"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "18"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "15"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "34"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "75"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "25"))#yes


hc06_clones <- c("6", "11",  "7", "8", "18", "15", "34", "25")

# Filter informative clones
seu$clones_filtered <- factor(seu$clonotype,
                              levels = c(levels(seu$clonotype), "non_clones"))
seu$clones_filtered[!seu$clones_filtered%in%hc06_clones] <- "non_clones"

# subset NKG2C+ clusters as this is how we sorted for WES
DimPlot(seu, group.by = "cluster", label = T)
seu_NKG2C <- subset(seu, cluster %in%c("Adaptive1", "Adaptive2", "Adaptive3"))


# Quantify total clone frequency
clone_counts_HC06 <- as.data.frame(table(seu_NKG2C$clones_filtered))
colnames(clone_counts_HC06) <- c("clone_id", "count")

# Convert to frequency:
clone_counts_HC06 <- clone_counts_HC06 %>%
  mutate(
    frequency = count / sum(count),
    percentage = 100 * frequency
  )

# Filter for clones and add donor id
clone_counts_HC06 <- clone_counts_HC06 %>% filter(clone_id%in%hc06_clones) %>% 
  mutate(donor = "HC06")


### CMVpos02 == P2 ###
seu <- readRDS("data/ATAC/inputs/CMVpos2_full_alleles.rds")

seu$cluster <- Idents(seu)


# Find clonotypes
DefaultAssay(seu) <- "alleles"
seu <- FindClonotypes(seu, metric = "euclidean")
seu <- FindClusters(seu, resolution = 1.3, group.singletons = F, algorithm = 1)
seu$clonotype <- Idents(seu)


# Find significantly enriched mutations to identify meaningful clonotypes
sigmut  <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut <- sigmut  %>% group_by(cluster) %>% top_n(n = 1, avg_log2FC)
filtered_clonotypes <- topmut$cluster[!topmut$cluster%in%c("singleton")]

# Clonotypes look well defined overall
DoHeatmap(seu, features = sigmut$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(seu, expression = clonotype %in% topmut$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()



# Identify clonotypes significantly associated with NKG2C+ adaptive population
# using Caleb's chi² approach
DimPlot(seu, label = T, group.by = "cluster")

NK_clonotypes <- data.frame(
  Clusters = seu$cluster,
  mito_cluster = seu$clonotype
)

# Include only clonotypes which are clearly identified by mutational profile
NK_clonotypes <- NK_clonotypes %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes))

# Calculate chi² statistic
NK_clonotypes_obs <- get_chisquare_stats(NK_clonotypes); NK_clonotypes$what <- "Observed"
NK_clonotypes_perm <- get_chisquare_stats(NK_clonotypes, 1); NK_clonotypes$what <- "Permuted"

sig_clones <- NK_clonotypes_obs %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)

# Inspect significantly associated clonotypes to check which ones are in adaptive
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "3"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "0"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "15"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "12"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "20"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "18"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "22"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "6"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "23"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "4"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "8"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "10"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "21"))#no


PC2_clones <- c("3", "15", "12","20",  "18", "22", "23")

# Filter informative clones
seu$clones_filtered <- factor(seu$clonotype,
                              levels = c(levels(seu$clonotype), "non_clones"))
seu$clones_filtered[!seu$clones_filtered%in%PC2_clones] <- "non_clones"

# subset NKG2C+ clusters as this is how we sorted for WES
DimPlot(seu, group.by = "cluster", label = T)
seu_NKG2C <- subset(seu, cluster %in%c("Adaptive1", "Adaptive2", "Adaptive3"))


# Quantify total clone frequency
clone_counts_PC02 <- as.data.frame(table(seu_NKG2C$clones_filtered))
colnames(clone_counts_PC02) <- c("clone_id", "count")

# Convert to frequency:
clone_counts_PC02 <- clone_counts_PC02 %>%
  mutate(
    frequency = count / sum(count),
    percentage = 100 * frequency
  )

# Filter for clones and add donor id
clone_counts_PC02 <- clone_counts_PC02 %>% filter(clone_id%in%PC2_clones) %>% 
  mutate(donor = "PC02")




### HC05 == CMVpos4 ###
seu <- readRDS("data/ATAC/inputs/CMVpos4_full_alleles")

seu$cluster <- Idents(seu)


# Find clonotypes
DefaultAssay(seu) <- "alleles"
seu <- FindClonotypes(seu, metric = "euclidean", k = 4)
seu <- FindClusters(seu, resolution = 1.2, group.singletons = F, algorithm = 1)
seu$clonotype <- Idents(seu)


# Find significantly enriched mutations to identify meaningful clonotypes
sigmut  <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut <- sigmut  %>% group_by(cluster) %>% top_n(n = 2, avg_log2FC)
filtered_clonotypes <- topmut$cluster[!topmut$cluster%in%c("singleton")]

# Clonotypes look well defined overall
DoHeatmap(seu, features = sigmut$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(seu, expression = clonotype %in% topmut$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()



# Identify clonotypes significantly associated with NKG2C+ adaptive population
# using Caleb's chi² approach
DimPlot(seu, label = T, group.by = "cluster")

NK_clonotypes <- data.frame(
  Clusters = seu$cluster,
  mito_cluster = seu$clonotype
)

# Include only clonotypes which are clearly identified by mutational profile
NK_clonotypes <- NK_clonotypes %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes))

# Calculate chi² statistic
NK_clonotypes_obs <- get_chisquare_stats(NK_clonotypes); NK_clonotypes$what <- "Observed"
NK_clonotypes_perm <- get_chisquare_stats(NK_clonotypes, 1); NK_clonotypes$what <- "Permuted"

sig_clones <- NK_clonotypes_obs %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)

# Inspect significantly associated clonotypes to check which ones are in adaptive
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "66"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "20"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "7"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "8"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "13"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "47"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "11"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "53"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "14"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "69"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "21"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "44"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "23"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "39"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "40"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "24"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "9"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "46"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "37"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "32"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "36"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "62"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "56"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "63"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "22"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "25"))#no
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "26"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "31"))#no

HC05_clones <- c("66", "20", "7", "8", "47", "53", "14", "69", "21", "40", "9",
                 "32", "36", "62", "26")

# Filter informative clones
seu$clones_filtered <- factor(seu$clonotype,
                              levels = c(levels(seu$clonotype), "non_clones"))
seu$clones_filtered[!seu$clones_filtered%in%HC05_clones] <- "non_clones"

# subset NKG2C+ clusters as this is how we sorted for WES
DimPlot(seu, group.by = "cluster", label = T)
seu_NKG2C <- subset(seu, cluster %in%c("Adaptive1", "Adaptive2", "Adaptive3"))


# Quantify total clone frequency
clone_counts_HC05 <- as.data.frame(table(seu_NKG2C$clones_filtered))
colnames(clone_counts_HC05) <- c("clone_id", "count")

# Convert to frequency:
clone_counts_HC05 <- clone_counts_HC05 %>%
  mutate(
    frequency = count / sum(count),
    percentage = 100 * frequency
  )

# Filter for clones and add donor id
clone_counts_HC05 <- clone_counts_HC05 %>% filter(clone_id%in%HC05_clones) %>% 
  mutate(donor = "HC05")

### Old donor: CMVpos1 == HC02 ###
seu <- readRDS("data/ATAC/inputs/CMVpos1_full.rds")

# Subset to mtASAP5
seu <- subset(seu, experiment == "mtASAP5")
DimPlot(seu)
seu$cluster <- Idents(seu)


# Find clonotypes
DefaultAssay(seu) <- "alleles"
seu <- FindClonotypes(seu, metric = "euclidean", k = 4)
seu <- FindClusters(seu, resolution = 1.2, group.singletons = F, algorithm = 1)
seu$clonotype <- Idents(seu)


# Find significantly enriched mutations to identify meaningful clonotypes
sigmut  <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.1, min.pct = 0.6)
topmut <- sigmut  %>% group_by(cluster) %>% top_n(n = 2, avg_log2FC)
filtered_clonotypes <- topmut$cluster[!topmut$cluster%in%c("singleton")]

# Clonotypes look well defined overall
DoHeatmap(seu, features = sigmut$gene, slot = "data", disp.max = 0.1, group.by = "clonotype",
          cells = WhichCells(seu, expression = clonotype %in% topmut$cluster & clonotype!= "singleton"))+
  scale_fill_gradientn(colours = c("white", "red"))+NoLegend()



# Identify clonotypes significantly associated with NKG2C+ adaptive population
# using Caleb's chi² approach
DimPlot(seu, label = T, group.by = "cluster")

NK_clonotypes <- data.frame(
  Clusters = seu$cluster,
  mito_cluster = seu$clonotype
)

# Include only clonotypes which are clearly identified by mutational profile
NK_clonotypes <- NK_clonotypes %>% dplyr::filter(mito_cluster%in%paste0(filtered_clonotypes))

# Calculate chi² statistic
NK_clonotypes_obs <- get_chisquare_stats(NK_clonotypes); NK_clonotypes$what <- "Observed"
NK_clonotypes_perm <- get_chisquare_stats(NK_clonotypes, 1); NK_clonotypes$what <- "Permuted"

sig_clones <- NK_clonotypes_obs %>% dplyr::filter(FDR< 0.05 ) %>% dplyr::select(clone_id)

# Inspect significantly associated clonotypes to check which ones are in adaptive
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "5"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "2"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "4"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "0"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "6"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "3"))#yes
DimPlot(seu, cells.highlight = WhichCells(seu, idents = "8"))#yes
HC02_clones <- sig_clones$clone_id

# Filter informative clones
seu$clones_filtered <- factor(seu$clonotype,
                              levels = c(levels(seu$clonotype), "non_clones"))
seu$clones_filtered[!seu$clones_filtered%in%HC02_clones] <- "non_clones"



# subset NKG2C+ clusters as this is how we sorted for WES
FeaturePlot(seu, features = "Anti-PE")
DimPlot(seu, group.by = "cluster", label = T)
seu_NKG2C <- subset(seu, cluster %in%c("Adaptive1", "Adaptive2", "Adaptive3", "Adaptive5"))


# Quantify total clone frequency
clone_counts_HC02 <- as.data.frame(table(seu_NKG2C$clones_filtered))
colnames(clone_counts_HC02) <- c("clone_id", "count")

# Convert to frequency:
clone_counts_HC02 <- clone_counts_HC02 %>%
  mutate(
    frequency = count / sum(count),
    percentage = 100 * frequency
  )
# Filter for clones and add donor id
clone_counts_HC02 <- clone_counts_HC02 %>% filter(clone_id%in%HC02_clones) %>% 
  mutate(donor = "HC02")

#### Combine all donors ####

clone_counts_all <- rbind(clone_counts_HC02, clone_counts_HC12,
                          clone_counts_HC09,clone_counts_HC05,
                          clone_counts_HC06,clone_counts_PC02)

clone_freq_summarized <- clone_counts_all %>% group_by(donor) %>% 
  summarize(cum_freq = sum(frequency))

write.csv(clone_counts_all,"data/ATAC/outs/mito_clone_size.csv")


### Combine with scifer results ###
clone_df <- read.csv("data/WES/outs/scifer_df.csv")
clone_df <- merge(clone_df, clone_freq_summarized)

plot_df <- clone_df %>%
  mutate(donor = factor(donor, levels = unique(donor))) %>%
  pivot_longer(
    cols = c(f, cum_freq),
    names_to = "measure",
    values_to = "frequency"
  ) %>%
  mutate(
    measure = factor(measure, levels = c("f", "cum_freq")),
    lower = if_else(measure == "f", f_low, NA_real_),
    upper = if_else(measure == "f", f_upp, NA_real_)
  )

ggplot(plot_df, aes(x = donor, y = frequency, fill = measure)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.2,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(f = "#4477AA", cum_freq = "#EE9944"),
    labels = c(f = "f", cum_freq = "Cumulative frequency")
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Donor", y = "Frequency", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")