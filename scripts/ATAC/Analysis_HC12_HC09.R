library(tidyverse)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(reticulate)
library(RColorBrewer)


# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
clusterpal <- brewer.pal(n = 9, "Blues")
clusterpal <- clusterpal[c(5,7,9)]
NKG2Cpal <- c("#40BAD5","#120136")
adaptivepal <- c("#880E4F","#7E191B", "#C0424E", "#E0115F","#BF0A30",  "#D21F3C",  "#FA8072","#CA3433",  "#CD5C5C" )
clusterpal <- c(clusterpal, adaptivepal)


# adding gene annotations
annotations <- readRDS('~/static/hg38_ucsc_annotations.rds')

import_path <- "input/ATAC"

# loading peaks
peaks <- read.table(paste0(import_path, 'cellranger_outs/peaks.bed'), col.names = c('chr', 'start', 'end')) %>%
  makeGRangesFromDataFrame()



# loading cells
metadata <- read.csv(paste0(import_path, 'cellranger_outs/singlecell.csv'))
metadata <- metadata %>% dplyr::filter(total > 500, barcode != "NO_BARCODE")
# loading fragments
fragments <- CreateFragmentObject(paste0(import_path, 'cellranger_outs/fragments.tsv.gz'), cells = metadata$barcode)
# creating feature counts
counts <- FeatureMatrix(fragments, peaks, cells = metadata$barcode)
# creating seurat object
alldata.atac <- CreateChromatinAssay(counts, fragments = fragments, features = peaks) %>%
  CreateSeuratObject(assay = 'ATAC', project = 'mtASAP6')
alldata.atac

# removing doublets
genotype_ids <- read.delim(paste0(import_path, 'vireo/donor_ids.tsv'), header = T)
genotype_ids <- genotype_ids[genotype_ids$cell %in% colnames(alldata.atac),]
genotype_ids <- genotype_ids[order(match(genotype_ids$cell, colnames(alldata.atac))), ]
name_index <- which(genotype_ids$donor_id == 'unassigned')
genotype_ids$donor_id[name_index] <- 'Unclear'
donor_id <- genotype_ids$donor_id
names(donor_id) <- genotype_ids$cell
alldata.atac <- AddMetaData(alldata.atac, metadata = donor_id, col.name = 'genotype')

amulet_doublets <- readLines(paste0(import_path, 'AMULET/MultipletBarcodes_01.txt'))
alldata.atac$amu_doublet <- FALSE
alldata.atac$amu_doublet[amulet_doublets] <- TRUE
table(alldata.atac$amu_doublet)

rm(genotype_ids, amulet_doublets)


#adding antibody tags and hashtags
folder <- paste0(import_path,'ADT/')
adt_counts <- t(readMM(file = paste0(folder, 'output.mtx')))
features <- read.delim(paste0(folder, 'output.genes.txt'), header = FALSE)
barcodes <- read.delim(paste0(folder, 'output.barcodes.txt'), header = FALSE)
colnames(adt_counts) <- barcodes$V1 %>% paste0('-1')
rownames(adt_counts) <- features$V1
adt_counts <- as(adt_counts, 'CsparseMatrix')
table(colnames(alldata.atac)%in%colnames(adt_counts))
adt_counts <- adt_counts[, colnames(adt_counts) %in% colnames(alldata.atac)]
alldata.atac <- subset(alldata.atac, cells = colnames(adt_counts))
alldata.atac[['ADT']] <- CreateAssayObject(adt_counts[1:15, ])
alldata.atac[['HTO']] <- CreateAssayObject(adt_counts[18:19, ])

rm(folder, adt_counts, features, barcodes)

# normalizing hashtags
alldata.atac <- NormalizeData(alldata.atac, assay = 'HTO', normalization.method = 'CLR', verbose = F)


# QC metrics

total_fragments <- CountFragments(paste0(import_path, 'cellranger_outs/fragments.tsv.gz'))
total_fragments <- total_fragments[total_fragments$CB %in% colnames(alldata.atac),]
rownames(total_fragments) <- total_fragments$CB
alldata.atac$nCount_fragments <- total_fragments[colnames(alldata.atac), "frequency_count"]

# Add annotations
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"


DefaultAssay(alldata.atac) <- 'ATAC'

Annotation(alldata.atac) <- annotations
alldata.atac <- NucleosomeSignal(alldata.atac)
alldata.atac <- TSSEnrichment(alldata.atac, process_n = 5000)
alldata.atac$blacklist_fraction <- FractionCountsInRegion(alldata.atac, assay = 'ATAC', regions = blacklist_hg38_unified)


# plotting qc metrics and setting appropriate filters

# Check first doublet and donor identification
FeatureScatter(alldata.atac, feature1 = "Hashtag9", feature2 = "Hashtag10", group.by = "genotype")
Idents(alldata.atac) <- alldata.atac$genotype

VlnPlot(
  object = alldata.atac,
  features = c('nCount_ATAC', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 2
)
VlnPlot(
  object = alldata.atac,
  features = c("blacklist_fraction", 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 2
)


high_quality <- WhichCells(alldata.atac, expression = 
                             nCount_ATAC > 1000 & 
                             TSS.enrichment >4 &
                             nucleosome_signal <= 2 &
                             blacklist_fraction < 0.002)
alldata.atac <- subset(alldata.atac, cells = high_quality)


# Assign donors
alldata.atac$donor <- 'Unclear'
alldata.atac$donor[alldata.atac$genotype=='donor0'] <- 'H12'
alldata.atac$donor[alldata.atac$genotype=='donor1'] <- 'H9'
alldata.atac <- subset(alldata.atac, idents = c('Unclear', 'doublet'), invert = T)


# Normalize ADT
alldata.atac <- NormalizeData(alldata.atac, assay = 'ADT', normalization.method = 'CLR', verbose = F)

# Import mito data
mito.data <- ReadMGATK(paste0(import_path, 'mgatk/final/'))
alldata.mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cells present in ATAC
alldata.mito <- subset(alldata.mito, cells = colnames(alldata.atac))

# Add to ATAC object
alldata.atac[["mito"]] <- alldata.mito
alldata.atac <- AddMetaData(alldata.atac, metadata = mito.data$depth, col.name = "mtDNA_depth")
rm(alldata.mito)


### Dimreduction, clustering ###
alldata.atac <- RunTFIDF(alldata.atac) %>%
  FindTopFeatures() %>%
  RunSVD() 

alldata.atac <- RunUMAP(alldata.atac, reduction = "lsi", dims = 2:30)
DimPlot(alldata.atac, group.by = "donor")

alldata.atac <- FindNeighbors(alldata.atac, reduction= "lsi", dims = 2:30, verbose = T)
alldata.atac <- FindClusters(alldata.atac, graph.name = 'ATAC_snn', resolution = 0.4, algorithm = 1, group.singletons = T)
DimPlot(alldata.atac, label = T)

# Add gene scores
gene.activities <- GeneActivity(alldata.atac)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
alldata.atac[['RNA_atac']] <- CreateAssayObject(counts = gene.activities)
alldata.atac <- NormalizeData(
  object = alldata.atac,
  assay = 'RNA_atac',
  normalization.method = 'LogNormalize',
  scale.factor = median(alldata.atac$nCount_RNA_atac)
)

# Exclude contaminations
DefaultAssay(alldata.atac) <- "ADT"
FeaturePlot(alldata.atac, features = c("adt_CD56", "CD5", "NKG2C-Strep",
                                       "KIR2DL1/S1/S3/S5", "KIR2DL2/L3",
                                       "KIR3DL1", "CD161", "NKp30", "CD85j"),
            min.cutoff = "q1", max.cutoff = "q99", order = T,)&NoLegend()&NoAxes()
FeaturePlot(alldata.atac, features = c("NCR1", "IL2RA", "IL7R", "RORC"),
            min.cutoff = "q1", max.cutoff = "q99", order = T)&NoLegend()&NoAxes()
NK.atac <- subset(alldata.atac, idents = c("8", "9"), invert = T)
NK.atac$coarse <- Idents(NK.atac)

# Call peaks
DefaultAssay(NK.atac) <- "ATAC"
peaks_MACS3 <- CallPeaks(NK.atac, group.by = "coarse",
                         macs2.path = "~/work/bin/miniforge3/envs/R_4.3.3/bin/macs3")


# Count Fragments

MACS3_counts <- FeatureMatrix(Fragments(NK.atac), features = peaks_MACS3)

NK.atac[["MACS3"]] <- CreateAssayObject(MACS3_counts)

# Separate donors
H12.atac <- subset(NK.atac, subset = donor == "H12")
H9.atac <- subset(NK.atac, subset = donor == "H9")


### H12
DefaultAssay(H12.atac) <- "MACS3"
H12.atac@reductions$lsi <- NULL
H12.atac <- RunTFIDF(H12.atac) %>%
  FindTopFeatures(min.cutoff = "q10") %>%
  RunSVD() 
H12.atac <- RunUMAP(H12.atac, reduction = "lsi", dims = 2:30)
DimPlot(H12.atac)

H12.atac <- FindNeighbors(H12.atac, reduction= "lsi", dims = 2:30, verbose = T)
H12.atac <- FindClusters(H12.atac, resolution = 0.4)
DimPlot(H12.atac, label = T)

DefaultAssay(H12.atac) <- "ADT"
FeaturePlot(H12.atac, features = c( "CD5", "NKG2C-Strep", "NKG2A", "CD2",
                                    "KIR2DL1/S1/S3/S5", "KIR2DL2/L3",
                                    "KIR3DL1",  "NKp30", "CD85j"),
            cols = c("grey", "red"),
            min.cutoff = "q5", max.cutoff = "q95", order = T)&NoLegend()&NoAxes()



# Cluster 2 contains NKG2A+ and KIR2DL1+ clones => subcluster
DefaultAssay(H12.atac) <- "MACS3"
temp <- subset(H12.atac, idents = "2")
FeaturePlot(temp, features = c("KIR2DL2/L3","CD57", "KLRG1", "NKG2A", "CD2", "NKG2C-Strep"),
            cols = c("grey", "red"),
            min.cutoff = "q5", max.cutoff = "q95", order = T)&NoLegend()&NoAxes()

temp <- FindNeighbors(temp, reduction = "lsi", dims = 2:30, k.param = 10)
temp <- FindClusters(temp, resolution = 0.2)
DimPlot(temp)
Idents(H12.atac, cells = WhichCells(temp, ident = "1")) <- "2a"

# Can't really resolve NKG2A+ cells as separate cluster
# but NKG2A+ NKG2C- NKp30- expression pattern is consistent with Tapestri
DefaultAssay(temp) <- "ADT"
FeatureScatter(temp, feature1 = "NKG2A", feature2 = "NKG2C-Strep", slot= "counts")
FeatureScatter(temp, feature1 = "NKG2A", feature2 = "KIR2DL2/L3")
FeatureScatter(temp, feature1 = "NKG2A", feature2 = "NKp30")
DimPlot(H12.atac)

### H9
DefaultAssay(H9.atac) <- "MACS3"
H9.atac@reductions$lsi <- NULL
H9.atac <- RunTFIDF(H9.atac) %>%
  FindTopFeatures(min.cutoff = "q10") %>%
  RunSVD() 
H9.atac <- RunUMAP(H9.atac, reduction = "lsi", dims = 2:30)
DimPlot(H9.atac)

DefaultAssay(H9.atac) <- "ADT"
FeaturePlot(H9.atac, features = c( "KLRG1", "NKG2C-Strep", "NKG2A", "CD2",
                                   "KIR2DL1/S1/S3/S5", "KIR2DL2/L3",
                                   "KIR3DL1",  "NKp30", "CD161"),
            cols = c("grey", "red"),
            min.cutoff = "q1", max.cutoff = "q99", order = T)&NoLegend()&NoAxes()
FeaturePlot(H9.atac, c("CADM1", "JAKMIP1", "ZBTB16"),
            cols = c("darkblue", "yellow"),
            min.cutoff = "q1", max.cutoff = "q99", order = T)&NoLegend()&NoAxes()
DefaultAssay(H9.atac) <- "MACS3"
H9.atac <- FindNeighbors(H9.atac, reduction= "lsi", dims = 2:30, verbose = T)
H9.atac <- FindClusters(H9.atac, resolution = 0.2)
DimPlot(H9.atac, label = T)

# Subcluster bright and small adaptive cluster included with dim
temp <- subset(H9.atac, idents = 1)
temp <- FindNeighbors(temp, reduction= "lsi", dims = 2:30, verbose = T)
temp <- FindClusters(temp, resolution = 0.4, algorithm = 1, group.singletons = T)
DimPlot(temp)
Idents(H9.atac, cells = WhichCells(temp, idents = "2")) <- "Bright"
Idents(H9.atac, cells = WhichCells(temp, idents = "3")) <- "3a"
markers <- FindMarkers(H9.atac, ident.1 = "3a", ident.2 = "2", assay = "RNA_atac")
DimPlot(H9.atac)



### Mitochondrial genotyping ###
# Stash cluster idents
H9.atac$annotation <- Idents(H9.atac)
H12.atac$annotation <- Idents(H12.atac)

# H9
VlnPlot(H9.atac, features = "mtDNA_depth")
table(H9.atac$mtDNA_depth > 5)
H9.atac <- subset(H9.atac, subset = mtDNA_depth > 5)

# Calling of variable sites and selection of variants with high confidence
variable.sites_H9 <- IdentifyVariants(H9.atac, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites_H9, concordance.threshold = 0.65, vmr.threshold = 0.001, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
# Variable features
high.conf_H9 <- subset(
  variable.sites_H9, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.001
)

# Compute variant frequency per cell
H9.atac <- AlleleFreq(
  object = H9.atac,
  variants = c(high.conf_H9$variant),
  assay = "mito"
)




# H12
VlnPlot(H12.atac, features = "mtDNA_depth")
table(H12.atac$mtDNA_depth > 5)
H12.atac <- subset(H12.atac, subset = mtDNA_depth > 5)

# Calling of variable sites and selection of variants with high confidence
variable.sites_H12 <- IdentifyVariants(H12.atac, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites_H12, concordance.threshold = 0.65, vmr.threshold = 0.01, min.cells = 5)

# Establish a filtered data frame of variants based on this processing
# Variable features
high.conf_H12 <- subset(
  variable.sites_H12, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

# Compute variant frequency per cell
H12.atac <- AlleleFreq(
  object = H12.atac,
  variants = c(high.conf_H12$variant, "991G>A", "12889G>A"),
  assay = "mito"
)

saveRDS(H9.atac, file = "H9_atac.rds")
saveRDS(H12.atac, file = "H12_atac.rds")
