library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(hdf5r)
source(file = "mgatk2_functions.R")


# Read in ATAC
seu <- readRDS("../ATAC/H12_atac.rds")

# Plot Protein markers parallel to Tapestri
DefaultAssay(seu) <- "ADT"
p1 <- FeaturePlot(seu, features = c("NKp30", "NKG2C-Strep"),ncol = 2,
                  min.cutoff = "q1", max.cutoff = "q95", order = T, coord.fixed = T)&
  scale_color_gradientn(colors = magma(9, alpha = 1, begin = 0, end = 1, direction = 1))&
  NoAxes()&NoLegend()


# Select informative variants
mito_muts <- rownames(seu[["alleles"]])[c(1,3,9,10,11,12,
                                          14,22,23,26,27,
                                          29,37,40,44,56)]
temp <- subset(seu, idents = c("0", "3"), invert = T)
DoHeatmap(temp, features = mito_muts, slot = "data",disp.min = 0, disp.max = 80)&scale_fill_gradient(low = "grey", high = "red")


# Get informative mito variants included in Tapestri
mito_variants <- read.csv("mito_whitelist.csv")
mito_variants <- mito_variants[,2]


result <- gsub("/", ">", sub("^chrM:", "", mito_variants))
result <- sub(":", "", result)
mito_h12 <- result[result%in%rownames(seu[["alleles"]])]

DefaultAssay(seu) <- "alleles"

# Load Tapestri
seu_tapestri <- readRDS("HC12/HC12_tapestri.rds")
mgatk_data <- read_mgatk_hdf5('HC12/mgatk2')
vcf <- read.table("mito_variants.vcf", 
                  comment.char = "#", 
                  header = FALSE,
                  col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))



variant_positions <- as.integer(sub("^([0-9]+).*", "\\1", 
                                    mito_h12))


variants <- identify_variants(mgatk_data)


variants_H12 <- variants %>% dplyr::filter(variant%in%mito_muts)

vafs <- calculate_allele_freq(mgatk_data, variants = variants_H12)
mito_assay <- CreateAssayObject(data = vafs)
mito_assay <- subset(mito_assay, cells = Cells(seu_tapestri))
seu_tapestri[["mito"]] <- mito_assay


# Demonstrate clonal specificity of selected variants
p1 <- VlnPlot(seu_tapestri, features = rownames(seu_tapestri[["mito"]]), ncol = 4)


# Assign variants to clones
mito_A <- c("12868G>A", "684T>C","5356T>C", "4848G>A","15246G>A", "16242C>T")
mito_B <- c("513G>A")


## Look for those variants in ATAC 
DefaultAssay(seu) <- "alleles"
VlnPlot(seu, features = mito_A)
cells_cloneA <- WhichCells(seu, expression = `12868G>A` > 0.75 | `684T>C` > 0.75 |
                             `4848G>A` >0.75 | `5356T>C` > 0.5 | `15246G>A` > 0.75)


VlnPlot(seu, features = mito_B)
cells_cloneB <- WhichCells(seu, expression = `513G>A` > 0.75)


p3 <- FeaturePlot(seu, features = c(mito_A, mito_B), ncol = 3,
                  min.cutoff = "q1", max.cutoff = "q99", order = T, coord.fixed = T)&
  scale_color_gradient(low = "grey", high = "darkred")&
  NoAxes()&NoLegend()


# Plot mito-clones on ATAC
p4 <- DimPlot(seu, cells.highlight = cells_cloneA, cols.highlight = "darkblue")&NoLegend()&NoAxes()&coord_fixed()
p5 <- DimPlot(seu, cells.highlight = cells_cloneB, cols.highlight = "darkorange")&NoLegend()&NoAxes()&coord_fixed()

p6 <- p4+p5


## Subset clones to analyze clone-specific features
seu_clones <- subset(seu, cells = c(cells_cloneA, cells_cloneB))
Idents(seu_clones, cells = cells_cloneA) <- "HC12_A"
Idents(seu_clones, cells = cells_cloneB) <- "HC12_B"

seu_clones$clone <- Idents(seu_clones)
seu_clones$clone <- factor(seu_clones$clone, levels = c("HC12_A", "HC12_B"))
Idents(seu_clones) <- seu_clones$clone
DefaultAssay(seu_clones) <- "ATAC"


DefaultAssay(seu_clones) <- "ATAC"
p1 <- CoveragePlot(seu_clones, region = "KIR2DL1", extend.upstream = 1000,
                   expression.assay = "ADT", features = "KIR2DL1/S1/S3/S5", links = F)&
  scale_fill_manual(values = c("darkblue", "darkorange"))
p2 <- CoveragePlot(seu_clones, region = "KIR2DL3", extend.upstream = 1000,, links = F,
                   extend.downstream= 1000, expression.assay = "ADT", features = "KIR2DL2/L3")&
  scale_fill_manual(values = c("darkblue",  "darkorange"))
