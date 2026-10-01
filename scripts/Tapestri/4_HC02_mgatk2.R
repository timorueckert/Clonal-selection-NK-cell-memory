library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
source(file = "mgatk2_functions.R")

# From GreenleafLab/ArchR
blueyellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2",
                "5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")


### scATAC preprocessing
seu <- readRDS("../ATAC/CMVpos1_full")

# Subset to mtASAP5
seu <- subset(seu, experiment == "mtASAP5")

# Update Fragment File path
temp <- Fragments(seu)[[25]]
temp <- UpdatePath(temp, "../ATAC/mtASAP5/fragments.tsv.gz")
Fragments(seu) <- NULL
Fragments(seu) <- temp
rm(temp)


# Plot Protein markers parallel to Tapestri
DefaultAssay(seu) <- "ADT"
p1 <- FeaturePlot(seu, features = c("NKp30", "Anti-PE"),ncol = 2,
                  min.cutoff = "q1", max.cutoff = "q95", order = T, coord.fixed = T)&
  scale_color_gradientn(colors = magma(9, alpha = 1, begin = 0, end = 1, direction = 1))&
  NoAxes()&NoLegend()


# Get informative mito variants included in Tapestri
# Select informative variants
DefaultAssay(seu) <- "alleles"
FeaturePlot(seu, features = rownames(seu),
            cols = c("grey", "darkred"),
            order = T)&NoLegend()&NoAxes()
asap_variants <- rownames(seu[["alleles"]])[c(1,2,3,6,7,9,11,15)]
mito_variants <- read.csv("mito_whitelist.csv")
mito_variants <- mito_variants[,2]
mito_variants <- c(mito_variants, "chrM:7395T>G")

result <- gsub("/", ">", sub("^chrM:", "", mito_variants))
result <- sub(":", "", result)
mito_h02 <- result[result%in%asap_variants]



# Load Tapestri
seu_tapestri <- readRDS("HC02/HC02_tapestri.rds")
mgatk_data <- read_mgatk_hdf5('HC02/mgatk2')
vcf <- read.table("mito_variants.vcf", 
                  comment.char = "#", 
                  header = FALSE,
                  col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))



variant_positions <- as.integer(sub("^([0-9]+).*", "\\1", 
                                    mito_h02))


variants <- identify_variants(mgatk_data)


variants_H02 <- variants %>% dplyr::filter(variant%in%mito_h02)
manual_variants <- data.frame(7395,"T>G", "7395T>G", 0,0,0,0,0)
colnames(manual_variants) <- colnames(variants_H02)
variants_H02 <- rbind(variants_H02, manual_variants)


vafs <- calculate_allele_freq(mgatk_data, variants = variants_H02)
mito_assay <- CreateAssayObject(data = vafs)
mito_assay <- subset(mito_assay, cells = Cells(seu_tapestri))
seu_tapestri[["mito"]] <- mito_assay



# Assign variants to clones
p1 <- VlnPlot(seu_tapestri, features = rownames(seu_tapestri[["mito"]]), ncol = 4)

mito_A <- c("13710A>G")
mito_C <- c("15221G>A", "6565A>G")

# "1485G>A" is subclone of larger "7395T>G" variant => combine
FeatureScatter(seu,  "7395T>G", "1485G>A")
mito_D <- c("1485G>A", "7395T>G")


## Use these variants to define clones in ATAC 
DefaultAssay(seu) <- "alleles"
VlnPlot(seu, features = mito_A)
cells_cloneA <- WhichCells(seu, expression = `13710A>G` > 0.75)

VlnPlot(seu, features = mito_D)
cells_cloneD <- WhichCells(seu, expression = `1485G>A` > 0.2| `7395T>G` > 0.4)


p3 <- FeaturePlot(seu, features = c("13710A>G", "1485G>A"), ncol = 3,
            min.cutoff = "0.05", max.cutoff = "0.5", order = T, coord.fixed = T)&
  scale_color_gradient(low = "grey", high = "darkred")&
  NoAxes()&NoLegend()


# Plot mito-clones on ATAC
p4 <- DimPlot(seu, cells.highlight = cells_cloneA, cols.highlight = "darkblue")&NoLegend()&NoAxes()&coord_fixed()
p6 <- DimPlot(seu, cells.highlight = cells_cloneD, cols.highlight = "lightpink")&NoLegend()&NoAxes()&coord_fixed()

# Plot phenotypic markers
p1 <- FeaturePlot(seu, features = c("Anti-PE", "ADT_KIR2DL1-S1-S3-S5", "ADT_KIR3DL1"),
                  order = T, min.cutoff = "q5", max.cutoff = "q95", ncol = 3
)&NoLegend()&NoAxes()&
  scale_color_viridis(option = "magma")&coord_fixed()


# Subset for coverage plots
seu_clones <- subset(seu, cells = c(cells_cloneA,cells_cloneD))
Idents(seu_clones, cells = cells_cloneA) <- "A"
Idents(seu_clones, cells = cells_cloneD) <- "D"
seu_clones$clone <- Idents(seu_clones)

seu_clones$clone <- factor(seu_clones$clone, levels = c("A", "D"))
Idents(seu_clones) <- seu_clones$clone
DefaultAssay(seu_clones) <- "MACS2"

p1 <- CoveragePlot(seu_clones, region = "KIR2DL1", extend.upstream = 1000,
                   expression.assay = "ADT", features = "KIR2DL1-S1-S3-S5", links = F)&
  scale_fill_manual(values = c("darkblue", "lightpink"))
p2 <- CoveragePlot(seu_clones, region = "chr19-54816468-54820000", extend.upstream = 1000,, links = F,
                   extend.downstream= 10000, expression.assay = "ADT", features = "KIR3DL1")&
  scale_fill_manual(values = c("darkblue",  "lightpink"))



