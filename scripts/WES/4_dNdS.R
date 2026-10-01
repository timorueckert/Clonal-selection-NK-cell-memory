#dN/dS analysis
library(dndscv)
library(dplyr)
library(Biobase)

# Import data
input <- read.csv("data/WES/outs/2609_full_cohort_filtered.csv")

# Filter for mutations with VAF > 0.06 in memory population
df <- input %>% dplyr::filter(cell_type == "NKG2C", tumor_af > 0.06)


# Re-format for dN/dS
input_df <- data.frame(sampleID = df$donor_id, chr = df$chrom,
                       pos = df$pos, ref = df$ref, mut = df$alt)
input_df$chr <- sub("^chr", "", input_df$chr)

input_df$chr <- factor(input_df$chr, levels = unique(input_df$chr))
head(input_df)

# Remove problematic sites:
# Mutations observed in contiguous sites within a sample
adjacent_positions <- input_df %>%
  arrange(sampleID, chr, pos) %>%
  group_by(sampleID, chr) %>%
  filter(coalesce(pos - lag(pos) == 1, FALSE) | coalesce(lead(pos) - pos == 1, FALSE)) %>%
  ungroup()

# Same mutations observed in different sampleIDs
recurrent_across_donors <- input_df %>%
  add_count(chr, pos, ref, mut, name = "n_donors") %>%
  filter(n_donors > 1) %>%
  select(-n_donors)

input_df <- input_df %>%
  anti_join(adjacent_positions, by = c("sampleID", "chr", "pos")) %>%
  anti_join(recurrent_across_donors, by = c("chr", "pos", "ref", "mut"))

message(sprintf(
  "dN/dS input: %d SNVs (dropped %d in contiguous pairs, %d recurrent across donors)",
  nrow(input_df), nrow(adjacent_positions), nrow(recurrent_across_donors)))

# Export for reproducibility
write.csv(input_df, "data/WES/outs/2609_dNdS_input.csv")

## Import cancer genes from CGC
cgc_genes <- read.csv("data/WES/input/Census_allFri Apr 11 18_29_11 2025.csv")
head(cgc_genes)
driver_genes <- cgc_genes$Gene.Symbol


# run analysis
dndsout_all = dndscv(input_df, refdb = "hg38",
                     outmats = T)

print(dndsout_all$globaldnds)

# Analysis on CGC genes
gene_list <- intersect(driver_genes, dndsout_all$genemuts$gene_name)
cgc_dnds <- genesetdnds(dndsout_all, gene_list = gene_list)
cgc_dnds

# Plot
plot_df1 <- dndsout_all$globaldnds
plot_df1$geneset <- "All"
plot_df1[1,1] <- "missense"
plot_df1[2,1] <- "nonsense"
plot_df1[3,1] <- "splice_site"
plot_df1[4,1] <- "truncating"
plot_df1[5,1] <- "all"

plot_df2 <- cgc_dnds$globaldnds_geneset[,c(1:3)]
plot_df2$geneset <- "Cancer Genes"
plot_df2$name <- c("missense",
                   "nonsense",
                   "splice_site",
                   "truncating",
                   "all")


plot_df <- rbind(plot_df1, plot_df2)

plot <- plot_df %>% filter(name %in%c("missense", "truncating", "all")) %>% arrange(name) %>% 
  ggplot(aes(x = name, y = mle, shape = geneset)) +
  geom_hline(yintercept = 1, linetype = "dotted")+
  geom_errorbar(aes(ymin = cilow, ymax = cihigh),
                width = 0, position = position_dodge(width = 0.4)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  xlab("") + ylab("dN/dS") +
  scale_fill_manual(values = c("black", "darkred"))+
  theme_classic()+
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"))&Seurat::NoLegend()
