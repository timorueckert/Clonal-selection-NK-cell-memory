##### Visualization of mutation loads and VAF distributions

library(dplyr)
library(ggplot2)

# include only on-target mutations associated
# with annotated genes to enable better comparability
# of mutation load analysis irrespective of capture efficiency

muts_exome <- muts_filtered %>% 
  filter(on_target) %>% 
  filter(!is.na(SYMBOL))

# Plot number of mutations in the two populations
df <- muts_exome %>% 
  mutate(pass = tumor_af > 0.05, on_target) %>% 
  group_by(donor_id, cell_type, cohort) %>% 
  dplyr::summarise(n_muts = sum(pass))


plot1 <- df %>% 
  ggplot(aes(y = n_muts, x = cell_type, fill = cell_type))+
  geom_jitter(width = 0.2, size = 2, shape = 21, color = "black")+
  scale_fill_manual(values = c("#3e8bca", "darkred"))+
  ylab("Number of Variants")+
  theme_classic()+
  theme(text = element_text(size=16), axis.text.x= element_text(size = 14, color = "black"),
        axis.text.y= element_text(size = 14, color = "black"),
        axis.title.x =  element_text(size = 0))


# Plot mutations histogram for one representative donor and all donors (supplement)
library(ggpubr)
library(cowplot)

# Plot > 0.03 for better visualization of clonal peaks
filtered_df <- muts_exome %>% 
  dplyr::filter(!chrom%in%c("chrX", "chrY", "chrM"), tumor_af > 0.03)


# Include only memory and conventional
filtered_df$cell_type <- factor(filtered_df$cell_type,
                                levels = c("NKG2C", "NKG2A"))

# Reorder donors so that the ones with paired 
# are plotted first conventional populations
donor_levels <- unique(filtered_df$donor_id)
donor_levels_sorted <- donor_levels[c(1:10,66:68,81,24)]
donor_levels_sorted <- unique(c(donor_levels_sorted, donor_levels))
filtered_df$donor_id <- factor(filtered_df$donor_id,
                               levels = donor_levels_sorted)

# Plot all
plot2 <- filtered_df %>% ggplot(aes(x = tumor_af, fill = cell_type))+
  geom_histogram(color = "black", binwidth = 0.02)+
  scale_fill_manual(values = c( "darkred", "#3e8bca"))+
  theme_classic()+
  xlab("Variant Allele Frequency")+
  coord_cartesian(xlim = c(0.03, 0.5)) +
  scale_y_continuous(expand = c(0,0))+
  ylab("Count")+
  theme(text = element_text(size=16), axis.text.x= element_text(size = 14, color = "black"),
        axis.text.y= element_text(size = 14, color = "black"))+
  facet_wrap(~donor_id)


# Plot representative
plot3 <- filtered_df %>% dplyr::filter(donor_id == "HC02") %>%  ggplot(aes(x = tumor_af, fill = cell_type))+
  geom_histogram(color = "black", binwidth = 0.02)+
  scale_fill_manual(values = c( "darkred", "#3e8bca"))+
  theme_classic()+
  xlab("Variant Allele Frequency")+
  coord_cartesian(xlim = c(0.03, 0.5)) +
  # scale_y_continuous(expand = c(0,0), breaks = c(0,5,10,15),
  #                    limits = c(0,15))+
  ylab("Count")+
  theme(text = element_text(size=16), axis.text.x= element_text(size = 14, color = "black"),
        axis.text.y= element_text(size = 14, color = "black"))&NoLegend()



#### Plotting of scifer results and immunephenotyping data ####

# Import Matthias' scifer results
scifer <- read.csv("data/WES/input/SciferResults WES vafmin=0.05.csv")

# Import flow data and merge
WES_flow <- read.csv("data/WES/input/flow_nk_freqs.csv", sep = ";", dec = ",")
colnames(WES_flow)[1] <- "donor"
head(WES_flow)

clone_df <- merge(scifer, WES_flow)

# Adjust for percentages
clone_df$clone_freq <- clone_df$f * clone_df$NK/100 * clone_df$NKG2C/100
clone_df$clone_freq_lo <- clone_df$f_low * clone_df$NK/100 * clone_df$NKG2C/100
clone_df$clone_freq_hi <- clone_df$f_upp* clone_df$NK/100 * clone_df$NKG2C/100

# Define donors with clearly quantifiable leading clone (size > 0.25)
clone_df <- clone_df %>% 
  mutate(clonal = ifelse(clone_df$f>0.25, 
                         "clonal",
                         "poly/oligoclonal"))
table(clone_df$clonal)

# Export
write.csv(clone_df, "data/WES/outs/scifer_df.csv")


# Plot clonality criterion
p1 <- ggplot(clone_df, aes(x = f*100, y = M)) +
  geom_point(aes(fill = clonal),
             shape = 21, color = "black",
             size = 2
  ) +
  theme_classic() +
  labs(
    y = "No. of founder mutations",
    x = "Clone size within memory NK cells (%)"
  ) +
  theme(
    text = element_text(size = 16),
  )+
  scale_fill_manual(values = c("darkred","grey"))


# Clone size within memory NK
# Use the same jitter position for points and error bars
jitter_pos <- position_jitter(width = 0.15, height = 0, seed = 123)

p2 <- ggplot(clone_df, aes(y = f, x = "Donors")) +
  geom_errorbar(
    aes(ymin = f_low, ymax = f_upp),
    width = 0,
    position = jitter_pos
  ) +
  geom_point(aes(fill = clonal),
             shape = 21, color = "black",
             position = jitter_pos,
             size = 2
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Clone size within memory NK cells (%)"
  ) +
  theme(
    text = element_text(size = 16),
  )+
  scale_fill_manual(values = c("darkred","grey"))



# Clone Size within Lymphocytes
p3 <- ggplot(clone_df, aes(y = clone_freq, x = "Donors")) +
  geom_errorbar(
    aes(ymin = clone_freq_lo, ymax = clone_freq_hi),
    width = 0,
    position = jitter_pos
  ) +
  geom_point(aes(fill = clonal),
             shape = 21, color = "black",
             position = jitter_pos,
             size = 2
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Clone size within lymphocytes (%)"
  ) +
  theme(
    text = element_text(size = 16),
  )+
  scale_fill_manual(values = c("darkred","grey"))


# NKG2C frequency
p5 <- ggplot(clone_df, aes(y = NKG2C, x = f*100)) +
  geom_errorbarh(
    aes(xmin = f_low*100, xmax = f_upp*100),
  ) +
  geom_point(aes(fill = clonal),
             shape = 21, color = "black",
             size = 2
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "NKG2C+ NK cells (% of NK cells)"
  ) +
  theme(
    text = element_text(size = 16),
  )+
  scale_fill_manual(values = c("darkred","grey"))



# Plot cumulative distributions
p8 <- ggplot(clone_df, aes(x = f)) +
  stat_ecdf(geom = "step", linewidth = 0.8, color = "darkred") +
  scale_x_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::label_percent(),
                     breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Clone frequency (% of memory NK)",
    y = "Cumulative proportion of donors"
  ) +
  theme_classic(base_size = 12)

p9 <- ggplot(clone_df, aes(x = clone_freq)) +
  stat_ecdf(geom = "step", linewidth = 0.8, color = "darkred") +
  scale_x_continuous(
    labels = scales::label_percent(),
    limits = c(0, 0.15)
  ) +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::label_percent(),
                     breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Clone frequency (% of lymphocytes)",
    y = "Cumulative proportion of donors"
  ) +
  theme_classic(base_size = 12)