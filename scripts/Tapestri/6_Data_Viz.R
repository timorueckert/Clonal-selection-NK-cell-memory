# Tapestri Data Viz
path <- "HC09_tapestri.rds"
# path <- "HC12_tapestri.rds"
# path <- "HC02_tapestri.rds"
clone_colors <- c("darkgrey","darkblue", "darkorange", "darkgreen", "lightpink", "gold1")

# Plot mutations
seu_tapestri <- readRDS(path)
Idents(seu_tapestri) <- seu_tapestri$ordered
DefaultAssay(seu_tapestri) <- "DNA_missing"
whitelist <- VariableFeatures(seu_tapestri)
mutations <- FindAllMarkers(seu_tapestri, only.pos = T, features = whitelist)

# Keep cell ordering for plotting of subclones within larger families
cell.order <- names(seu_tapestri$fine[order(seu_tapestri$fine)])

Idents(seu_tapestri) <- seu_tapestri$coarse

# Plot reordered mutations of clonal families
p1 <- DoHeatmap(seu_tapestri, features = mutations$gene, slot = "data",
                disp.min = -50, disp.max = 100, cells = cell.order, group.colors = clone_colors)+
  scale_fill_gradient2(low = "black",
                       mid = "#3B4D73",
                       high = "#D7ECEE",
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Plot proteins
DefaultAssay(seu_tapestri) <- "ADT"
seu_tapestri <- ScaleData(seu_tapestri,do.center = T,scale.max = 2)
p2 <- DoHeatmap(seu_tapestri, features = c("NKp30","NKG2C", "NKG2A", "CD57", "CD2",
                                     "KIR2DL2-L3-S2", "KIR2DL1-S1-S3-S5",
                                     "KIR3DL1", "CD161", "KLRG1"),
          disp.min = 0, disp.max = 2, cells = cell.order,
          slot = "scale.data")+
  scale_fill_gradientn(colors = magma(9, alpha = 1, begin = 0, end = 1, direction = 1))


#### Combined analyses of all donors ####

# Import and merge objects
seu_h9 <- readRDS("HC09_tapestri.rds")
seu_h12 <- readRDS("HC12_tapestri.rds")
seu_h2 <- readRDS("HC02_tapestri.rds")

seu_h9$donor <- "h9"
seu_h12$donor <- "h12"
seu_h2$donor <- "h2"

seu_h9$coarse <- factor(seu_h9$coarse, levels = c(levels(seu_h9$coarse), "polyclonal"))
seu_h9$coarse[seu_h9$coarse == "non-clones"] <- "polyclonal"

seu_h12$coarse <- factor(seu_h12$coarse, levels = c(levels(seu_h12$coarse), "polyclonal"))
seu_h12$coarse[seu_h12$coarse == "non-clones"] <- "polyclonal"

seu_tapestri <- merge(seu_h12, c(seu_h2, seu_h9))
seu_tapestri$coarse <- paste0(seu_tapestri$donor, seu_tapestri$coarse)

# Scale ADT assay across donors
DefaultAssay(seu_tapestri) <- "ADT"
seu_tapestri <- ScaleData(seu_tapestri, split.by = "donor", do.center = T, scale.max = 2)


# Clone size distribution for all three donors
df <- as.data.frame(table(seu_tapestri$coarse, seu_tapestri$donor))
colnames(df) <- c("clone", "donor", "n")

df <- df %>%
  dplyr::filter(n > 0)
df <- df %>%
  mutate(
    clone = factor(
      clone,
      levels = c(
        "h12polyclonal","h12A", "h12B", "h12C", "h12D", "h12E", 
        "h2polyclonal","h2A", "h2B", "h2C", "h2D", 
        "h9polyclonal","h9A", "h9B"
      )
    )
  )
ggplot(df, aes(x = donor, y = n, fill = clone)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Clone Size (%)",
  )+theme_classic()+
  scale_fill_manual(values=clone_colors[c(c(1:6),
                                          c(1:5),
                                          c(1:3))])

# Calculate cell size estimation
# Import flow data
WES_flow <- read.csv("flow_nk_freqs.csv", sep = ";", dec = ",")
colnames(WES_flow)[1] <- "donor"
head(WES_flow)

# Filter for three analyzed donors
population_size <- WES_flow %>% dplyr::filter(donor %in% c("HC02", "HC09", "HC12"))
# Convert NK % into frequency and put into df
NK_freq <- population_size$NK/100
names(NK_freq) <- population_size$donor
NK_freq <- data.frame(donor = names(NK_freq), NK_freq = NK_freq)

# Harmonize donor names
df$donor <- as.character(df$donor)
df$donor[df$donor=="h12"] <- "HC12"
df$donor[df$donor=="h2"] <- "HC02"
df$donor[df$donor=="h9"] <- "HC09"

# Merge dfs
clone_freq <- merge(df, NK_freq)

# Calculate clone freq NK
clone_freq <- clone_freq %>%
  group_by(donor) %>%
  mutate(clone_freq = n / sum(n)) %>%
  ungroup()

# Multiply with NK frequency to get frequency within lymphocytes
clone_freq <- clone_freq %>%
  group_by(donor) %>%
  mutate(clone_freq_ly = clone_freq * NK_freq) %>%
  ungroup()

# Multiply with adult human lymphocyte count to get absolute size estimate
clone_freq$clone_est <- clone_freq$clone_freq_ly*10^10

# Plot
conversion_factor <- 1e10

plot3 <- clone_freq %>% filter(!clone %in% c("h2polyclonal", "h9polyclonal", "h12polyclonal")) %>% 
ggplot(aes(x = donor, y = clone_freq_ly, fill = clone)) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    name = "Clone Frequency (% of Lymphocytes)",
    labels = label_percent(accuracy = 0.5),
    breaks = c(0, 2.5, 5, 7.5, 10) / 100,
    expand = c(0,0),
    sec.axis = sec_axis(
      ~ . * conversion_factor,
      name = "Estimated cell number",
      labels = label_number(scale = 1e-8, suffix = "x 10^8"),
      breaks = c(0, 2.5, 5, 7.5, 10) * 1e8
    )
  ) +
  labs(x = "Donor", fill = "Clone")+
  theme_classic()+
  scale_fill_manual(values=clone_colors[c(c(2:6),
                                          c(2:5),
                                          c(2:3))])&NoLegend()

### Comparison with clone sizes within NKG2C+ estimated by WES ###
DefaultAssay(seu_tapestri) <- "ADT"
VlnPlot(seu_tapestri, features = c("NKG2C", "NKG2A"),
        group.by = "donor", pt.size = 0.1)

# Subset NKG2C+
seu_2c <- subset(seu_tapestri, subset = NKG2C > 0.6 & NKG2A < 0.5)

tab <- table(seu_2c$coarse, seu_2c$donor)

freq_df <- as.data.frame(tab) %>%
  dplyr::rename(
    clone = Var1,
    donor = Var2,
    n = Freq
  ) %>%
  group_by(donor) %>%
  mutate(
    N = sum(n),
    frequency = n / N,
    
    # simple binomial standard error
    se = sqrt(frequency * (1 - frequency) / N),
    
    # Jeffreys 95% confidence interval
    ci_low = qbeta(0.025, n + 0.5, N - n + 0.5),
    ci_high = qbeta(0.975, n + 0.5, N - n + 0.5)
  ) %>%
  ungroup() %>%
  dplyr::filter(n > 0)

freq_df

# Import scifer data and combine with Tapestri results
scifer_df <- read.csv("scifer_df.csv")

freq_df <- data.frame(Tapestri = freq_df[freq_df$clone%in%c("h12A", "h2A", "h9A"),])
freq_df$donor <- c("HC12", "HC02", "HC09")
plot_df <- merge(freq_df, scifer_df, .by = "donor")

p1 <- plot_df %>% 
  ggplot(aes(x = f*100, y = Tapestri.frequency*100))+
  geom_abline(slope = 1, intercept = 0,linetype = "dashed", color = "grey")+
  geom_errorbar(
    aes(ymin = Tapestri.ci_low*100,
        ymax = Tapestri.ci_high*100),
    width= 0.01)+
  geom_errorbarh(
    aes(xmin = f_low*100,
        xmax = f_upp*100),
    height = 0.01
  )+
  geom_point(color = "darkred", size =1)+
  
  ylim(c(0,100))+xlim(c(0,100))+
  xlab("Frequency Estimate WES (%)")+
  ylab("Frequency Estimate Tapestri (%)")+
  theme_classic()
