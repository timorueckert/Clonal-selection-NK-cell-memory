## Phylogenetic analysis of subclones in infSCITE
seu_tapestri <- readRDS("data/Tapestri/outs/HC12_tapestri.rds")
# seu_tapestri <- readRDS("data/Tapestri/outs/HC02_tapestri.rds")
# seu_tapestri <- readRDS("data/Tapestri/outs/HC09_tapestri.rds")

# Extract NGT0 values for hierarchical clustering
Idents(seu_tapestri) <- seu_tapestri$fine
ngt_avg <- AverageExpression(seu_tapestri, features = VariableFeatures(seu_tapestri), assays = "NGT0", return.seurat = T)
ngt_avg_bin <- ngt_avg[["NGT0"]]@layers$data
rownames(ngt_avg_bin) <- rownames(ngt_avg)
colnames(ngt_avg_bin) <- colnames(ngt_avg)
# Binarize
ngt_avg_bin[ngt_avg_bin >0.8] <- 1L
ngt_avg_bin[ngt_avg_bin <0.8] <- 0

# Filter for clones and variants with at least 1 mutation
ngt_avg_bin_filtered <- ngt_avg_bin[rowSums(ngt_avg_bin)>0,]
#Export for analysis in infSCITE
ngt_bin_clean <- as.matrix(ngt_avg_bin_filtered)
rownames(ngt_bin_clean) <- NULL
colnames(ngt_bin_clean) <- NULL
ngt_bin_clean_muts <- rownames(ngt_avg_bin_filtered)
write.table(ngt_bin_clean, "data/Tapestri/outs/HC12_cluster_mutations.tsv",
            col.names = F, row.names = F, quote = F)
write.table(ngt_bin_clean_muts, quote = F, row.names = F, col.names = F,
            file = "data/Tapestri/outs/HC12_cluster.geneNames",
            sep = "\n")
dim(ngt_bin_clean)

# infSCITE prompt
# bash: ~/work/bin/infscite/infSCITE -i HC12_cluster_mutations.tsv -names HC12_cluster.geneNames -n 38 -m 18 -r 5 -l 300000 -fd 3.45e-3  -ad 0.1 -a -max_treelist_size 5


suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

# =========================================================
# 1. INPUT FILES
# =========================================================

gv_file <- "data/Tapestri/outs/HC12_cluster_mutations_ml0.gv"
newick_file <- "data/Tapestri/outs/HC12_cluster_mutations_ml0.newick"

# Import clone sizes
# Get Clone sizes from Seurat
clone_size <- table(seu_tapestri$fine)
clone_size <- as.matrix(clone_size)

# Calculate relative size
clone_size <- clone_size/sum(clone_size)
clone_fr <- clone_size[,1]

# infSCITE numbers samples from top to bottom (s_N) => add to clone size table
sample_id <- paste0("s_", 1:length(clone_fr))

# Add colors for plotting

# HC12
clone_colors <- c("darkgrey",
                  rep("darkblue", 10),
                  "darkorange",
                  rep("darkgreen", 2),
                  rep("lightpink", 2),
                  rep("gold1",2)
)
# # HC02,
# clone_colors <- c("darkgrey",
#                   rep("darkblue", 8),
#                   rep("darkorange", 3),
#                   rep("darkgreen", 2),
#                   rep("lightpink" , 3)
# )

# # HC09,
# clone_colors <- c("darkgrey",
#                   rep("darkblue", 18),
#                   rep("darkorange", 3)
# )
# Put into df
clone_anno <- data.frame(sample_label = sample_id,
                         clone_id = names(clone_fr),
                         clone_size = clone_fr,
                         clone_color = clone_colors)
# =========================================================
# 2. PARSE THE .GV FILE
# =========================================================
# We only extract:
# - mutation/internal node labels
# - sample attachment edges: <gv_node> -> s_<k>
#
# Assumption:
# - gv internal nodes are numbered 0,1,2,...
# - Newick labels are gv_node + 1
#   (as in your earlier example: Root=38 in gv becomes 39 in Newick)

parse_infscite_gv <- function(gv_file) {
  x <- readLines(gv_file, warn = FALSE)
  
  # all node definitions with labels
  node_lines <- x[str_detect(x, '^\\s*[^\\s]+\\[label=".*"\\];\\s*$')]
  
  node_df <- tibble(raw = node_lines) %>%
    transmute(
      node_id = str_match(raw, '^\\s*([^\\[]+)\\[label=')[, 2] |> str_trim(),
      label   = str_match(raw, 'label="(.*)"\\];\\s*$')[, 2]
    )
  
  # mutation/internal nodes
  mutation_nodes <- node_df %>%
    filter(!str_detect(node_id, "^s_\\d+$")) %>%
    mutate(
      gv_node = suppressWarnings(as.integer(node_id)),
      newick_label = as.character(gv_node + 1),
      gene = sub(":.*$", "", label)
    ) %>%
    arrange(gv_node)
  
  # sample nodes
  sample_nodes <- node_df %>%
    filter(str_detect(node_id, "^s_\\d+$")) %>%
    transmute(
      sample_node = node_id,
      sample_index0 = as.integer(str_extract(sample_node, "\\d+")),
      sample_index1 = sample_index0 + 1,
      sample_label = label
    ) %>%
    arrange(sample_index0)
  
  # edges
  edge_lines <- x[str_detect(x, '^\\s*[^\\s]+\\s*->\\s*[^\\s]+;\\s*$')]
  
  edge_df <- tibble(raw = edge_lines) %>%
    transmute(
      parent = str_match(raw, '^\\s*([^\\s]+)\\s*->')[, 2] |> str_trim(),
      child  = str_match(raw, '->\\s*([^;]+);\\s*$')[, 2] |> str_trim()
    )
  
  # sample attachments only: mutation node -> sample node
  attachments <- edge_df %>%
    filter(!str_detect(parent, "^s_"), str_detect(child, "^s_")) %>%
    transmute(
      gv_node = as.integer(parent),
      sample_node = child,
      newick_label = as.character(gv_node + 1)
    ) %>%
    left_join(sample_nodes, by = "sample_node") %>%
    left_join(
      mutation_nodes %>%
        select(gv_node, mut_label = label, gene),
      by = "gv_node"
    ) %>%
    arrange(sample_index0)
  
  list(
    mutation_nodes = mutation_nodes,
    sample_nodes = sample_nodes,
    attachments = attachments
  )
}

gv <- parse_infscite_gv(gv_file)

# =========================================================
# 3. CHECK ORDERING
# =========================================================

cat("Number of sample attachments in .gv:", nrow(gv$attachments), "\n")
cat("Number of rows in clone annotation table:", nrow(clone_anno), "\n")


# =========================================================
# 4. JOIN CLONE ANNOTATION BY ORDER
# =========================================================
# This is the key assumption:
# row 1 of clone_anno corresponds to s_1
# row 2 corresponds to s_2
# ...
#
# Since sample nodes in gv are s_0, s_1, ...
# and we created sample_index1 = sample_index0 + 1,
# sorting by sample_index0 aligns them to annotation row order.

# ---------------------------------------------------------
# Join clone annotations by order
# ---------------------------------------------------------
clone_map <- gv$attachments %>%
  mutate(order_index = row_number()) %>%
  bind_cols(clone_anno)

# ---------------------------------------------------------
# Collapse annotations for nodes that map to the same plotted node
# ---------------------------------------------------------
clone_map_collapsed <- clone_map %>%
  group_by(newick_label) %>%
  summarise(
    clone_size = sum(clone_size, na.rm = TRUE),
    clone_id = paste(clone_id, collapse = "/"),
    clone_color = if (dplyr::n_distinct(clone_color) == 1) dplyr::first(clone_color) else "grey50",
    gv_nodes = paste(unique(gv_node), collapse = "/"),
    mut_label = paste(unique(mut_label), collapse = "/"),
    gene = paste(unique(gene), collapse = "/"),
    n_attachments = n(),
    .groups = "drop"
  ) %>%
  mutate(
    clone_label = paste0(clone_id, " (", signif(clone_size, 3), ")")
  )

print(clone_map_collapsed)

# =========================================================
# 5. READ NEWICK TREE
# =========================================================

tree <- read.tree(newick_file)

# optional safety check
tree_labels <- tree$tip.label
cat("Tip labels in Newick:", length(tree_labels), "\n")

# =========================================================
# 6. BUILD ANNOTATION TABLE FOR GGTREE
# =========================================================

p0 <- ggtree(tree)

tree_df <- p0$data %>%
  mutate(label = as.character(label))

annot_df <- tree_df %>%
  left_join(
    gv$mutation_nodes %>%
      select(newick_label, mut_label = label, gene),
    by = c("label" = "newick_label")
  ) %>%
  left_join(
    clone_map_collapsed %>%
      select(newick_label, clone_id, clone_size, clone_color, clone_label, n_attachments),
    by = c("label" = "newick_label")
  )

# =========================================================
# 7. PLOT
# =========================================================
# Collapse unitary points to only plot topology
#tree2 <- collapse.singles(tree)
tree$edge.length <- rep(1, nrow(tree$edge))
p <- ggtree(tree) %<+% annot_df +
  geom_tree(linewidth = 0.5, color = "black") +
  
  # mutation/gene labels on internal nodes
  # geom_text2(
  #   aes(
  #     subset = !isTip & !is.na(gene) & gene != "Root",
  #     label = gene
  #   ),
  #   hjust = -0.1,
  #   vjust = 1.1,
  #   size = 2.7,
  #   color = "black"
  # ) +

  # clone-bearing nodes
  geom_point2(
    aes(
      subset = !is.na(clone_id),
      size = clone_size,
      fill = clone_color
    ),
    shape = 21,
    color = "black",
    stroke = 0.4,
    show.legend = FALSE
  ) +
  
  # clone labels
  # geom_text2(
  #   aes(
  #     subset = !is.na(clone_id),
  #     label = clone_id,
  #     color = clone_color
  #   ),
  #   hjust = -0.2,
  #   vjust = -0.4,
  #   size = 3,
  #   show.legend = FALSE
  # ) +
  geom_rootedge(rootedge = 3)+
  theme_tree2()+
  scale_fill_identity() +
  scale_color_identity() +
  scale_size_area(max_size = 70*max(clone_size))+
  scale_x_continuous(limits = c(0,16), expand = c(0.1))

print(p)
