#!/usr/bin/env Rscript

library(treeio)
library(ggtree)
library(ggplot2)

tree_file <- snakemake@input[["tree"]]
qc_file <- snakemake@input[["qc"]]
species <- snakemake@wildcards[["sp"]]
output_file <- snakemake@output[["tree_colored"]]

clean_id <- function(x) {
    x <- sub("\\.ref$", "", x)
    sub("\\.fa$", "", x)
}

tree <- read.tree(tree_file)
qc <- read.csv(qc_file, stringsAsFactors = FALSE)
failed_species <- qc[qc$species == species & qc$exclude == TRUE, , drop = FALSE]
outliers_id <- failed_species$genome

data <- data.frame(
    label = tree$tip.label,
    id = clean_id(tree$tip.label),
    status = ifelse(clean_id(tree$tip.label) %in% outliers_id, "Outlier", "Kept"),
    stringsAsFactors = FALSE
)

n_tips <- nrow(data)
n_out <- sum(data$status == "Outlier")

p <- ggtree(tree, size = 0.2) %<+% data +
    geom_tippoint(aes(color = status)) +
    scale_color_manual(values = c(Kept = "grey70", Outlier = "red")) +
    theme_tree2() +
    xlab("Substitutions per site") +
    labs(color = NULL, title = sprintf("%s  (%d outliers / %d genomes)",
                                       species, n_out, n_tips)) +
    theme(legend.position = "right")

ggsave(output_file, p, width = 10, height = max(8, n_tips * 0.015),
       dpi = 200, limitsize = FALSE)
