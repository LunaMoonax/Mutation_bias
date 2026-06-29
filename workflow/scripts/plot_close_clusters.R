#!/usr/bin/env Rscript

library(ape)
library(ggplot2)

tree_file <- snakemake@input[["tree"]]
ani_qc_file <- snakemake@input[["ani_qc"]]
species <- snakemake@wildcards[["sp"]]
plot_output <- snakemake@output[["plot"]]

clean_id <- function(x) {
    x <- sub("\\.ref$", "", x)
    sub("\\.fa$", "", x)
}

tree <- read.tree(tree_file)
n_tips <- length(tree$tip.label)
is_terminal <- tree$edge[, 2] <= n_tips
terminal_length <- tree$edge.length[is_terminal]
terminal_tip <- clean_id(tree$tip.label[tree$edge[is_terminal, 2]])

ani_qc <- read.csv(ani_qc_file, stringsAsFactors = FALSE)
species_qc <- ani_qc[ani_qc$species == species, , drop = FALSE]

df <- data.frame(
    genome = terminal_tip,
    terminal_branch = terminal_length,
    stringsAsFactors = FALSE
)

df <- merge (df, species_qc[, c("genome", "median_ani")],
             by = "genome", all.x = TRUE)

p <- ggplot(df, aes(x = median_ani, y = terminal_branch)) +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Median ANI", y = "Terminal branch length (subs/site)",
         title = sprintf("%s  (%d genomes)", species, nrow(df))) +
    theme_minimal()

ggsave(plot_output, p, width = 8, height = 6, dpi = 200)
