#!/usr/bin/env Rscript

# select "recent" SNPs: those reconstructed onto a terminal branch with 1..N events

library(ape)

tree_file <- snakemake@input[["tree"]]
vcf_file <- snakemake@input[["vcf"]]
embl_file <- snakemake@input[["embl"]]
vcf_out_file <- snakemake@output[["vcf_filtered"]]
table_file <- snakemake@output[["branch_table"]]
stats_file <- snakemake@output[["stats"]]
zero_table_file <- snakemake@output[["zero_lookahead"]]
tree_fig_file <- snakemake@output[["tree_fig"]]
N <- snakemake@params[["snp_threshold"]]

# parse gubbins' ancestral-reconstruction EMBL into one row per mutation event (pos + branch)
read_embl <- function(embl_file) {
    lines <- readLines(embl_file)

    is_var <- grepl("^FT\\s+variation\\s+\\d+\\s*$", lines)
    var_idx <- which(is_var)
    if (length(var_idx) == 0) stop("no variation records found in ", embl_file)

    pos <- as.integer(sub("^FT\\s+variation\\s+(\\d+)\\s*$", "\\1", lines[var_idx]))

    ends <- c(var_idx[-1] - 1L, length(lines))

    node <- character(length(var_idx))
    for (i in seq_along(var_idx)) {
        block <- lines[var_idx[i]:ends[i]]
        nl <- grep('/node="', block, value = TRUE)[1]
        node[i] <- if (is.na(nl)) NA_character_ else sub('.*/node="([^"]+)".*', "\\1", nl)
    }

    data.frame(pos = pos, node = node, stringsAsFactors = FALSE)
}

embl <- read_embl(embl_file)
n_events_all <- nrow(embl)
n_pos_all    <- length(unique(embl$pos))

# keep only events whose position survived recombination filtering (is still in the VCF)
vcf_lines <- readLines(vcf_file)
is_header <- startsWith(vcf_lines, "#")
vcf_pos <- as.integer(sub("^[^\t]+\t(\\d+)\t.*", "\\1", vcf_lines[!is_header]))

embl <- embl[embl$pos %in% vcf_pos, ]
n_events_drop   <- n_events_all - nrow(embl)
n_events_usable <- nrow(embl)
n_pos_usable    <- length(unique(embl$pos))

# split "parent->child" node label; internal nodes are named "internal_*" by gubbins
embl$child <- sub(".*->", "", embl$node)
embl$terminal <- !startsWith(embl$child, "internal_")
embl$parent <- sub("->.*", "", embl$node)

# count SNP events per branch
counts <- aggregate(pos ~ node + child + terminal + parent, data = embl, FUN = length)
names(counts)[names(counts) == "pos"] <- "n_snps"

tree <- read.tree(tree_file)
tips <- tree$tip.label

n_tips_tree <- Ntip(tree)

# look up a tip's immediate parent's node label via the tree's edge matrix
parent_of_tip <- function(tip_names) {
  tip_no    <- match(tip_names, tree$tip.label)
  edge_i    <- match(tip_no, tree$edge[, 2])
  parent_no <- tree$edge[edge_i, 1]
  tree$node.label[parent_no - n_tips_tree]
}

# tips with zero recorded mutation events are missing from `counts` entirely; add them back
seen         <- counts$child[counts$terminal]
zeros        <- setdiff(tips, seen)
zero_parents <- if (length(zeros) > 0) parent_of_tip(zeros) else character(0)

if (length(zeros) > 0) {
  counts <- rbind(counts, data.frame(
    node     = paste0(zero_parents, "->", zeros),
    child    = zeros,
    terminal = TRUE,
    parent   = zero_parents,
    n_snps   = 0L,
    stringsAsFactors = FALSE
  ))
}

counts <- counts[order(-counts$terminal, counts$n_snps), ]
write.table(counts, table_file, sep = "\t", quote = FALSE, row.names = FALSE)

# the actual selection: terminal branches with 1..N events are "recent"
short_nodes <- counts$node[counts$terminal & counts$n_snps >= 1 & counts$n_snps <= N]

recent_events <- embl[embl$node %in% short_nodes, ]
recent_pos <- unique(recent_events$pos)

# write out the filtered VCF (one line per position, regardless of how many
# branches independently reconstructed a mutation there)
keep <- is_header
keep[!is_header] <- vcf_pos %in% recent_pos
writeLines(vcf_lines[keep], vcf_out_file)
n_written <- sum(keep & !is_header)

snps_into <- setNames(counts$n_snps, counts$child)

# diagnostic table: for each 0-event tip, how many events did its parent branch have
zero_lookahead <- data.frame(tip = character(0), parent_label = character(0),
                             parent_events = integer(0), stringsAsFactors = FALSE)

if (length(zeros) > 0) {
  pev <- snps_into[zero_parents]
  pev[is.na(pev)] <- 0L

  zero_lookahead <- data.frame(
    tip           = zeros,
    parent_label  = zero_parents,
    parent_events = as.integer(pev),
    stringsAsFactors = FALSE
  )
}

write.table(zero_lookahead, zero_table_file, sep = "\t",
            quote = FALSE, row.names = FALSE)

term <- counts[counts$terminal, ]

n_pos_terminal <- length(unique(embl$pos[embl$terminal]))
n_pos_internal <- length(unique(embl$pos[!embl$terminal]))

max_show <- 15L

# table A: terminal branches, by exact event count k -- what the output VCF
# would contain if the threshold were set to k
terminal_dist_at_k <- function(k) {
  nodes_upto <- counts$node[counts$terminal & counts$n_snps >= 1 & counts$n_snps <= k]
  list(branches_at   = sum(term$n_snps == k),
       cum_branches  = sum(term$n_snps >= 1 & term$n_snps <= k),
       cum_uniq_snps = length(unique(embl$pos[embl$node %in% nodes_upto])))
}

# table B: 0-event tips, by their parent branch's exact event count k -- what
# recovering those tips (using the parent branch's SNPs instead) would add
zero_dist_at_k <- function(k) {
  uniq_parents <- unique(zero_lookahead$parent_label[zero_lookahead$parent_events >= 1 &
                                                       zero_lookahead$parent_events <= k])
  parent_nodes <- counts$node[counts$child %in% uniq_parents]
  list(tips_at       = sum(zero_lookahead$parent_events == k),
       cum_tips      = sum(zero_lookahead$parent_events >= 1 & zero_lookahead$parent_events <= k),
       cum_uniq_snps = length(unique(embl$pos[embl$node %in% parent_nodes])))
}

# write the stats report (both to stdout and to file) summarizing the threshold's effect
con <- file(stats_file, "w")
out <- function(s) { cat(s, "\n", sep = ""); writeLines(s, con) }

out("")
out(sprintf("species                          : %s", snakemake@wildcards[["sp"]]))
out(sprintf("threshold N (current)            : %d", N))
out("")

out("=== mutation events ===")
out(sprintf("total (raw, EMBL)                : %d events, %d unique positions", n_events_all, n_pos_all))
out(sprintf("usable (post recomb/filtering)   : %d events, %d unique positions", n_events_usable, n_pos_usable))
out(sprintf("  dropped, not in filtered VCF   : %d events", n_events_drop))
out("")

out("=== terminal branches ===")

out(sprintf("tree tips                        : %d", n_tips_tree))
out(sprintf("  with 0 events                  : %d", sum(term$n_snps == 0)))
out(sprintf("  with >=1 event                 : %d", sum(term$n_snps >= 1)))
out(sprintf("unique positions on terminal     : %d", n_pos_terminal))
out(sprintf("unique positions on internal     : %d", n_pos_internal))
out("")

out(sprintf("=== selected for downstream, at current N = %d ===", N))
out(sprintf("terminal tips selected (1..N events)   : %d", terminal_dist_at_k(N)$cum_branches))
out(sprintf("unique positions written to VCF        : %d", n_written))
out("")

out("=== Terminal branches -- SNP-count distribution ===")
out("")
out(sprintf("  %3s | %10s %12s %12s", "snps", "branches", "sum branches", "sum unique snps"))
for (k in 0:max_show) {
  r <- terminal_dist_at_k(k)
  out(sprintf("  %3d | %10d %12d %12d", k, r$branches_at, r$cum_branches, r$cum_uniq_snps))
}
n_over <- sum(term$n_snps > max_show)
if (n_over > 0) out(sprintf("  %3s | %10d %12s %12s", paste0(">", max_show), n_over, "-", "-"))
out("")

out("=== 0-event tips -- recovery via parent branch (test) ===")
out("")
out(sprintf("0-event terminal tips                  : %d", nrow(zero_lookahead)))
out(sprintf("  parent branch also has 0 events : %d",
            sum(zero_lookahead$parent_events == 0)))
out("")
out(sprintf("  %3s | %10s %12s %12s", "snps", "tips", "sum tips", "sum unique snps"))
for (k in 0:max_show) {
  r <- zero_dist_at_k(k)
  out(sprintf("  %3d | %10d %12d %12d", k, r$tips_at, r$cum_tips, r$cum_uniq_snps))
}
n_over_z <- sum(zero_lookahead$parent_events > max_show)
if (n_over_z > 0) out(sprintf("  %3s | %10d %12s %12s", paste0(">", max_show), n_over_z, "-", "-"))

close(con)

# plot the tree, colouring tips by whether they were selected as "recent"
short_tips <- counts$child[counts$node %in% short_nodes]

tip_col <- ifelse(tree$tip.label %in% short_tips, "firebrick",
           ifelse(tree$tip.label %in% zeros,      "grey85", "grey50"))

png(tree_fig_file, width = 2000, height = 2000, res = 150)
plot(tree, show.tip.label = FALSE,
     edge.color = "grey60", edge.width = 0.4, no.margin = TRUE)
tiplabels(pch = 19, col = tip_col, cex = 0.35)
legend("topright",
       legend = c(sprintf("recent (1-%d SNPs)", N),
                  "0 SNPs (no mutations)",
                  sprintf("> %d SNPs", N)),
       col = c("firebrick", "grey85", "grey50"),
       pch = 19, bty = "n", cex = 0.9)
title(sprintf("%s  -  %d of %d tips selected",
              snakemake@wildcards[["sp"]], length(short_tips), length(tree$tip.label)))
dev.off()