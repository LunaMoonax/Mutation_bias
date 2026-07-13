rule select_snp_tree:
    input:
        tree = "results/dataprep/{sp}/recombination/{sp}.node_labelled.final_tree.tre",
        vcf = "results/dataprep/{sp}/recombination_free.vcf",
        embl = "results/dataprep/{sp}/recombination/{sp}.branch_base_reconstruction.ref_coords.embl"
    output:
        vcf_filtered = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/recent_snps.vcf",
        branch_table = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/branch_table.txt",
        stats = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/stats.txt",
        zero_lookahead = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/zero_lookahead.txt",
        tree_fig = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/tree_fig.png"
    params:
        snp_threshold = lambda wildcards: int(wildcards.threshold)
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/select_snp_from_tree.R"