rule plot_close_clusters:
    input:
        tree = "results/QC/tree/{sp}.nwk",
        ani_qc = "results/QC/ani/ani_qc.csv"
    output: 
        plot = "results/QC/close_clusters/{sp}_distribution.png"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_close_clusters.R"