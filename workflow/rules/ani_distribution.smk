rule skani:
    input:
        fastas = lambda wc: glob.glob(f"{GENOMES}/{wc.sp}/fasta_files/*.fa")
    output:
        ani = "results/QC/skani/{sp}_ani_edge_list.txt"
    params:
        listfile = lambda wc: f"results/skani/{wc.sp}_list.txt"
    conda:
        "../envs/skani.yaml"
    threads: config["skani"]["threads"]
    shell:
        """
        printf '%s\n' {input.fastas} > {params.listfile}
        skani triangle -l {params.listfile} -E -t {threads} -o {output.ani}
        """
        
rule ani_distributions:
    input:
        edge_list   = expand("results/QC/skani/{sp}_ani_edge_list.txt", sp=SPECIES),
        genome_list = expand("results/QC/skani/{sp}_list.txt", sp=SPECIES),
    output:
        qc = "results/QC/ani/ani_qc.csv",
        ani_dist = "results/QC/ani/ani_distribution.png",
        ani_violin = "results/QC/ani/ani_violin.png"
    params:
        output_dir = "results/QC/ani/",
        floor = config["ani_distribution_plot"]["floor"],
        mad_k = config["ani_distribution_plot"]["mad_k"],
        min_drop = config["ani_distribution_plot"]["min_drop"]
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/ani_distribution_plot.py"
    