rule skani:
    input:
        fastas = lambda wc: glob.glob(f"{GENOMES}/{wc.sp}/fasta_files/*.fa")
    output:
        ani = "results/skani/{sp}_ani_edge_list.txt"
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
        edge_list   = expand("results/skani/{sp}_ani_edge_list.txt", sp=SPECIES),
        genome_list = expand("results/skani/{sp}_list.txt", sp=SPECIES),
    output:
        passed = expand("results/ani/{sp}_passed_genomes.txt", sp=SPECIES),
        failed = "results/ani/failed_genomes.txt",
        ani_dist = "results/ani/ani_distribution.png",
        ani_violin = "results/ani/ani_violin.png"
    params:
        output_dir = "results/ani/"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/ani_distribution_plot.py"
    