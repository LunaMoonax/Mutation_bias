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
        
    
    