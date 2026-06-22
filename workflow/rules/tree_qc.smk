rule pick_tree_ref:
    input:
        qc = "results/ani/ani_qc.csv"
    output:
        ref = "results/tree/{sp}_ref.txt"
    params:
        genomes_dir = GENOMES
    script:
        "../scripts/pick_tree_ref.py"

rule parsnp_allign:
    input:
        ref = "results/tree/{sp}_ref.txt",
        genome_dir = f"{GENOMES}/{{sp}}/fasta_files"
    output:
        align = "results/tree/parsnp/{sp}_parsnp/parsnp.snps.mblocks"
    params:
        output_dir = directory("results/tree/parsnp/{sp}_parsnp")
    threads: config["parsnp"]["threads"]
    conda:
        "../envs/parsnp.yaml"
    shell:
        """
        parsnp -p {threads} -d {input.genome_dir} \
            -r $(cat {input.ref}) --skip-phylogeny -c --fo -o {params.output_dir}
        """

rule fasttree:
    input:
        align = "results/tree/parsnp/{sp}_parsnp/parsnp.snps.mblocks"
    output:
        tree = "results/tree/{sp}.nwk"
    conda:
        "../envs/fasttree.yaml"
    shell:
        """
        fasttree -nt -gtr -gamma {input.align} > {output.tree}
        """
    
    