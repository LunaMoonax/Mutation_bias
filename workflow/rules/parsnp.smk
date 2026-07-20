rule parsnp:
    input:
        ref = "results/sourmash/ref/{sp}_ref.fna",
        genome_dir = f"{GENOMES}/{{sp}}/fasta_files"
    output:
        align = temp("results/parsnp/{sp}/parsnp.xmfa"),
        vcf = temp("results/parsnp/{sp}/parsnp.vcf")
    params:
        output_dir = "results/parsnp/{sp}"
    threads: config["parsnp"]["threads"]
    conda:
        "../envs/parsnp.yaml"
    shell:
        """
        parsnp -p {threads} -d {input.genome_dir} \
            -r {input.ref} --skip-phylogeny -c --fo -o {params.output_dir} --vcf
        """