rule snpeff_database:
    input:
        gbff = "results/sourmash/ref/{sp}_ref.gbff"
    output:
        db = os.path.join(config["reference_genomes"]["db_dir"], "snpeff/data/{sp}/snpEffectPredictor.bin")
    params:
        dir = os.path.join(config["reference_genomes"]["db_dir"], "snpeff/data"),
        cfg  = os.path.join(config["reference_genomes"]["db_dir"], "snpeff/{sp}_snpeff.config")
    conda:
        "../envs/snpeff.yaml"
    shell:
        """
        mkdir -p {params.dir}/{wildcards.sp}
        cp {input.gbff} {params.dir}/{wildcards.sp}/genes.gbk

        echo "{wildcards.sp}.genome : {wildcards.sp}" > {params.cfg}

        snpEff build -genbank -dataDir $(realpath {params.dir}) \
            -c {params.cfg} -v {wildcards.sp}
        """

rule run_snpeff:
    input:
        vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq_class}.vcf",
        db = os.path.join(config["reference_genomes"]["db_dir"], "snpeff/data/{sp}/snpEffectPredictor.bin")
    output:
        vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/snpeff/{freq_class}.snpeff_annotated.vcf",
        summary = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/snpeff/{freq_class}_snpeff_summary.html"
    wildcard_constraints:
        freq_class = "singleton|doubleton|tripleton"
    params:
        dir = os.path.join(config["reference_genomes"]["db_dir"], "snpeff/data"),
        cfg = os.path.join(config["reference_genomes"]["db_dir"], "snpeff/{sp}_snpeff.config")
    conda:
        "../envs/snpeff.yaml"
    shell:
        """
        snpEff ann -dataDir $(realpath {params.dir}) \
            -c {params.cfg} -stats {output.summary} \
            {wildcards.sp} {input.vcf} > {output.vcf}
        """