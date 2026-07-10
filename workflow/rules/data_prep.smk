rule select_singletons:
    input:
        vcf = "results/dataprep/{sp}/tree_snp/recent_snps.vcf"
    output:
        singleton_vcf = "results/dataprep/{sp}/singleton.vcf",
        singleton_table = "results/dataprep/{sp}/singleton_table.tsv"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/singleton_selection.py"