rule select_frequency_variants:
    input:
        vcf = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/recent_snps.vcf"
    output:
        singleton_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/singleton.vcf",
        singleton_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/singleton_table.tsv",
        doubleton_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/doubleton.vcf",
        doubleton_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/doubleton_table.tsv",
        tripleton_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/tripleton.vcf",
        tripleton_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/tripleton_table.tsv"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/frequency_selection.py"