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

rule annotate_frequency_tables:
    input:
        singleton_ann_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/snpeff/singleton.snpeff_annotated.vcf",
        singleton_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/singleton_table.tsv",
        doubleton_ann_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/snpeff/doubleton.snpeff_annotated.vcf",
        doubleton_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/doubleton_table.tsv",
        tripleton_ann_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/snpeff/tripleton.snpeff_annotated.vcf",
        tripleton_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/tripleton_table.tsv",
        ref = "results/sourmash/ref/{sp}_ref.fna"
    output:
        singleton_ann_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/singleton_table.ann.tsv",
        doubleton_ann_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/doubleton_table.ann.tsv",
        tripleton_ann_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/tripleton_table.ann.tsv"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/annotate_frequency_tables.py"