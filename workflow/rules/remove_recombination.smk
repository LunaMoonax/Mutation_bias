rule filter_vcf:
    input:
        vcf = "results/parsnp/{sp}/parsnp.vcf"
    output:
        vcf = "results/dataprep/{sp}/filtered.vcf"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/filter_vcf.py"