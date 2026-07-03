rule filter_vcf:
    input:
        vcf = "results/parsnp/{sp}/parsnp.vcf"
    output:
        vcf = "results/dataprep/{sp}/filtered.vcf"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/filter_vcf.py"

rule gubbins:
    input:
        xmfa = "results/parsnp/{sp}/parsnp.xmfa"
    output:
        aln  = "results/dataprep/{sp}/recombination/alignment.fasta",
        recomb = "results/dataprep/{sp}/recombination/{sp}.recombination_predictions.gff"
    params:
        prefix = "results/dataprep/{sp}/recombination/gubbins"
    threads: config["gubbins"]["threads"]
    conda:
        "../envs/gubbins.yaml"
    shell:
        """
        harvesttools -x {input.xmfa} -M {output.aln}

        aln=$(realpath {output.aln})
        outdir=$(realpath -m results/dataprep/{wildcards.sp}/recombination)

        workdir=$(mktemp -d /tmp/gubbins_{wildcards.sp}_XXXXXX)
        trap 'rm -rf "$workdir"' EXIT
        cd "$workdir"

        run_gubbins.py --mar --prefix {wildcards.sp} --threads {threads} \
            --first-model JC --tree-builder fasttree --iterations 3 "$aln"

        cp -a {wildcards.sp}.* "$outdir"/
        """