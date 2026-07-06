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

rule remap_gubbins_coords:
    input:
        xmfa = "results/parsnp/{sp}/parsnp.xmfa",
        gff = "results/dataprep/{sp}/recombination/{sp}.recombination_predictions.gff"
    output:
        gff = "results/dataprep/{sp}/recombination/{sp}.recombination_predictions.ref_coords.gff"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/remap_gubbins_coords.py"

rule sliding_window:
    input:
        vcf = "results/dataprep/{sp}/filtered.vcf"
    output:
        flagged = "results/dataprep/{sp}/recombination/window_flagged.bed"
    params:
        window = config["sliding_window"]["window"],
        threshold = config["sliding_window"]["threshold"]
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/rec_sliding_window.py"

rule compare:
    input:
        window_bed = "results/dataprep/{sp}/recombination/window_flagged.bed",
        gubbins_gff = "results/dataprep/{sp}/recombination/{sp}.recombination_predictions.ref_coords.gff",
        vcf = "results/dataprep/{sp}/filtered.vcf"
    output:
        summary = "results/dataprep/{sp}/recombination/comparing_summary.txt"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/window_vs_gubbins.py"

rule remove_recombination:
    input:
        gff = "results/dataprep/{sp}/recombination/{sp}.recombination_predictions.ref_coords.gff",
        vcf = "results/dataprep/{sp}/filtered.vcf"
    output:
        rec_free_vcf = "results/dataprep/{sp}/recombination_free.vcf"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/remove_recombinant_snps.py"
