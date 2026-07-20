rule select_frequency_variants:
    input:
        vcf = "results/dataprep/{sp}/tree_snp/snp_per_branch_{threshold}/recent_snps.vcf"
    output:
        expand("results/dataprep/{{sp}}/frequency_spectrum/snp_per_branch_{{threshold}}/{freq}.vcf",
               freq=["singleton", "doubleton", "tripleton"]),
        expand("results/dataprep/{{sp}}/frequency_spectrum/snp_per_branch_{{threshold}}/{freq}_table.tsv",
               freq=["singleton", "doubleton", "tripleton"]),
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/frequency_selection.py"

rule annotate_frequency_table:
    input:
        ann_vcf = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/snpeff/{freq}.snpeff_annotated.vcf",
        table   = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq}_table.tsv",
        ref     = "results/sourmash/ref/{sp}_ref.fna"
    output:
        ann_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq}_table.ann.tsv"
    wildcard_constraints:
        freq = "singleton|doubleton|tripleton"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/annotate_frequency_tables.py"

rule generate_substitutons:
    input:
        ref = "results/sourmash/ref/{sp}_ref.fna"
    output:
        vcf = "results/dataprep/{sp}/counts/all_substitutions.vcf"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/generate_all_substitutions.py"

rule count_opportunities:
    input:
        ann_vcf = "results/dataprep/{sp}/counts/all_substitutions.ann.vcf",
        ref = "results/sourmash/ref/{sp}_ref.fna"
    output:
        table = "results/dataprep/{sp}/counts/count_table.tsv"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/count_opportunities.py"

rule build_canonical:
    input:
        opp_table = "results/dataprep/{sp}/counts/count_table.tsv",
        table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq}_table.ann.tsv"
    output:
        raw_cann_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq}_raw_canonical.tsv"
    wildcard_constraints:
        freq = "singleton|doubleton|tripleton"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/build_canonical_table.py"

rule build_species_symmetry:
    input:
        tables = expand(
            "results/dataprep/{{sp}}/frequency_spectrum/snp_per_branch_" + str(max(THRESHOLDS)) + "/{freq}_raw_canonical.tsv",
            freq=["singleton", "doubleton", "tripleton"])
    output:
        sym_results = "results/dataprep/{sp}/counts/species_strand_symmetry.tsv",
        decision    = "results/dataprep/{sp}/counts/species_collapse_decision.txt"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/build_species_symmetry.py"

rule collapse_normalise:
    input:
        raw_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq}_raw_canonical.tsv",
        decision  = "results/dataprep/{sp}/counts/species_collapse_decision.txt"
    output:
        final_table = "results/dataprep/{sp}/frequency_spectrum/snp_per_branch_{threshold}/{freq}_canonical.tsv"
    wildcard_constraints:
        freq = "singleton|doubleton|tripleton"
    conda:
        "../envs/py_dataprep.yaml"
    script:
        "../scripts/build_collapse.py"