rule reference_database:
    output:
        db = os.path.join(config["reference_genomes"]["db_dir"], "gtdb-rs226-k31.dna.zip")
    shell:
        """
        curl -fL \
            https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-rs226-k31.dna.zip \
            -o {output.db}
        """

rule pick_ref:
    input:
        qc = "results/QC/ani/ani_qc.csv"
    output:
        ref = "results/sourmash/ref/{sp}_ref.txt"
    params:
        genomes_dir = GENOMES
    script:
        "../scripts/pick_ref.py"

rule sketching_sample:
    input:
        fasta_file = "results/sourmash/ref/{sp}_ref.txt"
    output:
        sig = "results/sourmash/sig_files/{sp}_genome.sig"
    params:
        k = config["sourmash"]["k"],
        scale = config["sourmash"]["scale"]
    conda:
        "../envs/sourmash.yaml"
    shell:
        """
        sourmash sketch dna -p scaled={params.scale},k={params.k} \
            $(cat {input.fasta_file}) -o {output.sig} 
        """

rule sourmash_search:
    input:
        db = os.path.join(config["reference_genomes"]["db_dir"], "gtdb-rs226-k31.dna.zip"),
        sig = "results/sourmash/sig_files/{sp}_genome.sig"
    output:
        res = "results/sourmash/{sp}_results.csv"
    conda:
        "../envs/sourmash.yaml"
    shell:
        """
        sourmash search {input.sig} {input.db} --containment -o {output.res}
        """

rule ncbi_datasets:
    input:
        csv = "results/sourmash/{sp}_results.csv"
    output:
        data = "results/sourmash/ncbi/{sp}_data.tsv"
    params:
        ntop = config["sourmash"]["data_top"]
    conda:
        "../envs/ncbi.yaml"
    shell:
        """
        awk -F',' -v n={params.ntop} 'NR>1 && NR<=n+1 {{print $4}}' {input.csv} \
            | awk '{{print $1}}' > {output.data}.acc

        datasets summary genome accession --inputfile {output.data}.acc --as-json-lines \
            | dataformat tsv genome \
                --fields accession,organism-name,assminfo-level,assmstats-number-of-contigs \
            > {output.data}
        """
    