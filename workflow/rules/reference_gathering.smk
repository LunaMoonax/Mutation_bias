rule reference_database:
    output:
        db = os.path.join(config["reference_genomes"]["db_dir"], "gtdb-rs226-k31.dna.zip"),
        tax = os.path.join(config["reference_genomes"]["db_dir"], "gtdb-rs226.lineages.csv")
    shell:
        """
        curl -fL \
            https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-rs226-k31.dna.zip \
            -o {output.db}
        
        curl -fL \
            https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.csv \
            -o {output.tax}
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
    