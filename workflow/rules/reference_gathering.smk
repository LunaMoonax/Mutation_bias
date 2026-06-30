rule reference_database:
    output:
        db = os.path.join(config["reference_genomes"]["db_dir"], "gtdb-rs226-k31.dna.zip")
    shell:
        """
        curl -fL \
            https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-rs226-k31.dna.zip \
            -o {output.db}
        """

rule pick_genome:
    input:
        qc = "results/QC/ani/ani_qc.csv"
    output:
        ref = "results/sourmash/genome/{sp}.txt"
    params:
        genomes_dir = GENOMES
    script:
        "../scripts/pick_ref.py"

rule sketching_sample:
    input:
        fasta_file = "results/sourmash/genome/{sp}.txt"
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
        ntop = config["ncbi"]["data_top"]
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

rule filter_pick_ref:
    input:
        sourmash_csv = "results/sourmash/{sp}_results.csv",
        ncbi_tsv = "results/sourmash/ncbi/{sp}_data.tsv"
    output:
        fallbacks = "results/sourmash/ref/{sp}_fallbacks.csv",
        ref = "results/sourmash/ref/{sp}_ref.acc"
    params:
        sim = config["sourmash"]["similarity"],
        contig = config["ncbi"]["contig_n"]
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/pick_sourmash_ref.py"

rule download_ref:
    input:
        ref = "results/sourmash/ref/{sp}_ref.acc"
    output:
        fasta = "results/sourmash/ref/{sp}_ref.fna",
        gff =   "results/sourmash/ref/{sp}_ref.gff",
        gbff  = "results/sourmash/ref/{sp}_ref.gbff"
    params:
        zip = "results/sourmash/ref/{sp}.zip",
        tmp = "results/sourmash/ref/{sp}_tmp"
    conda:
        "../envs/ncbi.yaml"
    shell:
        """
        datasets download genome accession $(cat {input.ref}) \
            --include genome,gff3,gbff \
            --filename {params.zip}

        unzip -o {params.zip} -d {params.tmp}

        cp {params.tmp}/ncbi_dataset/data/*/*.fna {output.fasta}
        cp {params.tmp}/ncbi_dataset/data/*/*.gff {output.gff}
        cp {params.tmp}/ncbi_dataset/data/*/*.gbff {output.gbff}

        rm -rf {params.tmp} {params.zip}
        """
    