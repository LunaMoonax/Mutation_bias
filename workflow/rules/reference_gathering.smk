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
    