rule reference_database:
    output:
        db = os.path.join(config["reference_genomes"]["db_dir"], "refseq_db.tar")
    params:
        db_dir = config["reference_genomes"]["db_dir"]
    shell:
        """
        curl -fL https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k31.tar.gz \
            -o {params.db_dir}/refseq_db.tar.gz
        gunzip {params.db_dir}/refseq_db.tar.gz
        """
    