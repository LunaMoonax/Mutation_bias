#!/usr/bin/env python3

import pandas as pd

qc_file      = snakemake.input.qc
genomes_dir  = snakemake.params.genomes_dir
species      = snakemake.wildcards.sp
out_file     = snakemake.output.ref


qc = pd.read_csv(qc_file)

passed = qc[(qc["species"] == species) & (qc["exclude"] == False)]

species_median = passed["species_median"].iloc[0]
ref_idx = (passed["median_ani"] - species_median).abs().idxmin()
ref_genome = passed.loc[ref_idx, "genome"]
ref_path = f"{genomes_dir}/{species}/fasta_files/{ref_genome}.fa"

with open(out_file, "w") as f:
    f.write(f"{ref_path}\n")