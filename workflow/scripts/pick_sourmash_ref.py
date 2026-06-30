#!/usr/bin/env python3

import pandas as pd

sourmash_csv = snakemake.input.sourmash_csv
ncbi_tsv = snakemake.input.ncbi_tsv
fallbacks_file = snakemake.output.fallbacks
ref_file = snakemake.output.ref
similarity = snakemake.params.sim
contig_n = snakemake.params.contig
species = snakemake.wildcards.sp

sourmash = pd.read_csv(sourmash_csv)
ncbi = pd.read_csv(snakemake.input.ncbi_tsv, sep="\t")

sourmash["accession"] = sourmash["name"].str.split().str[0]

combined = sourmash.merge(ncbi, left_on="accession",
                          right_on="Assembly Accession", how="left")
have_meta = combined[combined["Assembly Accession"].notna()]

fallbacks = []
chosen_ref = None

for _, row in have_meta.iterrows():
    if row["similarity"] >= similarity:
        if row["Assembly Level"] == "Complete Genome" and row["Assembly Stats Number of Contigs"] <= contig_n:
            chosen_ref = row["accession"]
            break
        else:
            fallbacks.append({"species": species, "accession": row["accession"],
                                           "contig_number": row["Assembly Stats Number of Contigs"],
                                           "assembly_level": row["Assembly Level"],
                                           "reason": "Not complete genome or big contig number"})
    else:
        fallbacks.append({"species": species, "accession": row["accession"],
                                           "similarity": row["similarity"],
                                           "reason": "Not enough similarity"})

pd.DataFrame(fallbacks).to_csv(fallbacks_file, index=False)

with open(ref_file, "w") as f:
    if chosen_ref is not None:
        f.write(chosen_ref + "\n")
    else:
        raise ValueError(f"No reference passed QC for {species}")