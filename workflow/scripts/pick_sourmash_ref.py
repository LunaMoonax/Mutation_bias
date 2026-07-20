#!/usr/bin/env python3

# pick the NCBI reference genome for a species: closest sourmash hit that is a
# complete, low-contig-count assembly; log everything rejected as a fallback.
# If nothing clears the primary similarity bar, retry once at a relaxed
# fallback bar and record that the bar was lowered. If nothing passes even the
# fallback, the species is excluded -- recorded in {sp}_excluded.txt, written
# directly (not a declared Snakemake output) so it survives even though the
# job itself fails and Snakemake cleans up the rule's declared outputs.

import os
import pandas as pd

sourmash_csv = snakemake.input.sourmash_csv
ncbi_tsv = snakemake.input.ncbi_tsv
fallbacks_file = snakemake.output.fallbacks
ref_file = snakemake.output.ref
similarity = snakemake.params.sim
fallback_similarity = snakemake.params.sim_fallback
contig_n = snakemake.params.contig
species = snakemake.wildcards.sp

ref_dir = os.path.dirname(ref_file)
threshold_note_file = os.path.join(ref_dir, f"{species}_threshold_note.txt")
excluded_file = os.path.join(ref_dir, f"{species}_excluded.txt")

sourmash = pd.read_csv(sourmash_csv)
ncbi = pd.read_csv(ncbi_tsv, sep="\t")

sourmash["accession"] = sourmash["name"].str.split().str[0]

combined = sourmash.merge(ncbi, left_on="accession",
                          right_on="Assembly Accession", how="left")
have_meta = combined[combined["Assembly Accession"].notna()]

# walk sourmash hits in order (best match first) at the given similarity bar;
# return the first that qualifies (Complete Genome, <= contig_n contigs), and
# append every rejected candidate (with its reason) to `fallbacks`
def find_reference(min_similarity, fallbacks):
    for _, row in have_meta.iterrows():
        if row["similarity"] >= min_similarity:
            if row["Assembly Level"] == "Complete Genome" and row["Assembly Stats Number of Contigs"] <= contig_n:
                return row["accession"]
            fallbacks.append({"species": species, "accession": row["accession"],
                               "contig_number": row["Assembly Stats Number of Contigs"],
                               "assembly_level": row["Assembly Level"],
                               "reason": "Not complete genome or big contig number"})
        else:
            fallbacks.append({"species": species, "accession": row["accession"],
                               "similarity": row["similarity"],
                               "reason": "Not enough similarity"})
    return None

fallbacks = []
chosen_ref = find_reference(similarity, fallbacks)
threshold_used = similarity
lowered = False

# nothing cleared the primary bar -- retry once at the relaxed fallback bar.
# start `fallbacks` fresh: the fallback pass re-evaluates every candidate, so
# reusing the primary pass's list would duplicate every candidate that
# appears in both passes
if chosen_ref is None:
    fallbacks = []
    chosen_ref = find_reference(fallback_similarity, fallbacks)
    if chosen_ref is not None:
        threshold_used = fallback_similarity
        lowered = True

pd.DataFrame(fallbacks).to_csv(fallbacks_file, index=False)

if chosen_ref is None:
    with open(excluded_file, "w") as f:
        f.write(f"species\t{species}\n")
        f.write("excluded\tTrue\n")
        f.write(f"reason\tNo reference passed QC at similarity>={similarity} or fallback>={fallback_similarity}\n")
    raise ValueError(
        f"No reference passed QC for {species} "
        f"(tried similarity>={similarity} and fallback>={fallback_similarity})"
    )

with open(ref_file, "w") as f:
    f.write(chosen_ref + "\n")

with open(threshold_note_file, "w") as f:
    f.write(f"species\t{species}\n")
    f.write(f"threshold_used\t{threshold_used}\n")
    f.write(f"threshold_lowered\t{lowered}\n")
