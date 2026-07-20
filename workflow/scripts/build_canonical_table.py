#!/usr/bin/env python3

# merge observed mutation counts with opportunity counts into one per-cell table

import pandas as pd

table = snakemake.input.table
opp_table = snakemake.input.opp_table

raw_cann_table = snakemake.output.raw_cann_table

KEYS = ["MUTATION_CLASS", "FIVE PRIME", "THREE PRIME", "CODON POS", "EFFECT"]

opp = pd.read_csv(opp_table, sep="\t", comment="#", dtype={k: str for k in KEYS})

# count observed mutations per cell
obs = pd.read_csv(table, sep="\t", comment="#", dtype={k: str for k in KEYS})
num = obs.groupby(KEYS).size().reset_index(name="N_MUTATIONS")

# every observed cell must already exist in the opportunity table (fail loudly if not)
canon = opp.merge(num, how="outer", on=KEYS)
assert not canon["OPPORTUNITIES"].isna().any(), "observed cell with no opportunity — key mismatch"
canon["OPPORTUNITIES"] = canon["OPPORTUNITIES"].astype(int)
canon["N_MUTATIONS"] = canon["N_MUTATIONS"].fillna(0).astype(int)
canon = canon[KEYS + ["N_MUTATIONS", "OPPORTUNITIES"]]
canon.to_csv(raw_cann_table, sep="\t", index=False)