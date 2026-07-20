#!/usr/bin/env python3

# apply the species-level collapse decision: fold 12 mutation classes down to the
# 6 pyrimidine-referenced ones if symmetric, otherwise keep all 12 as-is

import pandas as pd

raw_table_path = snakemake.input.raw_table
decision_path  = snakemake.input.decision
out_path       = snakemake.output.final_table

KEYS = ["MUTATION_CLASS", "FIVE PRIME", "THREE PRIME", "CODON POS", "EFFECT"]

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

PYRIMIDINE_REF = {"C", "T"}

# reverse-complement a mutation class string, e.g. "G>A" -> "C>T"
def revcomp_class(mut):
    ref, alt = mut.split(">")
    return f"{COMPLEMENT[ref]}>{COMPLEMENT[alt]}"

# pull the "decision" value out of the species_collapse_decision.txt key-value file
def read_decision(path):
    with open(path) as f:
        for line in f:
            if line.startswith("decision\t"):
                return line.strip().split("\t")[1]
    raise ValueError(f"no decision line in {path}")

canon = pd.read_csv(raw_table_path, sep="\t", dtype={k: str for k in KEYS})
canon["N_MUTATIONS"] = canon["N_MUTATIONS"].astype(int)
canon["OPPORTUNITIES"] = canon["OPPORTUNITIES"].astype(int)
decision = read_decision(decision_path)

if decision == "collapse_to_6":
    def folded_key(row):
        cls = row["MUTATION_CLASS"]
        if cls.split(">")[0] in PYRIMIDINE_REF:
            return (cls, row["FIVE PRIME"], row["THREE PRIME"],
                    row["CODON POS"], row["EFFECT"])

        # reverse-complementing also reverses read direction, so the
        # flanking bases swap sides in addition to being complemented
        return (revcomp_class(cls), COMPLEMENT[row["THREE PRIME"]],
                COMPLEMENT[row["FIVE PRIME"]], row["CODON POS"], row["EFFECT"])

    folded = canon.apply(folded_key, axis=1, result_type="expand")
    folded.columns = KEYS
    canon = pd.concat([folded, canon[["N_MUTATIONS", "OPPORTUNITIES"]]], axis=1)
    canon = (canon.groupby(KEYS, as_index=False)
             .agg(N_MUTATIONS=("N_MUTATIONS", "sum"),
                  OPPORTUNITIES=("OPPORTUNITIES", "sum")))
    spectrum = "collapsed_6"
else:
    spectrum = "full_12"

# rate per opportunity, NaN (not divide-by-zero) where a cell has no opportunities
canon["RATE"] = canon["N_MUTATIONS"] / canon["OPPORTUNITIES"].where(canon["OPPORTUNITIES"] > 0)

canon["SPECTRUM"] = spectrum

canon = canon[KEYS + ["N_MUTATIONS", "OPPORTUNITIES", "RATE", "SPECTRUM"]]
canon = canon.sort_values(KEYS, key=lambda c: c.astype(str)).reset_index(drop=True)
canon.to_csv(out_path, sep="\t", index=False)