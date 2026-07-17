#!/usr/bin/env python3

import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

tables = snakemake.input.tables

sym_results = snakemake.output.sym_results
decision = snakemake.output.decision

KEYS = ["MUTATION_CLASS", "FIVE PRIME", "THREE PRIME", "CODON POS", "EFFECT"]

PAIRS = [("C>A", "G>T"), ("C>G", "G>C"), ("C>T", "G>A"),
         ("T>A", "A>T"), ("T>C", "A>G"), ("T>G", "A>C")]

MIN_TOTAL_N = 50

combined = pd.concat(
    [pd.read_csv(t, sep="\t", dtype={k: str for k in KEYS}) for t in tables],
    ignore_index=True,
)
canon = combined.groupby(KEYS, as_index=False).agg(
    N_MUTATIONS=("N_MUTATIONS", "sum"),
    OPPORTUNITIES=("OPPORTUNITIES", "first"),
)

by_class = (canon.groupby("MUTATION_CLASS").agg(N=("N_MUTATIONS", "sum"),
            OPP=("OPPORTUNITIES", "sum")))

rows = []
for fwd, rev in PAIRS:
    fwd_n = int(by_class.loc[fwd, "N"] if fwd in by_class.index else 0)
    rev_n = int(by_class.loc[rev, "N"] if rev in by_class.index else 0)
    fwd_opp = int(by_class.loc[fwd, "OPP"] if fwd in by_class.index else 0)
    rev_opp = int(by_class.loc[rev, "OPP"] if rev in by_class.index else 0)

    total_n = fwd_n + rev_n
    total_opp = fwd_opp + rev_opp

    if total_n == 0 or total_opp == 0:
        p_null = float("nan")
        p_val = float("nan")
    else:
        p_null = fwd_opp / total_opp
        p_val = binomtest(fwd_n, total_n, p=p_null, alternative="two-sided").pvalue

    fwd_rate = fwd_n / fwd_opp if fwd_opp else float("nan")
    rev_rate = rev_n / rev_opp if rev_opp else float("nan")

    rows.append({
        "pair": f"{fwd}/{rev}",
        "forward": fwd, "reverse": rev,
        "fwd_n": fwd_n, "rev_n": rev_n,
        "fwd_opp": fwd_opp, "rev_opp": rev_opp,
        "fwd_rate": fwd_rate, "rev_rate": rev_rate,
        "p_null": p_null,
        "p_value": p_val,
    })

sym = pd.DataFrame(rows)

testable = sym["p_value"].notna()
sym["p_bonferroni"] = float("nan")
if testable.any():
    sym.loc[testable, "p_bonferroni"] = multipletests(
        sym.loc[testable, "p_value"], method="bonferroni")[1]

sym["verdict"] = sym["p_bonferroni"].apply(
    lambda p: "untestable" if pd.isna(p) else
        ("symmetric" if p >= 0.05 else "asymmetric"))

sym.to_csv(sym_results, sep="\t", index=False)

total_n_pooled = int(sym["fwd_n"].sum() + sym["rev_n"].sum())
if total_n_pooled < MIN_TOTAL_N:
    decision_label = "insufficient_data"
else:
    any_asymmetric = (sym["verdict"] == "asymmetric").any()
    decision_label = "keep_12" if any_asymmetric else "collapse_to_6"

with open(decision, "w") as f:
    f.write("decision\t" + decision_label + "\n")
    f.write("total_n_pooled\t" + str(total_n_pooled) + "\n")
    f.write("n_pairs_tested\t" + str(int(testable.sum())) + "\n")
    f.write("n_asymmetric\t" + str(int((sym['verdict'] == 'asymmetric').sum())) + "\n")
    f.write("n_symmetric\t" + str(int((sym['verdict'] == 'symmetric').sum())) + "\n")
    f.write("n_untestable\t" + str(int((sym['verdict'] == 'untestable').sum())) + "\n")
