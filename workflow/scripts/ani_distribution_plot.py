#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

FLOOR = snakemake.params.floor
MAD_K = snakemake.params.mad_k
MIN_DROP = snakemake.params.min_drop

def genome_id (p):
    n = os.path.basename(str(p).strip())
    for ext in (".gz", ".fasta", ".fna", ".fa"):
        if n.endswith(ext):
            n = n[: -len(ext)]
    return n

def specie_ani (edge_list, genome_list):
    df = pd.read_csv(edge_list, sep = "\t")
    genome_ref = df.columns[0]
    genome_qry = df.columns[1]
    ani = df.columns[2]

    with open(genome_list) as f:
        all_genomes = {genome_id(i) for i in f if i.strip()}
    long = pd.concat([
        df[[genome_ref, 
            ani]].rename(columns={genome_ref: "genome", ani: "ani"}),
        df[[genome_qry,
            ani]].rename(columns={genome_qry: "genome", ani: "ani"}),
    ])
    
    long["genome"] = long["genome"].map(genome_id)
    med = (long.groupby("genome", as_index=False)[["ani"]]
               .median().rename(columns={"ani": "median_ani"})
    )
    
    missing = sorted(all_genomes - set(med["genome"]))
    if missing:
        med = pd.concat([med, pd.DataFrame({"genome": missing, 
                                             "median_ani" : np.nan})])

    present_med = med["median_ani"].dropna()
    species_median = present_med.median()
    mad = (present_med - species_median).abs().median()
    fence = species_median - MAD_K * 1.4826 * mad

    cond = [
        med["median_ani"].isna(),
        med["median_ani"] < FLOOR,
        (med["median_ani"] < fence) &
        (species_median - med["median_ani"] >= MIN_DROP),
    ]
    choice = ["no_comparisons", "below_floor_95", "distant_outlier"]
    med["reason"] = np.select(cond, choice, default="")
    med["exclude"] = med["reason"] != ""
    med["mad_fence"] = fence
    med["species_median"] = species_median

    return (med.sort_values("median_ani").reset_index(drop=True), 
           df[ani].to_numpy())

def plot_ani_distribution (all_med, output_dir):
    d = all_med.dropna(subset=["median_ani"]).copy()
    d["rank"] = d.groupby("species")["median_ani"].rank(method="first")
    g = sns.relplot(
        data=d, x="rank", y="median_ani", hue="exclude",
        col="species", col_wrap=5, height=2.4, s=18, legend=False,
        palette={False: "black", True: "red"},
        edgecolor="none", alpha=0.6,
        facet_kws={"sharex": False, "sharey": False},
    )
    for sp, ax in g.axes_dict.items():
        ax.axhline(FLOOR, color="blue", ls="--", lw=1)
        fence = all_med.loc[all_med["species"] == sp, "mad_fence"].iloc[0]
        ax.axhline(fence, color="green", ls=":", lw=1)
        drop = all_med.loc[all_med["species"] == sp, "species_median"].iloc[0] - MIN_DROP
        ax.axhline(drop, color="orange", ls=":", lw=1)
    g.set_titles("{col_name}")
    g.figure.suptitle("Per-genome median ANI\n"
        "Blue dashed line: 95% ANI floor. Green dotted line: MAD fence. Orange dotted line: drop threshold.",
        y=1.02)
    g.set_axis_labels("genome (ordered by median ANI)", "median ANI")
    g.savefig(f"{output_dir}/ani_distribution.png",
              dpi=300, bbox_inches="tight")
    plt.close(g.figure)

def plot_ani_violin (all_ani, output_dir, cols_per_row=5):
    n = all_ani["species"].nunique()
    n_rows = (n + cols_per_row - 1) // cols_per_row
    g = sns.catplot(
        data=all_ani, x="species", y="ani", kind="violin",
        col="species", col_wrap=cols_per_row,
        cut=0, height=3.5, aspect=0.9,
        sharey=True, sharex=False,
    )
    g.set_titles("{col_name}")
    g.set_axis_labels("", "pairwise ANI")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelbottom=False)
    g.figure.suptitle("Distribution of pairwise ANI values per species", y=1.01)
    g.savefig(f"{output_dir}/ani_violin.png", dpi=300, bbox_inches="tight")
    plt.close(g.figure)
    
edge_files = list(snakemake.input.edge_list)
input_list = list(snakemake.input.genome_list)
output_dir = str(snakemake.params.output_dir)

lists = {os.path.basename(i).replace("_list.txt", ""): i for i in input_list}

all_med, all_ani = [], []

for edge_file in sorted(edge_files):
    species = os.path.basename(edge_file).replace("_ani_edge_list.txt", "")
    med, pw = specie_ani(edge_file, lists[species])
    med["species"] = species
    all_med.append(med)

    all_ani.append(pd.DataFrame({"species": species, "ani": pw}))
    print(f"{species:30s} kept {(~med['exclude']).sum():4d}  excluded {med['exclude'].sum():3d}")

plot_ani_distribution(pd.concat(all_med, ignore_index=True), output_dir)
plot_ani_violin(pd.concat(all_ani, ignore_index=True), output_dir)

all_med_df = pd.concat(all_med, ignore_index=True)
all_med_df.to_csv(f"{output_dir}/ani_qc.csv", index=False)