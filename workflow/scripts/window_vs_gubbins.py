#!/usr/bin/env python3

import bisect
from collections import Counter

def read_bed(bed_file):
    pos = set()
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            pos.update(range(int(fields[1]) + 1, int(fields[2]) + 1))
    return pos

def read_gff(gff_file):
    intervals = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            intervals.append((int(fields[3]), int(fields[4])))
    intervals.sort()
    return intervals

def read_vcf(vcf_file):
    pos = set()
    minor_count = {}
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            p = int(fields[1])
            pos.add(p)

            genomes = fields[10:]
            n_alt = sum(1 for g in genomes if g == "1")
            n_ref = sum(1 for g in genomes if g == "0")
            minor_count[p] = min(n_alt, n_ref)
    return pos, minor_count

def interval_index(p, intervals, starts):
    i = bisect.bisect_right(starts, p) - 1
    if i >= 0 and intervals[i][0] <= p <= intervals[i][1]:
        return i
    return -1

def freq_spectrum(pos_set, minor_count):
    c = Counter()
    for p in pos_set:
        mc = minor_count.get(p, 0)
        if mc == 1:
            c["singleton"] += 1
        elif mc == 2:
            c["doubleton"] += 1
        elif mc == 3:
            c["tripleton"] += 1
        else:
            c["other"] += 1
    return c

window_all  = read_bed(snakemake.input.window_bed)
gubbins_all = read_gff(snakemake.input.gubbins_gff)
all_snps, minor_count = read_vcf(snakemake.input.vcf)
gub_start   = [s for s, _ in gubbins_all]

gubbins_snp_set  = {p for p in all_snps if interval_index(p, gubbins_all, gub_start) >= 0}
both_set         = window_all & gubbins_snp_set
window_only_set  = window_all - gubbins_snp_set
gubbins_only_set = gubbins_snp_set - window_all
neither_set      = all_snps - window_all - gubbins_snp_set

n_win             = len(window_all)
total_snps        = len(all_snps)
gubbins_snp_total = len(gubbins_snp_set)
both              = len(both_set)
window_only       = len(window_only_set)
gubbins_only      = len(gubbins_only_set)
neither           = len(neither_set)
snp_union         = both + window_only + gubbins_only
jaccard_snp       = both / snp_union if snp_union else None
recall_snp        = both / gubbins_snp_total if gubbins_snp_total else None
precision_snp     = both / n_win if n_win else None

def pct(x):
    return f"{x:.1%}" if x is not None else "n/a"

with open(snakemake.output.summary, "w") as f:
    def out(s):
        print(s)
        f.write(s + "\n")
    out(f"species                       : {snakemake.wildcards.sp}")
    out(f"total SNP sites (VCF)         : {total_snps}")
    out(f"window-flagged SNPs           : {n_win}  ({pct(n_win/total_snps if total_snps else None)})")
    out(f"gubbins-region SNPs           : {gubbins_snp_total}  ({pct(gubbins_snp_total/total_snps if total_snps else None)})")
    out(f"  flagged by both methods     : {both}  ({pct(both/total_snps if total_snps else None)})")
    out(f"  flagged by window only      : {window_only}  ({pct(window_only/total_snps if total_snps else None)})")
    out(f"  flagged by gubbins only     : {gubbins_only}  ({pct(gubbins_only/total_snps if total_snps else None)})")
    out(f"  flagged by neither          : {neither}  ({pct(neither/total_snps if total_snps else None)})")
    out(f"precision (window -> gubbins) : {pct(precision_snp)}")
    out(f"recall (gubbins found by win) : {pct(recall_snp)}")
    out(f"Jaccard agreement             : {pct(jaccard_snp)}")

    out("")
    out("minor-allele frequency spectrum (singleton/doubleton/tripleton/other)")
    for label, pos_set in [
        ("all SNPs",       all_snps),
        ("window-flagged", window_all),
        ("gubbins-region", gubbins_snp_set),
        ("both methods",   both_set),
        ("window only",    window_only_set),
        ("gubbins only",   gubbins_only_set),
    ]:
        c = freq_spectrum(pos_set, minor_count)
        n = len(pos_set)
        out(f"  {label:16s}: n={n}"
            f"  singleton={c['singleton']} ({pct(c['singleton']/n if n else None)})"
            f"  doubleton={c['doubleton']} ({pct(c['doubleton']/n if n else None)})"
            f"  tripleton={c['tripleton']} ({pct(c['tripleton']/n if n else None)})"
            f"  other={c['other']} ({pct(c['other']/n if n else None)})")
