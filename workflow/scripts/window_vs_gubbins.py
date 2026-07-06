#!/usr/bin/env python3

import bisect

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

def read_vcf_positions(vcf_file):
    pos = set()
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            pos.add(int(fields[1]))
    return pos

def interval_index(p, intervals, starts):
    i = bisect.bisect_right(starts, p) - 1
    if i >= 0 and intervals[i][0] <= p <= intervals[i][1]:
        return i
    return -1

window_all  = read_bed(snakemake.input.window_bed)
gubbins_all = read_gff(snakemake.input.gubbins_gff)
all_snps    = read_vcf_positions(snakemake.input.vcf)
gub_start   = [s for s, _ in gubbins_all]

in_gubbins = sum(1 for p in window_all if interval_index(p, gubbins_all, gub_start) >= 0)
gubbins_snp_total = sum(1 for p in all_snps if interval_index(p, gubbins_all, gub_start) >= 0)

n_win         = len(window_all)
total_snps    = len(all_snps)
both          = in_gubbins
window_only   = n_win - both
gubbins_only  = gubbins_snp_total - both
neither       = total_snps - both - window_only - gubbins_only
snp_union     = both + window_only + gubbins_only
jaccard_snp   = both / snp_union if snp_union else None
recall_snp    = both / gubbins_snp_total if gubbins_snp_total else None
precision_snp = both / n_win if n_win else None

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
