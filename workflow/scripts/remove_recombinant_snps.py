#!/usr/bin/env python3

# drop any SNP that falls inside a gubbins-predicted recombination region

from collections import defaultdict

def read_gff(gff_file):
    intervals = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            intervals[fields[0]].append((int(fields[3]), int(fields[4])))
    return intervals

gubbins_regions = read_gff(snakemake.input.gff)

vcf_file = snakemake.input.vcf
recombination_free_vcf = snakemake.output.rec_free_vcf

n_total = n_kept = 0

with open(vcf_file) as fin, open(recombination_free_vcf, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue

        n_total += 1
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        pos = int(fields[1])

        in_recombination = any(start <= pos <= end for start, end in gubbins_regions.get(chrom, ()))

        if not in_recombination:
            fout.write(line)
            n_kept += 1

print(f"{n_kept} / {n_total} {n_kept/n_total:.1%} variants kept (outside Gubbins recombination regions)")