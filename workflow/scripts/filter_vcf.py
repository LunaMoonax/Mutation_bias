#!/usr/bin/env python3

BASES = set("ACGT")

vcf_file = snakemake.input.vcf
filtered_vcf = snakemake.output.vcf

n_total = n_kept = 0

with open(vcf_file) as fin, open(filtered_vcf, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue

        n_total += 1
        cols = line.rstrip("\n").split("\t")
        filt = cols[6]
        ref  = cols[3]
        alt  = cols[4]

        if filt != "PASS":
            continue
        if "," in alt:
            continue
        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref not in BASES or alt not in BASES:
            continue

        fout.write(line)
        n_kept += 1

print(f"{n_kept} / {n_total} variants kept (PASS + biallelic SNP)")