#!/usr/bin/env python3

from collections import defaultdict

def read_snps_positions(vcf_file):
    per_genome = defaultdict(lambda: defaultdict(list))

    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            genomes = fields[10:]

            n_alt = sum(1 for g in genomes if g == "1")
            n_ref = sum(1 for g in genomes if g == "0")

            if n_alt == 0 and n_ref == 0:
                continue

            minor_allele = "1" if n_alt < n_ref else "0"

            for i, g in enumerate(genomes):
                if g == minor_allele:
                    per_genome[i][chrom].append(pos)

    return per_genome

def flagged_snps(per_genome, window, threshold):
    flagged = defaultdict(list)

    for chrom_dict in per_genome.values():
        for chrom, pos in chrom_dict.items():
            pos.sort()
            left = 0
            for right in range(len(pos)):
                while pos[right] - pos[left] > window:
                    left += 1
                if right - left + 1 >= threshold:
                    for i in range(left, right + 1):
                        flagged[chrom].append(pos[i])

    return flagged

vcf_file = snakemake.input.vcf
window = snakemake.params.window
threshold = snakemake.params.threshold
output_file = snakemake.output.flagged

per_genome = read_snps_positions(vcf_file)
flagged = flagged_snps(per_genome, window, threshold)

with open(output_file, "w") as f:
    f.write("#chrom\tstart\tend\n")
    for chrom in sorted(flagged):
        start = prev = None
        for p in sorted(flagged[chrom]):
            if start is None:
                start = prev = p
            elif p <= prev + 1:
                prev = p
            else:
                f.write(f"{chrom}\t{start - 1}\t{prev}\n")
                start = prev = p
        if start is not None:
            f.write(f"{chrom}\t{start - 1}\t{prev}\n")


