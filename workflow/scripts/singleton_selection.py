#!/usr/bin/env python3

vcf_file = snakemake.input.vcf
singleton_vcf = snakemake.output.singleton_vcf
singleton_table_file = snakemake.output.singleton_table

def is_singleton(samples):
    alt_count = samples.count("1")
    if alt_count == 1:
        return True
    return False

singleton_table = []
header = ["POS", "SAMPLE", "REF", "ALT", "MUTATION_CLASS"]

n_total = n_singletons = 0
with open(vcf_file) as fin, open(singleton_vcf, "w") as fout:
    for line in fin:
        if line.startswith("#CHROM"):
            header_fields = line.rstrip("\n").split("\t")
            sample_names = header_fields[10:]
            fout.write(line)
            continue
        if line.startswith("#"):
            fout.write(line)
            continue
        
        n_total += 1
        fields = line.rstrip("\n").split("\t")
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        samples = fields[10:]

        if is_singleton(samples):
            fout.write(line)
            n_singletons += 1
            singleton_table.append([pos, sample_names[samples.index("1")], ref, alt, f"{ref}>{alt}"])

print(f"{n_singletons} / {n_total} {n_singletons/n_total:.1%} variants are singletons from recent SNPs")

with open(singleton_table_file, "w") as ftab:
    ftab.write("# " + str(n_singletons) + "/" + str(n_total) + " " + f"{n_singletons/n_total:.1%}" + " variants are singletons\n")
    ftab.write("\t".join(header) + "\n")
    for row in singleton_table:
        ftab.write("\t".join(row) + "\n")