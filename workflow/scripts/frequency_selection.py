#!/usr/bin/env python3

vcf_file = snakemake.input.vcf

singleton_vcf = open(snakemake.output.singleton_vcf, "w")
doubleton_vcf = open(snakemake.output.doubleton_vcf, "w")
tripleton_vcf = open(snakemake.output.tripleton_vcf, "w")
vcf_by_count = {1: singleton_vcf, 2: doubleton_vcf, 3: tripleton_vcf}

header = ["CHROM", "POS", "SAMPLES", "REF", "ALT", "MUTATION_CLASS"]
names = {1: "singleton", 2: "doubleton", 3: "tripleton"}
tables = {1: snakemake.output.singleton_table, 2: snakemake.output.doubleton_table,
          3: snakemake.output.tripleton_table}

singleton_rows, doubleton_rows, tripleton_rows = [], [], []
rows_by_count = {1: singleton_rows, 2: doubleton_rows, 3: tripleton_rows}

n_total = 0
counts = {1: 0, 2: 0, 3: 0}

with open(vcf_file) as fin:
    for line in fin:
        if line.startswith("#CHROM"):
            header_fields = line.rstrip("\n").split("\t")
            sample_names = header_fields[10:]
            for f in vcf_by_count.values():
                f.write(line)
            continue
        if line.startswith("#"):
            for f in vcf_by_count.values():
                f.write(line)
            continue

        n_total += 1
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        samples = fields[10:]

        n_carriers = samples.count("1")
        if n_carriers not in vcf_by_count:
            continue

        counts[n_carriers] += 1
        vcf_by_count[n_carriers].write(line)
        carriers = ",".join(n for n, s in zip(sample_names, samples) if s == "1")
        rows_by_count[n_carriers].append([chrom, pos, carriers, ref, alt, f"{ref}>{alt}"])

for f in vcf_by_count.values():
    f.close()

for n_carriers, rows in rows_by_count.items():
    pct = counts[n_carriers] / n_total if n_total else 0
    print(f"{counts[n_carriers]} / {n_total} ({pct:.1%}) are {names[n_carriers]}s")
    with open(tables[n_carriers], "w") as ftab:
        ftab.write(f"# {counts[n_carriers]}/{n_total} ({pct:.1%}) are {names[n_carriers]}s\n")
        ftab.write("\t".join(header) + "\n")
        for row in rows:
            ftab.write("\t".join(row) + "\n")
