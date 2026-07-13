#!/usr/bin/env python3

vcf_file = snakemake.input.vcf

CLASS_BY_ALT_COUNT = {1: "singleton", 2: "doubleton", 3: "tripleton"}

outputs = {
    "singleton": (snakemake.output.singleton_vcf, snakemake.output.singleton_table),
    "doubleton": (snakemake.output.doubleton_vcf, snakemake.output.doubleton_table),
    "tripleton": (snakemake.output.tripleton_vcf, snakemake.output.tripleton_table),
}

header = ["POS", "SAMPLES", "REF", "ALT", "MUTATION_CLASS"]

vcf_out_files = {cls: open(vcf_path, "w") for cls, (vcf_path, _) in outputs.items()}
tables = {cls: [] for cls in outputs}
counts = {cls: 0 for cls in outputs}
n_total = 0

with open(vcf_file) as fin:
    for line in fin:
        if line.startswith("#CHROM"):
            header_fields = line.rstrip("\n").split("\t")
            sample_names = header_fields[10:]
            for f in vcf_out_files.values():
                f.write(line)
            continue
        if line.startswith("#"):
            for f in vcf_out_files.values():
                f.write(line)
            continue

        n_total += 1
        fields = line.rstrip("\n").split("\t")
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        samples = fields[10:]

        cls = CLASS_BY_ALT_COUNT.get(samples.count("1"))
        if cls is None:
            continue

        vcf_out_files[cls].write(line)
        counts[cls] += 1
        carriers = ",".join(name for name, s in zip(sample_names, samples) if s == "1")
        tables[cls].append([pos, carriers, ref, alt, f"{ref}>{alt}"])

for f in vcf_out_files.values():
    f.close()

for cls, (_, table_path) in outputs.items():
    n = counts[cls]
    pct = n / n_total if n_total else 0
    print(f"{n} / {n_total} {pct:.1%} variants are {cls}s from recent SNPs")
    with open(table_path, "w") as ftab:
        ftab.write(f"# {n}/{n_total} {pct:.1%} variants are {cls}s\n")
        ftab.write("\t".join(header) + "\n")
        for row in tables[cls]:
            ftab.write("\t".join(row) + "\n")
