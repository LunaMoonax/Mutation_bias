#!/usr/bin/env python3

from Bio import SeqIO

singleton_table_path = snakemake.input.singleton_table
singleton_ann_vcf_path = snakemake.input.singleton_ann_vcf
doubleton_table_path  = snakemake.input.doubleton_table
doubleton_ann_vcf_path  = snakemake.input.doubleton_ann_vcf
tripleton_table_path  = snakemake.input.tripleton_table
tripleton_ann_vcf_path  = snakemake.input.tripleton_ann_vcf
ref_path = snakemake.input.ref
singleton_ann_table_path  = snakemake.output.singleton_ann_table
doubleton_ann_table_path  = snakemake.output.doubleton_ann_table
tripleton_ann_table_path  = snakemake.output.tripleton_ann_table

BASES = set("ACGT")

def load_ref(ref_path):
    return SeqIO.to_dict(SeqIO.parse(ref_path, "fasta"))

def get_codon_pos(cds_pos):
    if not cds_pos or cds_pos == "":
        return "intergenic"
    cds = cds_pos.split("/")[0]
    if not cds.isdigit():
        return "intergenic"
    return ((int(cds) - 1) % 3) + 1

def get_effect_class(eff):
    if eff in ("synonymous_variant", "start_retained", "stop_retained_variant"):
        return "syn"
    if eff in ("missense_variant", "initiator_codon_variant", "stop_gained",
               "stop_lost", "start_lost"):
        return "nonsyn"
    return "intergenic"

def parse_ann_vcf(ann_vcf_path):
    out = {}
    with open(ann_vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]

            ann_field = next((x for x in info.split(";") if x.startswith("ANN=")), None)
            if ann_field is None:
                continue
            first_ann = ann_field.split("=", 1)[1].split(",")[0]
            parts = first_ann.split("|")

            eff = parts[1]
            cds_pos = parts[12]
            out[(chrom, pos)] = (get_codon_pos(cds_pos), get_effect_class(eff))
    return out

def read_table(table_path):
    with open(table_path) as f:
        lines = f.readlines()
    comment = lines[0] if lines[0].startswith("#") else None
    start = 1 if comment else 0
    header = lines[start].rstrip("\n").split("\t")
    rows = [l.rstrip("\n").split("\t") for l in lines[start + 1:]]
    return comment, header, rows


def get_context(ref_seq, chrom, pos):
    seq = str(ref_seq[chrom].seq).upper()
    idx = pos - 1
    if idx < 1 or idx >= len(seq) - 1:
        return None
    five, base, three = seq[idx - 1], seq[idx], seq[idx + 1]
    if any(b not in BASES for b in (five, base, three)):
        return None
    return five, base, three

def annotate(table_path, ann_vcf_path, out_path, ref_seq):
    comment, header, rows = read_table(table_path)
    ann = parse_ann_vcf(ann_vcf_path)
    pos_idx = header.index("POS")
    chrom_idx = header.index("CHROM")

    new_header = header + ["FIVE PRIME", "THREE PRIME", "CODON POS", "EFFECT"]
    n_written = n_skipped = 0

    with open(out_path, "w") as fout:
        if comment:
            fout.write(comment)
        
        fout.write("\t".join(new_header) + "\n")

        for row in rows:
            pos = int(row[pos_idx])
            chrom = row[chrom_idx]

            context = get_context(ref_seq, chrom, pos)
            if context is None:
                n_skipped += 1
                continue
            five, _, three = context

            codon_pos, effect = ann.get((chrom, pos), ("NA", "NA"))

            fout.write("\t".join(row + [five, three, str(codon_pos), effect]) + "\n")
            n_written += 1
    
    print(f"{out_path}: {n_written} annotated, {n_skipped} skipped (contig edge)")

ref_seq = load_ref(ref_path)

for table_in, ann_vcf_in, table_out in [
    (singleton_table_path, singleton_ann_vcf_path, singleton_ann_table_path),
    (doubleton_table_path, doubleton_ann_vcf_path, doubleton_ann_table_path),
    (tripleton_table_path, tripleton_ann_vcf_path, tripleton_ann_table_path),
]:
    annotate(table_in, ann_vcf_in, table_out, ref_seq)