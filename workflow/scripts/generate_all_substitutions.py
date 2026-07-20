#!/usr/bin/env python3

# build an exhaustive VCF of every possible substitution at every ref position (opportunity space)

from Bio import SeqIO

ref_path = snakemake.input.ref
all_sub_vcf_path = snakemake.output.vcf

BASES = ("A", "C", "G", "T")

def load_ref(ref_path):
    return SeqIO.to_dict(SeqIO.parse(ref_path, "fasta"))

ref_seq = load_ref(ref_path)

with open(all_sub_vcf_path, "w") as fout:
    fout.write("##fileformat=VCFv4.2\n")
    fout.write("##Generated all-substitutions VCF from the reference FASTA\n")
    for chrom, rec in ref_seq.items():
        fout.write(f"##contig=<ID={chrom},length={len(rec.seq)}>\n")
    fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # for every position, write the 3 non-ref bases as separate substitution records
    for chrom, rec in ref_seq.items():
        seq = str(rec.seq).upper()
        for i, ref in enumerate(seq):
            if ref not in BASES:
                continue
            pos = i + 1
            for alt in BASES:
                if alt != ref:
                    fout.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t\n")
