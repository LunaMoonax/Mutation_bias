#!/usr/bin/env python3

from collections import Counter
from Bio import SeqIO

ann_vcf_path = snakemake.input.ann_vcf
ref_path = snakemake.input.ref
out_path = snakemake.output.table
species = snakemake.wildcards.sp

BASES = ("A", "C", "G", "T")

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
    terms = eff.split("&")
    if any(t in ("synonymous_variant", "start_retained", "stop_retained_variant") for t in terms):
        return "syn"
    if any(t in ("missense_variant", "initiator_codon_variant", "stop_gained",
                 "stop_lost", "start_lost") for t in terms):
        return "nonsyn"
    return "intergenic"

def get_context(seqs, chrom, pos):
    seq = seqs[chrom]
    idx = pos - 1
    if idx < 1 or idx >= len(seq) - 1:
        return None
    five, base, three = seq[idx - 1], seq[idx], seq[idx + 1]
    if any(b not in BASES for b in (five, base, three)):
        return None
    return five, base, three

ref_seq = load_ref(ref_path)
seqs = {chrom: str(rec.seq).upper() for chrom, rec in ref_seq.items()}
counts  = Counter()
skipped = 0

with open(ann_vcf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        ctx = get_context(seqs, chrom, pos)
        if ctx is None:
            skipped += 1
            continue
        five, base, three = ctx

        ann_field = next((x for x in info.split(";") if x.startswith("ANN=")), None)
        if ann_field is None:
            continue
        first_ann = ann_field.split("=", 1)[1].split(",")[0]
        parts = first_ann.split("|")
        codon   = get_codon_pos(parts[12])
        effect  = get_effect_class(parts[1])

        counts[(f"{ref}>{alt}", five, three, codon, effect)] += 1

with open(out_path, "w") as out:
    out.write("MUTATION_CLASS\tFIVE PRIME\tTHREE PRIME\tCODON POS\tEFFECT\tOPPORTUNITIES\n")
    for key, n in sorted(counts.items(), key=lambda kv: tuple(str(x) for x in kv[0])):
        mut, five, three, codon, effect = key
        out.write(f"{mut}\t{five}\t{three}\t{codon}\t{effect}\t{n}\n")


import sys
sys.stderr.write(f"[{species}] cells: {len(counts)}, edge/N skipped: {skipped}\n")