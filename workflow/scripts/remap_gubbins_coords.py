#!/usr/bin/env python3

# gubbins output (recombination GFF, ancestral-reconstruction EMBL) is in xmfa
# alignment coordinates; remap both back to the reference genome's own coordinates

# find the reference's index/id/length in the xmfa header (file ending ".ref")
def find_reference(xmfa_file):
    ref_index = None
    ref_seqid = None
    ref_length = None
    cur_index = None
    with open(xmfa_file) as f:
        for line in f:
            if line.startswith("##SequenceIndex"):
                cur_index = int(line.split()[1])
            elif line.startswith("##SequenceFile"):
                fname = line.split(None, 1)[1].strip()
                if fname.endswith(".ref"):
                    ref_index = cur_index
            elif line.startswith("##SequenceHeader") and cur_index == ref_index and ref_seqid is None:
                ref_seqid = line.split(None, 1)[1].strip().lstrip(">")
            elif line.startswith("##SequenceLength") and cur_index == ref_index and ref_length is None:
                ref_length = int(line.split()[1].rstrip("bp\n"))
            elif not line.startswith("#"):
                break
    return ref_index, ref_seqid, ref_length

# build alignment-column -> reference-position array from the reference's own
# blocks in the xmfa (None at gap columns)
def build_alignment_to_ref_map(xmfa_file, ref_index):
    prefix = f"> {ref_index}:"
    mapping = []
    in_ref_block = False
    ref_pos = None
    with open(xmfa_file) as f:
        for line in f:
            if line.startswith(">"):
                if line.startswith(prefix):
                    in_ref_block = True
                    coords = line[len(prefix):].split()[0]
                    ref_pos = int(coords.split("-")[0])
                else:
                    in_ref_block = False
                continue
            if line.startswith("="):
                in_ref_block = False
                continue
            if in_ref_block:
                for ch in line.strip():
                    if ch == "-":
                        mapping.append(None)
                    else:
                        mapping.append(ref_pos)
                        ref_pos += 1
    return mapping

# nearest non-gap ref position at/after i, for snapping interval starts
def nearest_forward(mapping, i):
    while i < len(mapping) and mapping[i] is None:
        i += 1
    return mapping[i] if i < len(mapping) else None

# nearest non-gap ref position at/before i, for snapping interval ends
def nearest_backward(mapping, i):
    while i >= 0 and mapping[i] is None:
        i -= 1
    return mapping[i] if i >= 0 else None

# exact single-position lookup (no snapping) for point positions like SNPs
def remap_pos(mapping, aln_pos):
    i = aln_pos - 1
    if 0 <= i < len(mapping):
        return mapping[i]
    return None

ref_index, ref_seqid, ref_length = find_reference(snakemake.input["xmfa"])
mapping = build_alignment_to_ref_map(snakemake.input["xmfa"], ref_index)

# remap the recombination-prediction GFF intervals, snapping to nearest ref base at gaps
with open(snakemake.input["gff"], "r") as fin, open(snakemake.output["gff"], "w") as fout:
    fout.write("##gff-version 3\n")
    fout.write(f"##sequence-region {ref_seqid} 1 {ref_length}\n")
    for line in fin:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        aln_start, aln_end = int(fields[3]), int(fields[4])
        ref_start = nearest_forward(mapping, aln_start - 1)
        ref_end   = nearest_backward(mapping, aln_end - 1)
        if ref_start is None or ref_end is None or ref_start > ref_end:
            continue
        fields[0] = ref_seqid
        fields[3] = str(ref_start)
        fields[4] = str(ref_end)
        fout.write("\t".join(fields) + "\n")

# remap the ancestral-reconstruction EMBL's variation positions; drop (and skip the
# rest of) any variation record that lands on a reference gap
n_total = n_kept = n_gap = 0
with open(snakemake.input["embl"], "r") as fin, open(snakemake.output["embl"], "w") as fout:
    skip_record = False
    for line in fin:
        if line.startswith("FT") and "variation" in line.split():
            n_total += 1
            aln_pos = int(line.split()[-1])
            ref_pos = remap_pos(mapping, aln_pos)
            if ref_pos is None:
                n_gap += 1
                skip_record = True
                continue
            skip_record = False
            n_kept += 1
            fout.write(line.replace(f" {aln_pos}", f" {ref_pos}"))
        elif skip_record and line.startswith("FT") and "variation" not in line.split():
            continue
        else:
            skip_record = False
            fout.write(line)

print(f"branch SNPs: {n_kept}/{n_total} remapped ({n_gap} on reference gaps, dropped)")