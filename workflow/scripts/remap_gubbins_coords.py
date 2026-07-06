#!/usr/bin/env python3

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

def nearest_forward(mapping, i):
    while i < len(mapping) and mapping[i] is None:
        i += 1
    return mapping[i] if i < len(mapping) else None

def nearest_backward(mapping, i):
    while i >= 0 and mapping[i] is None:
        i -= 1
    return mapping[i] if i >= 0 else None

ref_index, ref_seqid, ref_length = find_reference(snakemake.input.xmfa)
mapping = build_alignment_to_ref_map(snakemake.input.xmfa, ref_index)

with open(snakemake.input.gff) as fin, open(snakemake.output.gff, "w") as fout:
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
