#!/usr/bin/env python3

import subprocess
from collections import defaultdict
import sys
import os
import shutil


# ---------------- FASTA ----------------
def read_fasta(path):
    seqs = {}
    current_id = None
    buf = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(buf)
                current_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)

        if current_id:
            seqs[current_id] = "".join(buf)

    return seqs


def write_fasta(records, path):
    with open(path, "w") as f:
        for rid, seq in records.items():
            f.write(f">{rid}\n{seq}\n")


# ---------------- REVERSE COMPLEMENT ----------------
def reverse_complement(seq):
    comp = {
        "A": "T", "T": "A",
        "C": "G", "G": "C",
        "a": "t", "t": "a",
        "c": "g", "g": "c"
    }
    return "".join(comp.get(b, b) for b in reversed(seq))


# ---------------- MINIMAP2 ----------------
def run_minimap2(edges_fasta, reference, paf_file, minimap2="minimap2"):
    with open(paf_file, "w") as out:
        subprocess.run([minimap2, reference, edges_fasta],
                       stdout=out, check=True)


# ---------------- PAF PARSER ----------------
def parse_paf(paf_file, min_len=1000):
    alignments = defaultdict(list)

    with open(paf_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            edge = cols[0]  
            strand = cols[4]
            ref_start = int(cols[7])  
            ref_end = int(cols[8])    
            aln_len = int(cols[10])

            if aln_len < min_len:
                continue

            alignments[edge].append({
                "start": ref_start,
                "end": ref_end,
                "strand": strand,
                "len": aln_len
            })

    return alignments


# ---------------- EDGE ORIENTATION ----------------
def orient_edges(edges_fasta, paf_file, output_fasta):
    sequences = read_fasta(edges_fasta)
    orientations = defaultdict(set)

    with open(paf_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue

            edge = cols[0]  
            strand = cols[4]
            orientations[edge].add(strand)

    multi_oriented = [e for e, d in orientations.items() if len(d) > 1]
    single_oriented = [e for e, d in orientations.items() if len(d) == 1]

    if len(multi_oriented) != 1 or len(single_oriented) != 2:
        raise RuntimeError(
            f"Unexpected orientation pattern: multi={multi_oriented}, single={single_oriented}"
        )

    target_dir = "+"

    for record_id, seq in sequences.items():
        if record_id in single_oriented:
            current_dir = list(orientations[record_id])[0]
            if current_dir != target_dir:
                sequences[record_id] = reverse_complement(seq)

    write_fasta(sequences, output_fasta)
    return sequences


# ---------------- CLASSIFY ----------------
def classify_edges(sequences, alignments):
    lengths = sorted(sequences.items(), key=lambda x: len(x[1]))
    LSC_edge = lengths[-1][0]

    IR_candidates = [e for e in alignments if len(alignments[e]) >= 2]
    if len(IR_candidates) != 1:
        raise RuntimeError(f"Expected 1 IR edge, got {IR_candidates}")

    IR_edge = IR_candidates[0]

    SSC_candidates = [e for e in sequences if e not in (LSC_edge, IR_edge)]
    if len(SSC_candidates) != 1:
        raise RuntimeError(f"Expected 1 SSC edge, got {SSC_candidates}")

    SSC_edge = SSC_candidates[0]

    return LSC_edge, IR_edge, SSC_edge


# ---------------- IR ORIENTATION ----------------
def orient_ir(IR_edge, alignments, sequences, debug=False):
    if IR_edge not in alignments:
        raise KeyError(f"{IR_edge} not found in alignments")

    ir_aligns = list(alignments[IR_edge])

    if len(ir_aligns) < 2:
        raise RuntimeError(
            f"Expected at least 2 alignments for IR edge '{IR_edge}', got {len(ir_aligns)}"
        )

    ir_aligns.sort(key=lambda x: (x["start"], x["end"]))
    IRa = ir_aligns[0]

    if debug:
        print(f"\n[DEBUG] IR edge: {IR_edge}")
        print("[DEBUG] IR alignments sorted by reference position:")
        for i, aln in enumerate(ir_aligns, 1):
            print(
                f"  {i}: start={aln['start']}, end={aln['end']}, "
                f"strand={aln['strand']}, len={aln['len']}"
            )
        print(f"[DEBUG] Selected IRa: {IRa}")

    IR_seq = sequences[IR_edge]

    if IRa["strand"] == "-":
        if debug:
            print("[DEBUG] IRa on '-' strand -> reverse complementing")
        IR_seq = reverse_complement(IR_seq)
    else:
        if debug:
            print("[DEBUG] IRa on '+' strand -> keeping original orientation")

    IR_seq_RC = reverse_complement(IR_seq)

    return IR_seq, IR_seq_RC


# ---------------- MERGE ----------------
def merge_edges(sequences, alignments, outputdir, sample_name):
    mapped_edges = set(alignments.keys())
    sequences = {e: s for e, s in sequences.items() if e in mapped_edges}

    if len(sequences) < 3:
        raise RuntimeError(f"Too few mapped edges: {len(sequences)}")

    LSC_edge, IR_edge, SSC_edge = classify_edges(sequences, alignments)

    LSC_seq = sequences[LSC_edge]
    SSC_seq = sequences[SSC_edge]
    IR_seq, IR_seq_RC = orient_ir(IR_edge, alignments, sequences)

    chloroplast = LSC_seq + IR_seq + SSC_seq + IR_seq_RC

    os.makedirs(outputdir, exist_ok=True)
    outpath = os.path.join(outputdir, f"{sample_name}_chloroplast.fasta")

    with open(outpath, "w") as f:
        f.write(f">{sample_name}_chloroplast\n{chloroplast}\n")


# ---------------- MAIN ----------------
def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} edges.fa reference.fa sample_name")
        sys.exit(1)

    edges_fasta = sys.argv[1]
    reference = sys.argv[2]
    sample_name = sys.argv[3]

    workdir = os.path.dirname(os.path.abspath(edges_fasta))

    paf_file = os.path.join(workdir, "orientation.paf")
    output_fasta = os.path.join(workdir, "edges_oriented.fa")
    outputdir = os.path.join(workdir, "chloroplast_final_assembly")

    if not shutil.which("minimap2"):
        raise RuntimeError("minimap2 not found in PATH")

    run_minimap2(edges_fasta, reference, paf_file)

    alignments = parse_paf(paf_file)
    if not alignments:
        raise RuntimeError("No valid alignments found")

    sequences = orient_edges(edges_fasta, paf_file, output_fasta)
    merge_edges(sequences, alignments, outputdir, sample_name)


if __name__ == "__main__":
    main()
