#!/usr/bin/env python3

import os
import shutil
from datetime import datetime
from collections import Counter

from Bio import SeqIO
from io_utils import load_genetic_code
from Bio.Seq import Seq   

# ============================================================
# CODON COUNTS
# ============================================================
def calculate_codon_counts(seqs):

    counts = Counter()

    for seq in seqs:
        seq = str(seq).upper()

        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]

            if len(codon) == 3:
                counts[codon] += 1

    return counts


# ============================================================
# BEST CODON PER AMINOZUUR
# ============================================================
def get_best_codons(reference_fasta, codon2aa):

    sequences = [str(r.seq) for r in SeqIO.parse(reference_fasta, "fasta")]

    codon_counts = calculate_codon_counts(sequences)

    aa_groups = {}

    for codon, count in codon_counts.items():

        aa = codon2aa.get(codon)

        if aa is None:
            continue

        if aa not in aa_groups:
            aa_groups[aa] = {}

        aa_groups[aa][codon] = count

    best_codons = {}

    for aa, codons in aa_groups.items():

        # choose codons with highest percentages
        best = max(codons, key=codons.get)

        best_codons[aa] = best

    return best_codons


# ============================================================
# SEQUENCE OPTIMIZATION
# ============================================================
def optimize_sequence(seq, codon2aa, best_codons):

    seq = str(seq).upper()
    new_seq = []

    for i in range(0, len(seq), 3):

        codon = seq[i:i+3]

        if len(codon) != 3:
            new_seq.append(codon)
            continue

        aa = codon2aa.get(codon)

        if aa in best_codons:
            new_seq.append(best_codons[aa])
        else:
            new_seq.append(codon)

    return "".join(new_seq)


# ============================================================
# MAIN PIPELINE
# ============================================================
def optimize_fasta_pipeline(
    input_fasta,
    reference_choice="open",  # open, ncbi, target
    results_dir=None
):

    # ========================================================
    # SETUP
    # ========================================================
    if results_dir is None:
        timestamp = datetime.now().strftime("%y%m%d_%H%M")
        results_dir = f"results_{timestamp}"

    os.makedirs(results_dir, exist_ok=True)

    input_dir = os.path.join(results_dir, "input")
    os.makedirs(input_dir, exist_ok=True)

    shutil.copy(input_fasta, input_dir)

    # ========================================================
    # CHOOSE REFERENCE
    # ========================================================
    if reference_choice == "open":
        reference_fasta = "op_cds_chr1-4.fasta"

    elif reference_choice == "ncbi":
        reference_fasta = "NCBI_cds.fasta"

    elif reference_choice == "target":
        reference_fasta = input_fasta

    else:
        raise ValueError("reference_choice has to be 'open', 'ncbi' or 'target'. ")

    print(f"\n🔬 Running optimization ({reference_choice})")
    print(f"📁 Reference: {reference_fasta}")

    codon2aa = load_genetic_code()

    # ========================================================
    # BEST CODONS
    # ========================================================
    best_codons = get_best_codons(reference_fasta, codon2aa)

    print("\n✅ Best codons per amino acid:")
    for aa, codon in best_codons.items():
        print(f"{aa} → {codon}")

    # ========================================================
    # LOAD INPUT FASTA
    # ========================================================
    records = list(SeqIO.parse(input_fasta, "fasta"))

    optimized_records = []

    for rec in records:

        new_seq = optimize_sequence(rec.seq, codon2aa, best_codons)

        rec.seq = Seq(new_seq)   
        rec.id = f"{rec.id}_opt_{reference_choice}"
        rec.description = ""

        optimized_records.append(rec)

    # ========================================================
    # SAVE OUTPUT
    # ========================================================
    base = os.path.splitext(os.path.basename(input_fasta))[0]

    out_fasta = os.path.join(
        results_dir,
        f"optimized_{reference_choice}_{base}.fasta"
    )

    SeqIO.write(optimized_records, out_fasta, "fasta")

    print("\n✅ Done.")
    print("📄 Output:", out_fasta)

    return out_fasta