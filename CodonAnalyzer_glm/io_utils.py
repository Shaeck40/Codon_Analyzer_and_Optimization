#!/usr/bin/env python3

import os
import sys
import pandas as pd
from Bio import SeqIO


def resource_path(path):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, path)
    return os.path.join(os.path.abspath("."), path)


def load_genetic_code():
    path = os.path.join(os.path.abspath("."), "genetic_code.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing: {path}")
    return dict(zip(pd.read_csv(path)["codon"], pd.read_csv(path)["amino_acid"]))


def parse_fasta(path):
    path = resource_path(path)

    records = []
    for rec in SeqIO.parse(path, "fasta"):
        parts = rec.description.split()

        records.append({
            "cds_id": parts[0],
            "ProteinName": parts[1] if len(parts) > 1 else "",
            "GeneNames": parts[2] if len(parts) > 2 else "",
            "sequence": str(rec.seq)
        })

    return pd.DataFrame(records)