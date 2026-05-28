#!/usr/bin/env python3

import os
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO

from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

from io_utils import load_genetic_code, parse_fasta
from plotting import make_scatter


# ============================================================
# CODON USAGE
# ============================================================
def calculate_codon_usage(seqs):

    total = Counter()

    for seq in seqs:

        seq = str(seq).upper()

        total.update([
            seq[i:i+3]
            for i in range(0, len(seq) - 2, 3)
        ])

    tot = sum(total.values())

    return pd.DataFrame([
        {
            "codon": c,
            "count": v,
            "perc": (v / tot * 100) if tot > 0 else 0
        }
        for c, v in total.items()
    ])


# ============================================================
# MATCH TARGETS
# ============================================================
def match_targets(df, targets):

    df = df.copy()

    df["GeneNames"] = df.get(
        "GeneNames", ""
    ).astype(str).str.lower()

    df["ProteinName"] = df.get(
        "ProteinName", ""
    ).astype(str).str.lower()

    df["cds_id"] = df["cds_id"].astype(str).str.lower()

    targets = set(str(t).lower() for t in targets)

    return df[
        df["GeneNames"].isin(targets)
        | df["ProteinName"].isin(targets)
        | df["cds_id"].isin(targets)
    ]


# ============================================================
# RSCU
# ============================================================
def calculate_rscu(df, count_col):

    df = df.copy()

    aa_totals = df.groupby("AA")[count_col].transform("sum")
    n_codons = df.groupby("AA")["codon"].transform("count")

    expected = aa_totals / n_codons

    df["RSCU"] = np.where(
        expected > 0,
        df[count_col] / expected,
        np.nan
    )

    return df["RSCU"]


# ============================================================
# CHI-SQUARE CODON BIAS
# ============================================================
def chi_square_codons(codon_df, reference):

    rows = []
    pvals = []

    target_col = "count_target"
    ref_col = f"count_{reference}"

    for aa, group in codon_df.groupby("AA"):

        if len(group) < 2:
            continue

        target_counts = group[target_col].values
        ref_counts = group[ref_col].values

        # Skip empty groups
        if target_counts.sum() == 0 or ref_counts.sum() == 0:
            continue

        contingency = np.array([
            target_counts,
            ref_counts
        ])

        try:
            chi2, p, dof, expected = chi2_contingency(contingency)

        except ValueError:
            continue

        pvals.append(p)

        # log2 enrichment per codon
        for _, row in group.iterrows():

            target_rscu = row["RSCU_target"]
            ref_rscu = row[f"RSCU_{reference}"]

            if pd.isna(target_rscu) or pd.isna(ref_rscu):
                log2fc = np.nan
            else:
                log2fc = np.log2(
                    (target_rscu + 1e-9) /
                    (ref_rscu + 1e-9)
                )

            rows.append({
                "AA": aa,
                "codon": row["codon"],
                "reference": reference,
                "chi2_p_raw": p,
                "chi2_p_adj": np.nan,
                "target_count": row[target_col],
                "ref_count": row[ref_col],
                "target_perc_syn": row["perc_target_syn"],
                "ref_perc_syn": row[f"perc_{reference}_syn"],
                "RSCU_target": target_rscu,
                f"RSCU_{reference}": ref_rscu,
                "log2_RSCU_ratio": log2fc
            })

    if not rows:
        return pd.DataFrame()

    _, p_adj, _, _ = multipletests(
        pvals,
        method="fdr_bh"
    )

    # Assign same corrected p-value to all codons
    # within same amino acid
    aa_order = []

    for aa, group in codon_df.groupby("AA"):
        if len(group) >= 2:
            aa_order.append(aa)

    aa_to_adj = {
        aa: adj
        for aa, adj in zip(aa_order, p_adj)
    }

    for r in rows:
        r["chi2_p_adj"] = aa_to_adj.get(r["AA"], np.nan)

    return pd.DataFrame(rows)


# ============================================================
# SAFE SHEET
# ============================================================
def safe_sheet(df, name):

    if df is None or df.empty:

        return pd.DataFrame({
            "info": [f"No significant results for {name}"]
        })

    return df


# ============================================================
# MAIN PIPELINE
# ============================================================
def run_pipeline(targets_excel, ref_choice="both"):

    import shutil

    timestamp = datetime.now().strftime("%y%m%d_%H%M")

    results_dir = f"results_{timestamp}"

    os.makedirs(results_dir, exist_ok=True)

    input_dir = os.path.join(results_dir, "input")
    os.makedirs(input_dir, exist_ok=True)

    shutil.copy(targets_excel, input_dir)

    input_name = os.path.splitext(
        os.path.basename(targets_excel)
    )[0]

    codon2aa = load_genetic_code()

    cds_df = parse_fasta("op_cds_chr1-4.fasta")

    ncbi_seqs = [
        str(r.seq)
        for r in SeqIO.parse("NCBI_cds.fasta", "fasta")
    ]

    openpichia = calculate_codon_usage(
        cds_df["sequence"]
    )

    ncbi = calculate_codon_usage(
        ncbi_seqs
    )

    excel = pd.ExcelFile(targets_excel)

    for sheet in excel.sheet_names:

        df = excel.parse(sheet)

        targets = df.iloc[:, 0].astype(str)

        matched = match_targets(cds_df, targets)

        if len(matched) == 0:
            continue

        target = calculate_codon_usage(
            matched["sequence"]
        )

        codon = target.rename(columns={
            "count": "count_target",
            "perc": "perc_target"
        })

        # =====================================================
        # MERGE REFERENCES
        # =====================================================
        if ref_choice in ("open", "both"):

            codon = codon.merge(
                openpichia.rename(columns={
                    "count": "count_open",
                    "perc": "perc_open"
                }),
                on="codon",
                how="outer"
            )

        if ref_choice in ("ncbi", "both"):

            codon = codon.merge(
                ncbi.rename(columns={
                    "count": "count_ncbi",
                    "perc": "perc_ncbi"
                }),
                on="codon",
                how="outer"
            )

        codon = codon.fillna(0)

        codon["AA"] = codon["codon"].map(codon2aa)

        codon = codon.dropna(subset=["AA"])

        # =====================================================
        # SYNONYMOUS PERCENTAGES
        # =====================================================
        for prefix, count_col in [
            ("target", "count_target"),
            ("open", "count_open"),
            ("ncbi", "count_ncbi"),
        ]:

            if count_col not in codon.columns:
                continue

            aa_totals = codon.groupby("AA")[
                count_col
            ].transform("sum")

            codon[f"perc_{prefix}_syn"] = np.where(
                aa_totals > 0,
                codon[count_col] / aa_totals * 100,
                0
            )

        # =====================================================
        # RSCU
        # =====================================================
        codon["RSCU_target"] = calculate_rscu(
            codon,
            "count_target"
        )

        if "count_open" in codon.columns:

            codon["RSCU_open"] = calculate_rscu(
                codon,
                "count_open"
            )

        if "count_ncbi" in codon.columns:

            codon["RSCU_ncbi"] = calculate_rscu(
                codon,
                "count_ncbi"
            )

        # =====================================================
        # AMINO ACID SUMMARY
        # =====================================================
        aa = codon.groupby("AA")[[
            "count_target"
        ]].sum().reset_index()

        aa = aa.rename(columns={
            "AA": "amino_acid",
            "count_target": "Target"
        })

        if "count_open" in codon.columns:

            aa["OpenPichia"] = codon.groupby(
                "AA"
            )["count_open"].sum().values

        if "count_ncbi" in codon.columns:

            aa["NCBI"] = codon.groupby(
                "AA"
            )["count_ncbi"].sum().values

        # =====================================================
        # PLOTS (UNCHANGED)
        # =====================================================
        rename_codon = {
            "perc_target": "Target_perc"
        }

        if "perc_open" in codon.columns:
            rename_codon["perc_open"] = "OpenPichia_perc"

        if "perc_ncbi" in codon.columns:
            rename_codon["perc_ncbi"] = "NCBI_perc"

        codon_plot = codon.rename(columns=rename_codon)

        base_title = f"{input_name}_{sheet.replace(' ', '_')}"

        if "OpenPichia_perc" in codon_plot.columns:

            make_scatter(
                codon_plot,
                "Target_perc",
                "OpenPichia_perc",
                "codon",
                f"{base_title}_codon_Target_vs_OpenPichia",
                results_dir
            )

        if "NCBI_perc" in codon_plot.columns:

            make_scatter(
                codon_plot,
                "Target_perc",
                "NCBI_perc",
                "codon",
                f"{base_title}_codon_Target_vs_NCBI",
                results_dir
            )

        if "OpenPichia" in aa.columns:

            make_scatter(
                aa,
                "Target",
                "OpenPichia",
                "amino_acid",
                f"{base_title}_aa_Target_vs_OpenPichia",
                results_dir
            )

        if "NCBI" in aa.columns:

            make_scatter(
                aa,
                "Target",
                "NCBI",
                "amino_acid",
                f"{base_title}_aa_Target_vs_NCBI",
                results_dir
            )

        # =====================================================
        # CHI-SQUARE ANALYSIS
        # =====================================================
        chi_open = pd.DataFrame()
        chi_ncbi = pd.DataFrame()

        if "count_open" in codon.columns:

            chi_open = chi_square_codons(
                codon,
                "open"
            )

        if "count_ncbi" in codon.columns:

            chi_ncbi = chi_square_codons(
                codon,
                "ncbi"
            )

        # =====================================================
        # WRITE OUTPUT
        # =====================================================
        safe_sheet_name = sheet.replace(" ", "_")

        out_file = os.path.join(
            results_dir,
            f"{input_name}_{safe_sheet_name}_results.xlsx"
        )

        with pd.ExcelWriter(
            out_file,
            engine="openpyxl"
        ) as writer:

            codon.to_excel(
                writer,
                sheet_name="codons",
                index=False
            )

            aa.to_excel(
                writer,
                sheet_name="amino_acids",
                index=False
            )

            safe_sheet(
                chi_open,
                "chi_square_open"
            ).to_excel(
                writer,
                sheet_name="codon_bias_vs_open",
                index=False
            )

            safe_sheet(
                chi_ncbi,
                "chi_square_ncbi"
            ).to_excel(
                writer,
                sheet_name="codon_bias_vs_ncbi",
                index=False
            )

    print("Done.")

    return results_dir