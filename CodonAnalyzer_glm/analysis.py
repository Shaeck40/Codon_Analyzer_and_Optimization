#!/usr/bin/env python3

import os
from datetime import datetime
from collections import Counter
import shutil

import numpy as np
import pandas as pd
from Bio import SeqIO

import statsmodels.api as sm
import statsmodels.formula.api as smf
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
        total.update([seq[i:i+3] for i in range(0, len(seq)-2, 3)])
    tot = sum(total.values())
    return pd.DataFrame([
        {"codon": c, "count": v, "perc": (v/tot*100 if tot > 0 else 0)}
        for c, v in total.items()
    ])


# ============================================================
# MATCH TARGETS
# ============================================================
def match_targets(df, targets):
    df = df.copy()
    df["GeneNames"] = df.get("GeneNames","").astype(str).str.lower()
    df["ProteinName"] = df.get("ProteinName","").astype(str).str.lower()
    df["cds_id"] = df["cds_id"].astype(str).str.lower()

    targets = set(str(t).lower() for t in targets)

    return df[
        df["GeneNames"].isin(targets)
        | df["ProteinName"].isin(targets)
        | df["cds_id"].isin(targets)
    ]


# ============================================================
# AA GLM
# ============================================================
def glm_amino_acid_codon_distribution(glm_df, samples):

    rows, pvals = [], []
    refs = [s for s in samples if s != "target"]

    for ref in refs:
        sub = glm_df[glm_df["sample"].isin(["target", ref])]

        for aa, group in sub.groupby("AA"):

            if group["codon"].nunique() < 2:
                continue

            group = group.copy()
            group["sample"] = pd.Categorical(group["sample"], ["target", ref])
            group["offset"] = np.log(group.groupby("sample")["count"].transform("sum") + 1e-9)

            try:
                model = sm.GLM.from_formula(
                    "count ~ C(codon) + C(sample) + C(codon):C(sample)",
                    data=group,
                    family=sm.families.Poisson(),
                    offset=group["offset"]
                ).fit()
            except:
                continue

            terms = [t for t in model.pvalues.index if f":C(sample)[T.{ref}]" in t]
            if not terms:
                continue

            p_global = np.min([model.pvalues[t] for t in terms])

            rows.append({
                "AA": aa,
                "sample": ref,
                "p_raw": p_global,
                "n_codons": group["codon"].nunique()
            })
            pvals.append(p_global)

    if not rows:
        return pd.DataFrame(columns=["AA","sample","p_raw","p_adj","n_codons"])

    _, p_adj, _, _ = multipletests(pvals, method="fdr_bh")

    for r, adj in zip(rows, p_adj):
        r["p_adj"] = adj

    return pd.DataFrame(rows)


# ============================================================
# CODON GLM
# ============================================================
def glm_codons_poisson(glm_df, samples):

    rows, pvals = [], []
    refs = [s for s in samples if s != "target"]

    for ref in refs:
        sub = glm_df[glm_df["sample"].isin(["target", ref])]

        for aa, group in sub.groupby("AA"):

            if group["codon"].nunique() < 2:
                continue

            group = group.copy()
            group["sample"] = pd.Categorical(group["sample"], ["target", ref])
            group["offset"] = np.log(group.groupby("sample")["count"].transform("sum")+1e-9)

            try:
                model = smf.glm(
                    "count ~ C(codon) + C(sample) + C(codon):C(sample)",
                    data=group,
                    family=sm.families.Poisson(),
                    offset=group["offset"]
                ).fit()
            except:
                continue

            for term, pv in model.pvalues.items():
                if "codon" in term and f":C(sample)[T.{ref}]" in term:
                    cod = term.split("[T.")[1].split("]")[0]

                    rows.append({
                        "AA": aa,
                        "codon": cod,
                        "sample": ref,
                        "coef": model.params[term],
                        "p_raw": pv
                    })
                    pvals.append(pv)

    if not rows:
        return pd.DataFrame(columns=["AA","codon","sample","coef","p_raw","p_adj"])

    _, p_adj, _, _ = multipletests(pvals, method="fdr_by")

    for r, adj in zip(rows, p_adj):
        r["p_adj"] = adj

    return pd.DataFrame(rows)


# ============================================================
# OPTIMIZATION
# ============================================================
def optimize_sequence(seq, best_codons, codon2aa):
    new_seq = ""
    for i in range(0, len(seq)-2, 3):
        cod = seq[i:i+3]
        aa = codon2aa.get(cod)
        new_seq += best_codons.get(aa, cod)
    return new_seq


def optimize_fasta(input_fasta, output_fasta, best_codons, codon2aa):
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            new_seq = optimize_sequence(str(rec.seq), best_codons, codon2aa)
            out.write(f">{rec.id}_optimized\n{new_seq}\n")


# ============================================================
# MAIN PIPELINE
# ============================================================
def run_pipeline(targets_excel, ref_choice="both", fasta=None):

    ts = datetime.now().strftime("%y%m%d_%H%M")
    results_dir = f"results_{ts}"
    os.makedirs(results_dir, exist_ok=True)

    input_dir = os.path.join(results_dir, "input")
    os.makedirs(input_dir, exist_ok=True)
    shutil.copy(targets_excel, input_dir)

    name = os.path.splitext(os.path.basename(targets_excel))[0]

    codon2aa = load_genetic_code()
    cds_df = parse_fasta("op_cds_chr1-4.fasta")
    ncbi_seqs = [str(r.seq) for r in SeqIO.parse("NCBI_cds.fasta","fasta")]

    open_u = calculate_codon_usage(cds_df["sequence"])
    ncbi_u = calculate_codon_usage(ncbi_seqs)

    excel = pd.ExcelFile(targets_excel)

    for sheet in excel.sheet_names:

        df = excel.parse(sheet)
        targets = df.iloc[:,0]

        matched = match_targets(cds_df, targets)
        if len(matched)==0:
            continue

        target = calculate_codon_usage(matched["sequence"])
        codon = target.rename(columns={"count":"count_target","perc":"perc_target"})

        if ref_choice in ("open","both"):
            codon = codon.merge(open_u.rename(columns={"count":"count_open","perc":"perc_open"}), on="codon", how="outer")

        if ref_choice in ("ncbi","both"):
            codon = codon.merge(ncbi_u.rename(columns={"count":"count_ncbi","perc":"perc_ncbi"}), on="codon", how="outer")

        codon = codon.fillna(0)
        codon["AA"] = codon["codon"].map(codon2aa)
        codon = codon.sort_values(["AA","codon"])

        # SYN %
        for p,c in [("target","count_target"),("open","count_open"),("ncbi","count_ncbi")]:
            if c in codon.columns:
                tot = codon.groupby("AA")[c].transform("sum")
                codon[f"perc_{p}_syn"] = np.where(tot>0, codon[c]/tot*100,0)

        # ===== PLOTS =====
        if "perc_open" in codon.columns:
            make_scatter(codon,"perc_target","perc_open","codon",f"{sheet}_codon_open",results_dir)

        if "perc_ncbi" in codon.columns:
            make_scatter(codon,"perc_target","perc_ncbi","codon",f"{sheet}_codon_ncbi",results_dir)

        # ===== GLM =====
        samples=["target"]
        if "count_open" in codon.columns: samples.append("open")
        if "count_ncbi" in codon.columns: samples.append("ncbi")

        glm_df = pd.DataFrame([
            {"AA":r["AA"],"codon":r["codon"],"sample":s,"count":r.get(f"count_{s}",0)}
            for _,r in codon.iterrows() for s in samples
        ])

        glm_cod = glm_codons_poisson(glm_df, samples)
        glm_aa = glm_amino_acid_codon_distribution(glm_df, samples)

        out = os.path.join(results_dir,f"{name}_{sheet}.xlsx")

        with pd.ExcelWriter(out) as writer:

            codon.to_excel(writer, sheet_name="codons", index=False)

            # CLEAN OPEN
            if "open" in samples:
                merged = codon.merge(glm_cod[glm_cod["sample"]=="open"], on=["AA","codon"], how="left")
                df_open = pd.DataFrame({
                    "codon": merged["codon"],
                    "AA": merged["AA"],
                    "target_count": merged["count_target"],
                    "open_count": merged["count_open"],
                    "target_syn%": merged["perc_target_syn"],
                    "open_syn%": merged["perc_open_syn"],
                    "coef": merged["coef"],
                    "p_raw": merged["p_raw"],
                    "p_adj": merged["p_adj"]
                }).sort_values(["AA","codon"])

                df_open.to_excel(writer, sheet_name="codon_vs_open", index=False)

            # CLEAN NCBI
            if "ncbi" in samples:
                merged = codon.merge(glm_cod[glm_cod["sample"]=="ncbi"], on=["AA","codon"], how="left")
                df_ncbi = pd.DataFrame({
                    "codon": merged["codon"],
                    "AA": merged["AA"],
                    "target_count": merged["count_target"],
                    "ncbi_count": merged["count_ncbi"],
                    "target_syn%": merged["perc_target_syn"],
                    "ncbi_syn%": merged["perc_ncbi_syn"],
                    "coef": merged["coef"],
                    "p_raw": merged["p_raw"],
                    "p_adj": merged["p_adj"]
                }).sort_values(["AA","codon"])

                df_ncbi.to_excel(writer, sheet_name="codon_vs_ncbi", index=False)

            # AA
            if not glm_aa.empty:
                if "open" in samples:
                    glm_aa[glm_aa["sample"]=="open"].to_excel(writer, sheet_name="aa_vs_open", index=False)
                if "ncbi" in samples:
                    glm_aa[glm_aa["sample"]=="ncbi"].to_excel(writer, sheet_name="aa_vs_ncbi", index=False)

        # ===== OPTIMIZATION =====
        if fasta:
            best = codon.sort_values("perc_target").groupby("AA").last()["codon"].to_dict()
            out_fasta = os.path.join(results_dir,f"{name}_{sheet}_optimized.fasta")
            optimize_fasta(fasta,out_fasta,best,codon2aa)

    print("Done ✅")
    return results_dir



