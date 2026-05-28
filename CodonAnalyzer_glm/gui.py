#!/usr/bin/env python3

import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import threading
import traceback

LAST_RESULTS_DIR = None

# =========================================================
# SAFE IMPORT
# =========================================================
try:
    from analysis import run_pipeline
    from optimization import optimize_fasta_pipeline
except Exception as e:
    run_pipeline = None
    optimize_fasta_pipeline = None
    IMPORT_ERROR = str(e)
else:
    IMPORT_ERROR = None


# =========================================================
# THREAD SAFE LOGGER
# =========================================================
def log_write(widget, msg):
    widget.insert(tk.END, msg + "\n")
    widget.see(tk.END)


# =========================================================
# ANALYSIS
# =========================================================
def run_analysis(excel_path, ref_choice, log_widget):
    global LAST_RESULTS_DIR

    if IMPORT_ERROR:
        messagebox.showerror("Import error", IMPORT_ERROR)
        return

    if not excel_path:
        messagebox.showwarning("No file", "Choose input Excel file.")
        return

    try:
        log_write(log_widget, f"Running analysis ({ref_choice})...")

        outdir = run_pipeline(excel_path, ref_choice)
        LAST_RESULTS_DIR = outdir

        log_write(log_widget, "Analysis completed.")
        log_write(log_widget, f"Output folder: {outdir}")

    except Exception:
        log_write(log_widget, traceback.format_exc())


# =========================================================
# OPTIMIZATION
# =========================================================
def run_optimization(fasta_path, optimize_choice, log_widget):
    global LAST_RESULTS_DIR

    if IMPORT_ERROR:
        messagebox.showerror("Import error", IMPORT_ERROR)
        return

    if not fasta_path:
        messagebox.showwarning("No FASTA", "Choose input FASTA file.")
        return

    try:
        log_write(log_widget, f"Running optimization ({optimize_choice})...")

        results_dir = LAST_RESULTS_DIR

        out = optimize_fasta_pipeline(
            input_fasta=fasta_path,
            reference_choice=optimize_choice,
            results_dir=results_dir
        )

        log_write(log_widget, "Optimization completed.")
        log_write(log_widget, f"Output FASTA: {out}")

    except Exception:
        log_write(log_widget, traceback.format_exc())


# =========================================================
# THREAD WRAPPER
# =========================================================
def start(func, *args):
    threading.Thread(target=func, args=args, daemon=True).start()


# =========================================================
# FILE SELECT
# =========================================================
def select_file(entry, filetypes):
    f = filedialog.askopenfilename(filetypes=filetypes)
    if f:
        entry.delete(0, tk.END)
        entry.insert(0, f)


# =========================================================
# GUI
# =========================================================
def main():

    root = tk.Tk()
    root.title("Codon Usage Analyzer & Optimizer")
    root.geometry("750x650")

    tk.Label(root, text="Codon Usage Analyzer & Optimizer",
             font=("Arial", 16, "bold")).pack(pady=10)

    # ---------------- Excel ----------------
    tk.Label(root, text="Select target Excel file:").pack()
    excel_entry = tk.Entry(root, width=60)
    excel_entry.pack(pady=5)

    tk.Button(root, text="Browse",
              command=lambda: select_file(
                  excel_entry,
                  [("Excel files", "*.xlsx *.xls")]
              )).pack()

    # ---------------- FASTA ----------------
    tk.Label(root, text="Select FASTA file for optimization (optional):").pack(pady=10)
    fasta_entry = tk.Entry(root, width=60)
    fasta_entry.pack(pady=5)

    tk.Button(root, text="Browse",
              command=lambda: select_file(
                  fasta_entry,
                  [("FASTA files", "*.fasta *.fa")]
              )).pack()

    # ---------------- ANALYSIS OPTIONS ----------------
    ref_choice = tk.StringVar(value="both")

    tk.Label(root, text="Compare with:", font=("Arial", 12)).pack(pady=10)

    tk.Radiobutton(root, text="OpenPichia",
                   variable=ref_choice, value="open").pack()

    tk.Radiobutton(root, text="NCBI",
                   variable=ref_choice, value="ncbi").pack()

    tk.Radiobutton(root, text="Both",
                   variable=ref_choice, value="both").pack()

    # ---------------- OPTIMIZATION ----------------
    optimize_choice = tk.StringVar(value="target")

    tk.Label(root, text="Optimize using:", font=("Arial", 12)).pack(pady=10)

    tk.Radiobutton(root, text="Target",
                   variable=optimize_choice, value="target").pack()

    tk.Radiobutton(root, text="OpenPichia",
                   variable=optimize_choice, value="open").pack()

    tk.Radiobutton(root, text="NCBI",
                   variable=optimize_choice, value="ncbi").pack()

    # ---------------- LOG ----------------
    tk.Label(root, text="Log output:").pack(pady=10)

    log = scrolledtext.ScrolledText(root, width=90, height=5)
    log.pack()

    # ---------------- BUTTONS ----------------
    frame = tk.Frame(root)
    frame.pack(pady=15)

    tk.Button(
        frame,
        text="Run Analysis",
        font=("Arial", 12, "bold"),
        width=15,
        command=lambda: start(
            run_analysis,
            excel_entry.get(),
            ref_choice.get(),
            log
        )
    ).pack(side="left", padx=10)

    tk.Button(
        frame,
        text="Run Optimization",
        font=("Arial", 12, "bold"),
        width=15,
        command=lambda: start(
            run_optimization,
            fasta_entry.get(),        # only FASTA
            optimize_choice.get(),
            log
        )
    ).pack(side="left", padx=10)

    root.mainloop()


if __name__ == "__main__":
    main()
