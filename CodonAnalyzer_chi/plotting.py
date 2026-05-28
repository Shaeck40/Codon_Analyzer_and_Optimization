#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt


def make_scatter(df, xcol, ycol, labelcol, title, outdir):

    plt.figure(figsize=(8, 6))

    x = df[xcol]
    y = df[ycol]

    if labelcol == "codon":
        colors = df["codon"].astype("category").cat.codes
        plt.scatter(x, y, c=colors, cmap="tab20", alpha=0.8)
    else:
        plt.scatter(x, y, color="blue", alpha=0.8)

    if len(df) > 1:
        m, b = np.polyfit(x, y, 1)
        plt.plot(x, m*x + b, "--", color="black")

    for _, r in df.iterrows():
        plt.text(r[xcol], r[ycol], r[labelcol], fontsize=6)

    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.title(title)

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{title}.png")

    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()

    return path