"""Shared figure style and small helpers. Importing this module writes no files.

Only presentation and formatting live here. Group names, colours, comparisons and
every other experiment-specific choice belong in the caller's config, so that this
module carries no assumption about what is being measured.
"""
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK, GREY = "#252525", "#606060"
MM = 1 / 25.4


def apply_style():
    """Journal-neutral defaults: sans-serif, thin rules, editable vector text."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.linewidth": 0.65, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.top": False, "ytick.right": False, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.color": INK, "ytick.color": INK, "legend.frameon": False, "text.usetex": False,
        # Keep text as text in SVG and PDF so figures stay editable downstream.
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "figure.dpi": 150,
    })


def describe(values):
    """Mean, sample SD (ddof=1), median and the full SD interval, including a negative lower end."""
    m, sd = float(values.mean()), float(values.std(ddof=1))
    return {"n": len(values), "mean": m, "sd": sd, "median": float(np.median(values)),
            "mean_minus_sd": m - sd, "mean_plus_sd": m + sd, "zeros": int((values == 0).sum())}


def panel_heading(ax, letter, title):
    """Bold panel letter plus a plain title, both anchored outside the axes."""
    ax.set_label(letter)
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-17, 39),
                textcoords="offset points", fontsize=8, fontweight="bold",
                va="bottom", annotation_clip=False)
    ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 39),
                textcoords="offset points", fontsize=8, va="bottom", annotation_clip=False)


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
