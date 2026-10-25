"""Builds lab3/results/report.md (+ grid plots, + report.pdf) from the output of ./ls.

Usage: python3 report.py [--data-dir DIR] [--results DIR] [--no-pdf]
"""
import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parent
INSTANCES = ("TSPA", "TSPB")
NEW_PREFIX = ("Steepest", "Greedy,")

CSS = """
@page { size: A4; margin: 15mm; }
body { font-family: sans-serif; font-size: 10.5pt; line-height: 1.4; }
table { border-collapse: collapse; margin: 8px 0; font-size: 9pt; }
th, td { border: 1px solid #999; padding: 3px 8px; text-align: left; }
img { max-width: 100%; }
pre, code { font-size: 8.5pt; white-space: pre-wrap; word-break: break-all; }
h1 { font-size: 18pt; } h2 { font-size: 14pt; } h3 { font-size: 11.5pt; }
"""


def is_new(method):
    return method.startswith(NEW_PREFIX)


def load(results):
    runs = defaultdict(lambda: defaultdict(list))  # inst -> method -> [(value, ms, evals, ls_ms)]
    order = []
    with open(results / "raw.csv") as f:
        for r in csv.DictReader(f):
            if r["method"] not in order:
                order.append(r["method"])
            runs[r["instance"]][r["method"]].append(
                (int(r["value"]), float(r["time_us"]) / 1000, int(r["evals"]),
                 float(r["ls_time_us"]) / 1000))
    best = {}
    for line in (results / "best.txt").read_text().splitlines():
        inst, method, val, *tour = line.split(";")
        best[(inst, method)] = (int(val), [int(x) for x in tour])
    return runs, order, best


def table(runs, order, col, fmt, fmt_av, title_cols=INSTANCES, methods=None):
    """Cells 'av (min – max)'; best average and best minimum of each column in bold,
    names of the methods from the current assignment underlined."""
    methods = methods or order
    stats = {}
    for inst in title_cols:
        for m in methods:
            v = np.array([r[col] for r in runs[inst][m]], dtype=float)
            stats[inst, m] = (v.mean(), v.min(), v.max())
    win_av = {i: min(stats[i, m][0] for m in methods) for i in title_cols}
    win_mn = {i: min(stats[i, m][1] for m in methods) for i in title_cols}
    lines = ["| Method | " + " | ".join(title_cols) + " |", "|---|" + "---|" * len(title_cols)]
    for m in methods:
        name = f"<u>{m}</u>" if is_new(m) else m
        cells = []
        for i in title_cols:
            av, mn, mx = stats[i, m]
            a, lo, hi = fmt_av(av), fmt(mn), fmt(mx)
            if av == win_av[i]:
                a = f"**{a}**"
            if mn == win_mn[i]:
                lo = f"**{lo}**"
            cells.append(f"{a} ({lo} – {hi})")
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def grid_plot(inst, methods, best, xy, cost, out):
    ratio = np.ptp(xy[:, 1]) / np.ptp(xy[:, 0])
    fig, axes = plt.subplots(4, 2, figsize=(11, 4 * 5 * ratio + 1.5), constrained_layout=True)
    vmin, vmax = cost.min(), cost.max()
    for ax, m in zip(axes.flat, methods):
        val, tour = best[inst, m]
        t = tour + [tour[0]]
        ax.plot(xy[t, 0], xy[t, 1], "-", color="0.4", lw=0.9, zorder=1)
        sel = np.zeros(len(cost), dtype=bool)
        sel[tour] = True
        ax.scatter(xy[~sel, 0], xy[~sel, 1], c=cost[~sel], cmap="viridis", s=10, alpha=0.35,
                   vmin=vmin, vmax=vmax, zorder=2)
        sc = ax.scatter(xy[sel, 0], xy[sel, 1], c=cost[sel], cmap="viridis", s=30,
                        edgecolors="black", linewidths=0.5, vmin=vmin, vmax=vmax, zorder=3)
        ax.set_title(f"{m}\nobjective = {val}", fontsize=11)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(sc, ax=axes, label="node cost", shrink=0.8, pad=0.01)
    fig.suptitle(f"{inst} - best solutions found", fontsize=14)
    fig.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)


def checker(tour, xy, cost, k):
    assert len(tour) == k and len(set(tour)) == k
    total = 0
    for i in range(k):
        a, b = tour[i], tour[(i + 1) % k]
        dx, dy = xy[a] - xy[b]
        total += int(math.sqrt(dx * dx + dy * dy) + 0.5) + int(cost[a])
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(ROOT.parent / "lab1"))
    ap.add_argument("--results", default=str(ROOT / "results"))
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    results = Path(args.results)
    runs, order, best = load(results)
    new = [m for m in order if is_new(m)]

    sections = []
    grids = []
    for inst in INSTANCES:
        data = np.loadtxt(Path(args.data_dir) / f"{inst}.csv", delimiter=";")
        xy, cost = data[:, :2], data[:, 2].astype(int)
        k = math.ceil(len(cost) / 2)
        for m in order:   # independent check of every best solution
            val, tour = best[inst, m]
            assert checker(tour, xy, cost, k) == val, (inst, m)
        grid_plot(inst, new, best, xy, cost, results / f"{inst}_best_grid.png")
        grids.append(f"### {inst}\n\n![]({inst}_best_grid.png)\n")
        for m in new:
            val, tour = best[inst, m]
            grids.append(f"{inst}, {m}: {val}\n\n```\n{tour}\n```\n")

    ms = lambda x: f"{x:.3f}"
    ls_methods = [m for m in order if is_new(m)]
    text = (ROOT / "report_template.md").read_text()
    text = (text.replace("{{OBJ_TABLE}}", table(runs, order, 0, lambda x: f"{x:.0f}", lambda x: f"{x:.1f}"))
                .replace("{{TIME_TABLE}}", table(runs, order, 1, ms, ms))
                .replace("{{LSTIME_TABLE}}", table(runs, order, 3, ms, ms, methods=ls_methods))
                .replace("{{EVALS_TABLE}}", table(runs, order, 2, lambda x: f"{x:.0f}", lambda x: f"{x:.0f}", methods=ls_methods))
                .replace("{{GRIDS}}", "\n".join(grids)))
    (results / "report.md").write_text(text)
    print("wrote", results / "report.md")

    if not args.no_pdf:
        body = MarkdownIt("commonmark", {"html": True}).enable("table").render(text)
        html = results / "report.html"
        html.write_text(f"<!doctype html><meta charset='utf-8'><style>{CSS}</style>{body}")
        subprocess.run(["google-chrome", "--headless", "--disable-gpu", "--no-sandbox",
                        "--no-pdf-header-footer", f"--print-to-pdf={results / 'report.pdf'}",
                        html.as_uri()], check=True, capture_output=True)
        html.unlink()
        print("wrote", results / "report.pdf")


if __name__ == "__main__":
    main()
