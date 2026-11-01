"""Builds lab4/results/report.md (+ grid plots, + report.pdf) from the output of ./cand.

Objective values of the methods from assignments 1-3 are taken from lab3/results/raw.csv (they do not
depend on the machine); running times are reported only for the methods measured here.

Usage: python3 report.py [--data-dir DIR] [--results DIR] [--no-pdf]
"""
import argparse
import csv
import math
import shutil
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
NEW_PREFIX = "Candidates"
LAB3_RAW = ROOT.parent / "lab3" / "results" / "raw.csv"
CHROME = ["google-chrome", "chromium", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"]

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


def load_raw(path, runs, order, skip=()):
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["method"] in skip:
                continue
            if r["method"] not in order:
                order.append(r["method"])
            runs[r["instance"]][r["method"]].append(
                (int(r["value"]), float(r["time_us"]) / 1000, int(r["evals"]),
                 float(r["ls_time_us"]) / 1000))


def load(results):
    here = defaultdict(lambda: defaultdict(list))   # inst -> method -> [(value, ms, evals, ls_ms)]
    here_order = []
    load_raw(results / "raw.csv", here, here_order)
    prev = defaultdict(lambda: defaultdict(list))
    prev_order = []
    load_raw(LAB3_RAW, prev, prev_order, skip=set(here_order))
    best = {}
    for line in (results / "best.txt").read_text().splitlines():
        inst, method, val, *tour = line.split(";")
        best[(inst, method)] = (int(val), [int(x) for x in tour])
    return here, here_order, prev, prev_order, best


def table(runs, methods, col, fmt, fmt_av):
    """Cells 'av (min – max)'; best average and best minimum of each column in bold,
    names of the methods from the current assignment underlined."""
    stats = {}
    for inst in INSTANCES:
        for m in methods:
            v = np.array([r[col] for r in runs[inst][m]], dtype=float)
            stats[inst, m] = (v.mean(), v.min(), v.max())
    win_av = {i: min(stats[i, m][0] for m in methods) for i in INSTANCES}
    win_mn = {i: min(stats[i, m][1] for m in methods) for i in INSTANCES}
    lines = ["| Method | " + " | ".join(INSTANCES) + " |", "|---|" + "---|" * len(INSTANCES)]
    for m in methods:
        name = f"<u>{m}</u>" if is_new(m) else m
        cells = []
        for i in INSTANCES:
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
    rows = math.ceil(len(methods) / 2)
    fig, axes = plt.subplots(rows, 2, figsize=(11, rows * 5 * ratio + 1.5), constrained_layout=True)
    vmin, vmax = cost.min(), cost.max()
    for ax in axes.flat[len(methods):]:
        ax.axis("off")
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


def to_pdf(text, results):
    chrome = next((c for c in CHROME if shutil.which(c) or Path(c).exists()), None)
    if chrome is None:
        print("no Chrome found, PDF skipped")
        return
    body = MarkdownIt("commonmark", {"html": True}).enable("table").render(text)
    html = results / "report.html"
    html.write_text(f"<!doctype html><meta charset='utf-8'><style>{CSS}</style>{body}")
    subprocess.run([chrome, "--headless", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
                    f"--print-to-pdf={results / 'report.pdf'}", html.as_uri()],
                   check=True, capture_output=True)
    html.unlink()
    print("wrote", results / "report.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(ROOT.parent / "data"))
    ap.add_argument("--results", default=str(ROOT / "results"))
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    results = Path(args.results)
    here, here_order, prev, prev_order, best = load(results)

    grids = []
    for inst in INSTANCES:
        data = np.loadtxt(Path(args.data_dir) / f"{inst}.csv", delimiter=";")
        xy, cost = data[:, :2], data[:, 2].astype(int)
        k = math.ceil(len(cost) / 2)
        for m in here_order:   # independent check of every best solution
            val, tour = best[inst, m]
            assert checker(tour, xy, cost, k) == val, (inst, m)
        grid_plot(inst, here_order, best, xy, cost, results / f"{inst}_best_grid.png")
        grids.append(f"### {inst}\n\n![]({inst}_best_grid.png)\n")
        for m in here_order:
            val, tour = best[inst, m]
            grids.append(f"{inst}, {m}: {val}\n\n```\n{tour}\n```\n")

    # objective table: previous methods (lab3 raw data) + methods run here
    allruns = defaultdict(dict)
    for inst in INSTANCES:
        allruns[inst].update(prev[inst])
        allruns[inst].update(here[inst])
    obj_methods = prev_order + here_order

    ms = lambda x: f"{x:.3f}"
    num = lambda x: f"{x:.0f}"
    text = (ROOT / "report_template.md").read_text()
    text = (text.replace("{{OBJ_TABLE}}", table(allruns, obj_methods, 0, num, lambda x: f"{x:.1f}"))
                .replace("{{TIME_TABLE}}", table(here, here_order, 1, ms, ms))
                .replace("{{LSTIME_TABLE}}", table(here, here_order, 3, ms, ms))
                .replace("{{EVALS_TABLE}}", table(here, here_order, 2, num, num))
                .replace("{{GRIDS}}", "\n".join(grids)))
    (results / "report.md").write_text(text)
    print("wrote", results / "report.md")
    if not args.no_pdf:
        to_pdf(text, results)


if __name__ == "__main__":
    main()
