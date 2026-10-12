import math
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
N_RANDOM = 200
SEED = 42


# data
def load_instance(name):
    data = np.loadtxt(ROOT / f"{name}.csv", delimiter=";", dtype=np.int64)
    return data[:, :2].astype(float), data[:, 2]


def distance_matrix(xy):
    diff = xy[:, None, :] - xy[None, :, :]
    return np.floor(np.sqrt((diff ** 2).sum(-1)) + 0.5).astype(np.int64)


def objective(tour, D, cost):
    t = np.asarray(tour)
    return int(D[t, np.roll(t, -1)].sum() + cost[t].sum())


# methods
def random_solution(D, cost, k, rng):
    return rng.sample(range(len(cost)), k)


def nn_end(D, cost, k, start):
    path = [start]
    free = np.ones(len(cost), dtype=bool)
    free[start] = False
    while len(path) < k:
        delta = D[path[-1]] + cost
        delta[~free] = np.iinfo(np.int64).max
        j = int(delta.argmin())
        path.append(j)
        free[j] = False
    return path


def nn_any(D, cost, k, start):
    path = [start]
    free = np.ones(len(cost), dtype=bool)
    free[start] = False
    while len(path) < k:
        cand = np.flatnonzero(free)
        p = np.array(path)
        best_front = D[p[0], cand] + cost[cand]
        best_back = D[p[-1], cand] + cost[cand]
        options = [best_front[None, :]]
        if len(p) > 1:
            a, b = p[:-1], p[1:]
            inner = D[a][:, cand] + D[b][:, cand] - D[a, b][:, None] + cost[cand][None, :]
            options.append(inner)
        options.append(best_back[None, :])
        delta = np.vstack(options)  
        pos, ci = np.unravel_index(delta.argmin(), delta.shape)
        j = int(cand[ci])
        path.insert(int(pos), j)
        free[j] = False
    return path


def greedy_cycle(D, cost, k, start):
    free = np.ones(len(cost), dtype=bool)
    free[start] = False
    delta = D[start] + cost
    delta[~free] = np.iinfo(np.int64).max
    second = int(delta.argmin())
    cycle = [start, second]
    free[second] = False
    while len(cycle) < k:
        cand = np.flatnonzero(free)
        a = np.array(cycle)
        b = np.roll(a, -1)
        delta = D[a][:, cand] + D[b][:, cand] - D[a, b][:, None] + cost[cand][None, :]
        e, ci = np.unravel_index(delta.argmin(), delta.shape)
        j = int(cand[ci])
        cycle.insert(int(e) + 1, j)
        free[j] = False
    return cycle


METHODS = {
    "Random": None,
    "NN (end only)": nn_end,
    "NN (any position)": nn_any,
    "Greedy cycle": greedy_cycle,
}


# checker
def check(tour, xy, cost, k):
    assert len(tour) == k and len(set(tour)) == k, "wrong number of / repeated nodes"
    total = 0
    for i in range(k):
        a, b = tour[i], tour[(i + 1) % k]
        dx, dy = xy[a] - xy[b]
        total += int(math.sqrt(dx * dx + dy * dy) + 0.5) + int(cost[a])
    return total


# experiment
def plot(name, method, tour, xy, cost, value):
    fig, ax = plt.subplots(figsize=(8, 6))
    t = list(tour) + [tour[0]]
    ax.plot(xy[t, 0], xy[t, 1], "-", color="0.4", lw=0.9, zorder=1)
    sel = np.zeros(len(cost), dtype=bool)
    sel[tour] = True
    ax.scatter(xy[~sel, 0], xy[~sel, 1], c=cost[~sel], cmap="viridis", s=12,
               alpha=0.35, vmin=cost.min(), vmax=cost.max(), zorder=2)
    sc = ax.scatter(xy[sel, 0], xy[sel, 1], c=cost[sel], cmap="viridis", s=40,
                    edgecolors="black", linewidths=0.5, vmin=cost.min(), vmax=cost.max(), zorder=3)
    fig.colorbar(sc, ax=ax, label="node cost")
    ax.set_title(f"{name} - {method} (objective = {value})")
    ax.set_aspect("equal")
    fig.tight_layout()
    fname = OUT / f"{name}_{method.replace(' ', '_').replace('(', '').replace(')', '')}.png"
    fig.savefig(fname, dpi=120)
    plt.close(fig)
    return fname.name


def main():
    OUT.mkdir(exist_ok=True)
    rng = random.Random(SEED)
    lines = []
    summary = ["| Instance | Method | Min | Max | Avg |", "|---|---|---|---|---|"]
    for name in ("TSPA", "TSPB"):
        xy, cost = load_instance(name)
        D = distance_matrix(xy)  # optimisation methods see only D and cost
        n = len(cost)
        k = math.ceil(n / 2)
        for method, fn in METHODS.items():
            if fn is None:
                tours = [random_solution(D, cost, k, rng) for _ in range(N_RANDOM)]
            else:
                tours = [fn(D, cost, k, s) for s in range(n)]
            vals = [objective(t, D, cost) for t in tours]
            bi = int(np.argmin(vals))
            best = tours[bi]
            checked = check(best, xy, cost, k)
            assert checked == vals[bi]
            img = plot(name, method, best, xy, cost, vals[bi])
            summary.append(f"| {name} | {method} | {min(vals)} | {max(vals)} | {np.mean(vals):.1f} |")
            lines.append(f"### {name} - {method}\n\nObjective: {vals[bi]} "
                         f"(checker: {checked})\n\n![]({img})\n\n```\n{best}\n```\n")
            print(f"{name:5} {method:18} min={min(vals):7} max={max(vals):7} "
                  f"avg={np.mean(vals):9.1f}  checker OK")
    (OUT / "results.md").write_text("\n".join(summary) + "\n\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
