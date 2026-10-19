import math
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "lab1"))
import greedy  # noqa: E402  (lab1: data loading, objective, checker, plots, previous methods)

OUT = ROOT / "results"
greedy.OUT = OUT


def insertion_deltas(D, cost, seq, cand, cycle):
    """Rows = insertion positions, columns = candidates, value = objective increase.

    Row r means "insert at index r of seq" (path) or "insert into edge r" (cycle).
    """
    s = np.array(seq)
    if cycle:
        a, b = s, np.roll(s, -1)
        return D[a][:, cand] + D[b][:, cand] - D[a, b][:, None] + cost[cand][None, :]
    rows = [D[s[0], cand] + cost[cand]]
    if len(s) > 1:
        a, b = s[:-1], s[1:]
        rows.extend(D[a][:, cand] + D[b][:, cand] - D[a, b][:, None] + cost[cand][None, :])
    rows.append(D[s[-1], cand] + cost[cand])
    return np.vstack(rows)


def regret_construct(D, cost, k, start, cycle, w_regret, w_best):
    """Greedy construction choosing the node with the largest
    w_regret * (2nd best insertion - best insertion) - w_best * best insertion."""
    free = np.ones(len(cost), dtype=bool)
    free[start] = False
    seq = [start]
    if cycle:
        delta = D[start] + cost
        delta[~free] = np.iinfo(np.int64).max
        second = int(delta.argmin())
        seq.append(second)
        free[second] = False
    while len(seq) < k:
        cand = np.flatnonzero(free)
        delta = insertion_deltas(D, cost, seq, cand, cycle)
        two = np.partition(delta, 1, axis=0)
        best, regret = two[0], two[1] - two[0]
        score = w_regret * regret - w_best * best
        ties = np.flatnonzero(score == score.max())
        ci = ties[np.argmin(best[ties])]
        pos = int(delta[:, ci].argmin())
        seq.insert(pos + 1 if cycle else pos, int(cand[ci]))
        free[cand[ci]] = False
    return seq


def make(cycle, w_regret, w_best):
    return lambda D, cost, k, s: regret_construct(D, cost, k, s, cycle, w_regret, w_best)


PREVIOUS = dict(greedy.METHODS)
NEW = {
    "2-regret (NN any)": make(False, 1, 0),
    "2-regret (greedy cycle)": make(True, 1, 0),
    "Weighted (NN any)": make(False, 1, 1),
    "Weighted (greedy cycle)": make(True, 1, 1),
}
WEIGHTS = [(1, 1), (2, 1), (1, 2), (1, 4)]


def run(fn, D, cost, k, n, rng):
    if fn is None:
        return [greedy.random_solution(D, cost, k, rng) for _ in range(greedy.N_RANDOM)]
    return [fn(D, cost, k, s) for s in range(n)]


def main():
    OUT.mkdir(exist_ok=True)
    rng = random.Random(greedy.SEED)
    summary = ["| Instance | Method | Min | Max | Avg |", "|---|---|---|---|---|"]
    weights = ["| Instance | Variant | w_regret | w_best | Min | Max | Avg |",
               "|---|---|---|---|---|---|---|"]
    lines = []
    for name in ("TSPA", "TSPB"):
        xy, cost = greedy.load_instance(name)
        D = greedy.distance_matrix(xy)
        n = len(cost)
        k = math.ceil(n / 2)
        for method, fn in {**PREVIOUS, **NEW}.items():
            tours = run(fn, D, cost, k, n, rng)
            vals = [greedy.objective(t, D, cost) for t in tours]
            summary.append(f"| {name} | {method} | {min(vals)} | {max(vals)} | {np.mean(vals):.1f} |")
            print(f"{name:5} {method:24} min={min(vals):7} max={max(vals):7} avg={np.mean(vals):9.1f}")
            if method not in NEW:
                continue
            bi = int(np.argmin(vals))
            best = tours[bi]
            checked = greedy.check(best, xy, cost, k)
            assert checked == vals[bi]
            img = greedy.plot(name, method, best, xy, cost, vals[bi])
            lines.append(f"### {name} - {method}\n\nObjective: {vals[bi]} "
                         f"(checker: {checked})\n\n![]({img})\n\n```\n{best}\n```\n")
        for cycle in (False, True):
            variant = "greedy cycle" if cycle else "NN any"
            for wr, wb in WEIGHTS:
                vals = [greedy.objective(regret_construct(D, cost, k, s, cycle, wr, wb), D, cost)
                        for s in range(n)]
                weights.append(f"| {name} | {variant} | {wr} | {wb} | {min(vals)} | {max(vals)} "
                               f"| {np.mean(vals):.1f} |")
                print(f"{name:5} weights {variant:12} {wr}:{wb}  min={min(vals):7} avg={np.mean(vals):9.1f}")
    (OUT / "results.md").write_text("\n".join(summary) + "\n\n" + "\n".join(weights) + "\n\n"
                                    + "\n".join(lines))


if __name__ == "__main__":
    main()
