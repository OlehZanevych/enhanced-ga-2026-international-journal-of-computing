#!/usr/bin/env python3
"""analyze.py — statistical analysis of the experiments and generation of the
LaTeX tables / pgfplots figures used in the article.

Usage:
    python3 analyze.py [--results results] [--figures OUT_FIG] [--tables OUT_TAB]

Requires numpy and scipy. Reads results/{main,equal_time}/v*_*/seed*/
{summary,iterations}.csv written by run_experiments.sh and writes
  OUT_TAB/*.tex       LaTeX table bodies
  OUT_FIG/*.pgf       pgfplots figures (\\input{} in the article)
  OUT_TAB/numbers.json  every number quoted in the text
  OUT_TAB/report.md     human-readable summary
"""
import argparse
import csv
import glob
import json
import math
import os
from collections import defaultdict

import numpy as np
from scipy import stats

MODES = ["no_ls", "ls", "ls_w", "amls"]
LABEL = {"no_ls": "GA", "ls": "GA+LS", "ls_w": "GA+LS-W", "amls": "GA+AMLS"}
COLOR = {"no_ls": "blue", "ls": "red", "ls_w": "orange!90!black", "amls": "green!50!black"}
MARK = {"no_ls": "o", "ls": "square", "ls_w": "triangle", "amls": "diamond"}
VARIANTS = [0, 1, 2, 3, 4]
VSIZE = {0: 30, 1: 60, 2: 120, 3: 250, 4: 500}


# ── loading ───────────────────────────────────────────────────────────────────
def load_summaries(exp_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(exp_dir, "v*_*", "seed*", "summary.csv"))):
        with open(f) as fh:
            r = next(csv.DictReader(fh))
        r["_dir"] = os.path.dirname(f)
        rows.append(r)
    out = defaultdict(dict)       # (variant, mode) -> seed -> row
    for r in rows:
        out[(int(r["variant"]), r["mode"])][int(r["seed"])] = r
    return out


def load_trace(run_dir):
    it, best, hard, t = [], [], [], []
    with open(os.path.join(run_dir, "iterations.csv")) as fh:
        rd = csv.reader(fh)
        next(rd)
        for row in rd:
            it.append(int(row[0])); best.append(float(row[1])); hard.append(float(row[2])); t.append(float(row[5]))
    return np.array(it), np.array(best), np.array(hard), np.array(t)


def col(runs, key, seeds=None):
    seeds = sorted(runs) if seeds is None else seeds
    return np.array([float(runs[s][key]) for s in seeds])


# ── statistics ────────────────────────────────────────────────────────────────
def wilcoxon(x, y):
    """Two-sided Wilcoxon signed-rank test, zero differences discarded."""
    d = np.asarray(x) - np.asarray(y)
    nz = int(np.count_nonzero(d))
    if nz == 0:
        return {"n": 0, "T": 0.0, "p": 1.0}
    res = stats.wilcoxon(x, y, zero_method="wilcox", alternative="two-sided")
    return {"n": nz, "T": float(res.statistic), "p": float(res.pvalue)}


def a12(x, y):
    """Vargha–Delaney A: probability that a run of x is better (lower) than a run of y."""
    x = np.asarray(x); y = np.asarray(y)
    gt = (x[:, None] < y[None, :]).sum()
    eq = (x[:, None] == y[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(x) * len(y)))


def mcnemar_exact(a, b):
    """Exact McNemar test for paired binary outcomes a, b (True = success)."""
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    n01 = int(np.sum(a & ~b)); n10 = int(np.sum(~a & b))
    n = n01 + n10
    if n == 0:
        return {"a_only": n01, "b_only": n10, "p": 1.0}
    p = min(1.0, 2 * stats.binom.cdf(min(n01, n10), n, 0.5))
    return {"a_only": n01, "b_only": n10, "p": float(p)}


def holm(pvals):
    """Holm–Bonferroni adjusted p-values (same order as input)."""
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * pvals[i])
        adj[i] = min(1.0, running)
    return adj.tolist()


def stars(p):
    return "$^{***}$" if p < 0.001 else "$^{**}$" if p < 0.01 else "$^{*}$" if p < 0.05 else "ns"


def fmt_p(p):
    if p < 1e-3:
        e = int(math.floor(math.log10(p)))
        return f"${p / 10**e:.1f}\\cdot10^{{{e}}}$"
    return f"${p:.4f}$"


def fmed(x):
    """Median of integer counts: integer if whole, else one decimal (no banker's rounding)."""
    return f"{x:.0f}" if float(x).is_integer() else f"{x:.1f}"


def fnum(x, nd=1):
    return f"{x:.{nd}f}"


# ── tables ────────────────────────────────────────────────────────────────────
def table_main(main, N):
    lines = []
    for v in VARIANTS:
        modes = [m for m in MODES if (v, m) in main]
        means = {m: col(main[(v, m)], "final_best_score").mean() for m in modes}
        best_mean = min(means.values())
        for k, m in enumerate(modes):
            r = main[(v, m)]
            f = col(r, "final_best_score"); h = col(r, "hard_penalty"); w = col(r, "soft_penalty")
            feas = int(col(r, "feasible").sum()); conv = int(col(r, "converged").sum())
            ms = f"{f.mean():.1f} \\pm {f.std(ddof=1):.1f}" if len(f) > 1 else f"{f.mean():.1f}"
            ms = f"$\\mathbf{{{ms}}}$" if abs(means[m] - best_mean) < 1e-9 and len(modes) > 1 else f"${ms}$"
            pre = (f"\\multirow{{{len(modes)}}}{{*}}{{{VSIZE[v]}}}" if k == 0 else "")
            lines.append(f"  {pre} & {LABEL[m]} & {feas}/{len(f)} & {conv}/{len(f)} & {ms} & "
                         f"{np.median(f):.1f} & {f.min():.1f} & {f.max():.1f} & {h.mean():.1f} & {w.mean():.1f} \\\\")
        lines.append("\\midrule" if v != VARIANTS[-1] else "")
    return "\n".join(l for l in lines if l)


COMPARISONS = [("amls", "no_ls"), ("amls", "ls_w"), ("ls", "no_ls"), ("ls_w", "ls")]


def stats_tests(main):
    """Pairwise comparisons on final score (paired by seed) + McNemar on feasibility/optimality."""
    out = []
    for (x, y) in COMPARISONS:
        for v in VARIANTS:
            if (v, x) not in main or (v, y) not in main:
                continue
            seeds = sorted(set(main[(v, x)]) & set(main[(v, y)]))
            fx = col(main[(v, x)], "final_best_score", seeds); fy = col(main[(v, y)], "final_best_score", seeds)
            w = wilcoxon(fx, fy)
            out.append({
                "x": x, "y": y, "variant": v, "V": VSIZE[v], "n_seeds": len(seeds),
                "wilcoxon_n": w["n"], "T": w["T"], "p": w["p"],
                "A12": a12(fx, fy),
                "wins": int(np.sum(fx < fy)), "ties": int(np.sum(fx == fy)), "losses": int(np.sum(fx > fy)),
                "mean_x": float(fx.mean()), "mean_y": float(fy.mean()),
                "feas": mcnemar_exact(col(main[(v, x)], "feasible", seeds) > 0, col(main[(v, y)], "feasible", seeds) > 0),
                "opt": mcnemar_exact(col(main[(v, x)], "converged", seeds) > 0, col(main[(v, y)], "converged", seeds) > 0),
            })
    adj = holm([o["p"] for o in out])
    for o, a in zip(out, adj):
        o["p_holm"] = a
    return out


def table_stats(tests):
    lines = []
    cur = None
    for o in tests:
        key = (o["x"], o["y"])
        if key != cur:
            if cur is not None:
                lines.append("\\midrule")
            cur = key
            lines.append(f"\\multicolumn{{8}}{{l}}{{\\textit{{{LABEL[o['x']]} vs.\\ {LABEL[o['y']]}}}}} \\\\")
        p_disp = fmt_p(o["p_holm"]) if o["wilcoxon_n"] > 0 else "---"
        sig = stars(o["p_holm"]) if o["wilcoxon_n"] > 0 else "---"
        imp = (1 - o["mean_x"] / o["mean_y"]) * 100 if o["mean_y"] > 0 else float("nan")
        imp_s = f"{imp:.1f}\\%" if not math.isnan(imp) else "---"
        lines.append(f"  {o['V']} & {imp_s} & {o['wilcoxon_n']} & {p_disp} & {sig} & {o['A12']:.2f} & "
                     f"{o['wins']}/{o['ties']}/{o['losses']} & {fmt_p(o['feas']['p']) if o['feas']['a_only'] + o['feas']['b_only'] else '---'} \\\\")
    return "\n".join(lines)


def table_speed(main):
    lines = []
    ga_ms = {}
    for v in VARIANTS:
        if (v, "no_ls") in main:
            ga_ms[v] = 1000 * col(main[(v, "no_ls")], "avg_iter_time_sec").mean()
    for v in VARIANTS:
        modes = [m for m in MODES if (v, m) in main]
        for k, m in enumerate(modes):
            r = main[(v, m)]
            seeds = sorted(r)
            ffi = np.array([float(r[s]["first_feasible_iter"]) for s in seeds if int(r[s]["feasible"]) == 1])
            ffs = np.array([float(r[s]["first_feasible_sec"]) for s in seeds if int(r[s]["feasible"]) == 1])
            conv_it = np.array([float(r[s]["iterations_run"]) for s in seeds if int(r[s]["converged"]) == 1])
            ms = 1000 * col(r, "avg_iter_time_sec").mean()
            share = (col(r, "ls_time_sec").sum() / col(r, "total_time_sec").sum() * 100) if m != "no_ls" else 0.0
            pre = (f"\\multirow{{{len(modes)}}}{{*}}{{{VSIZE[v]}}}" if k == 0 else "")
            ffi_s = fmed(np.median(ffi)) if len(ffi) else "---"
            ffs_s = f"{1000 * np.median(ffs):.1f}" if len(ffs) else "---"
            conv_s = fmed(np.median(conv_it)) if len(conv_it) else "---"
            ovh = f"{ms / ga_ms[v]:.2f}$\\times$" if v in ga_ms else "---"
            share_s = f"{share:.0f}\\%" if m != "no_ls" else "---"
            lines.append(f"  {pre} & {LABEL[m]} & {ffi_s} & {ffs_s} & {conv_s} & {ms:.3f} & {ovh} & {share_s} \\\\")
        lines.append("\\midrule" if v != VARIANTS[-1] else "")
    return "\n".join(l for l in lines if l)


def table_ls_effect(main):
    lines = []
    for v in VARIANTS:
        cells = []
        for m in ["ls", "ls_w", "amls"]:
            if (v, m) not in main:
                cells.append("--- & ---"); continue
            r = main[(v, m)]
            n = col(r, "ls_sampled").sum()
            b = col(r, "ls_conflicted_before").sum() / n if n else float("nan")
            a = col(r, "ls_conflicted_after").sum() / n if n else float("nan")
            cells.append(f"{b:.2f} & {a:.2f}")
        lines.append(f"  {VSIZE[v]} & " + " & ".join(cells) + " \\\\")
    return "\n".join(lines)


def equal_time_tests(eq):
    """Wilcoxon tests of every LS variant against GA in E2, Holm-adjusted."""
    tests = []
    for v in [3, 4]:
        if (v, "no_ls") not in eq:
            continue
        seeds = sorted(eq[(v, "no_ls")])
        ga = col(eq[(v, "no_ls")], "final_best_score", seeds)
        for m in MODES[1:]:
            if (v, m) not in eq:
                continue
            f = col(eq[(v, m)], "final_best_score", seeds)
            w = wilcoxon(f, ga)
            tests.append({"variant": v, "mode": m, "n": w["n"], "p": w["p"], "A12": a12(f, ga),
                          "feas": mcnemar_exact(col(eq[(v, m)], "feasible", seeds) > 0,
                                                col(eq[(v, "no_ls")], "feasible", seeds) > 0)})
    for t, a in zip(tests, holm([t["p"] for t in tests])):
        t["p_holm"] = a
    return tests


def table_equal_time(eq, budgets, tests):
    lines = []
    tmap = {(t["variant"], t["mode"]): t for t in tests}
    for v in [3, 4]:
        modes = [m for m in MODES if (v, m) in eq]
        if not modes:
            continue
        means = {m: col(eq[(v, m)], "final_best_score").mean() for m in modes}
        best_mean = min(means.values())
        for k, m in enumerate(modes):
            r = eq[(v, m)]
            f = col(r, "final_best_score")
            ms = f"{f.mean():.1f} \\pm {f.std(ddof=1):.1f}" if len(f) > 1 else f"{f.mean():.1f}"
            ms = f"$\\mathbf{{{ms}}}$" if abs(means[m] - best_mean) < 1e-9 else f"${ms}$"
            its = col(r, "iterations_run")
            t = tmap.get((v, m))
            p_s = f"{fmt_p(t['p_holm'])} {stars(t['p_holm'])}" if t and t["n"] else "---"
            a_s = f"{t['A12']:.2f}" if t else "---"
            pre = (f"\\multirow{{{len(modes)}}}{{*}}{{{VSIZE[v]} ({budgets.get(v, float('nan')):.1f}\\,s)}}" if k == 0 else "")
            lines.append(f"  {pre} & {LABEL[m]} & {fmed(np.median(its))} & {int(col(r, 'feasible').sum())}/{len(f)} & "
                         f"{int(col(r, 'converged').sum())}/{len(f)} & {ms} & {np.median(f):.1f} & {p_s} & {a_s} \\\\")
        if v == 3:
            lines.append("\\midrule")
    return "\n".join(lines)


# ── figures (pgfplots) ────────────────────────────────────────────────────────
def mean_curve(runs, key_idx, n_points, max_iter):
    """Mean of best-so-far curves over a common iteration grid (runs that stopped
    early at f = 0 keep their final value)."""
    grid = np.unique(np.round(np.geomspace(1, max_iter, n_points)).astype(int))
    mats = []
    for s, r in runs.items():
        it, best, hard, t = load_trace(r["_dir"])
        y = best if key_idx == "best" else hard
        if len(it) == 0:
            mats.append(np.full(len(grid), float(r["final_best_score"]))); continue
        idx = np.searchsorted(it, grid, side="right") - 1
        vals = np.where(idx >= 0, y[np.clip(idx, 0, len(y) - 1)], y[0])
        vals = np.where(grid > it[-1], y[-1], vals)
        mats.append(vals)
    M = np.array(mats)
    return grid, M


def feasible_fraction_curve(runs, n_points, max_iter):
    """Fraction of runs whose best-so-far schedule is feasible (H = 0) at each
    generation of a common grid. Because f is a weighted penalty, the best
    schedule can become infeasible again (a schedule with a minor clash may have
    a lower f than a feasible one with windows), so the curve need not be monotone;
    its last value equals the number of runs that end feasible."""
    grid = np.unique(np.round(np.geomspace(1, max_iter, n_points)).astype(int))
    grid, H = mean_curve(runs, "hard", n_points, max_iter)
    return grid, (H == 0).mean(axis=0)


def pgf_axis(opts, plots):
    return "\\begin{tikzpicture}\n\\begin{axis}[\n  " + ",\n  ".join(opts) + "\n]\n" + "\n".join(plots) + "\n\\end{axis}\n\\end{tikzpicture}\n"


def coords(xs, ys, fx="{:.4g}", fy="{:.4g}"):
    return " ".join(f"({fx.format(x)},{fy.format(y)})" for x, y in zip(xs, ys))


def fig_convergence(main, v, path, width, legend=True, max_iter=5000):
    plots = []
    for m in MODES:
        if (v, m) not in main:
            continue
        grid, M = mean_curve(main[(v, m)], "best", 60, max_iter)
        y = M.mean(axis=0) + 1.0
        plots.append(f"\\addplot[{COLOR[m]}, thick, mark={MARK[m]}, mark repeat=6, mark size=1.3pt] coordinates {{{coords(grid, y)}}};")
        if legend:
            plots.append(f"\\addlegendentry{{{LABEL[m]}}}")
    opts = [f"width={width}", "height=4.6cm", "xmode=log", "ymode=log",
            "xlabel={Generation}", "ylabel={Mean best $f$ + 1}", f"title={{$V = {VSIZE[v]}$}}",
            "xmin=1", f"xmax={max_iter}", "grid=major", "grid style={gray!30}",
            "legend style={font=\\scriptsize, at={(0.03,0.03)}, anchor=south west}",
            "tick label style={font=\\scriptsize}", "label style={font=\\small}", "title style={font=\\small}"]
    with open(path, "w") as fh:
        fh.write(pgf_axis(opts, plots))


def fig_feasible(main, v, path, width, legend=True, max_iter=5000):
    plots = []
    for m in MODES:
        if (v, m) not in main:
            continue
        grid, fr = feasible_fraction_curve(main[(v, m)], 80, max_iter)
        plots.append(f"\\addplot[{COLOR[m]}, thick, const plot] coordinates {{{coords(grid, 100 * fr, fy='{:.1f}')}}};")
        if legend:
            plots.append(f"\\addlegendentry{{{LABEL[m]}}}")
    opts = [f"width={width}", "height=4.6cm", "xmode=log", "xlabel={Generation}",
            "ylabel={Runs with feasible $\\sigma^*$ (\\%)}", f"title={{$V = {VSIZE[v]}$}}", "ymin=0", "ymax=102",
            "xmin=1", f"xmax={max_iter}", "grid=major", "grid style={gray!30}",
            "legend style={font=\\scriptsize, at={(0.03,0.97)}, anchor=north west}",
            "tick label style={font=\\scriptsize}", "label style={font=\\small}", "title style={font=\\small}"]
    with open(path, "w") as fh:
        fh.write(pgf_axis(opts, plots))


def fig_boxplot(main, v, path, width):
    plots = []
    labels = []
    pos = 0
    for m in MODES:
        if (v, m) not in main:
            continue
        pos += 1
        labels.append(LABEL[m])
        f = col(main[(v, m)], "final_best_score")
        q1, med, q3 = np.percentile(f, [25, 50, 75])
        iqr = q3 - q1
        lw = f[f >= q1 - 1.5 * iqr].min(); uw = f[f <= q3 + 1.5 * iqr].max()
        outl = f[(f < lw) | (f > uw)]
        oc = " ".join(f"({pos},{o:g})" for o in outl)
        plots.append(f"\\addplot[boxplot prepared={{lower whisker={lw:g}, lower quartile={q1:g}, median={med:g}, "
                     f"upper quartile={q3:g}, upper whisker={uw:g}}}, fill={COLOR[m]}!35, draw=black, "
                     f"mark=*, mark size=0.8pt] coordinates {{{oc}}};")
    opts = [f"width={width}", "height=4.8cm", "boxplot/draw direction=y", "ylabel={Final $f$}",
            f"title={{$V = {VSIZE[v]}$}}", "xtick={" + ",".join(str(i + 1) for i in range(pos)) + "}",
            "xticklabels={" + ",".join(labels) + "}", "x tick label style={font=\\scriptsize, rotate=30, anchor=east}",
            "grid=major", "grid style={gray!30}", "ymin=-1", "tick label style={font=\\scriptsize}",
            "label style={font=\\small}", "title style={font=\\small}"]
    with open(path, "w") as fh:
        fh.write(pgf_axis(opts, plots))


def fig_walltime(eq, v, path, width, budget, legend=True):
    plots = []
    grid = np.linspace(0, budget, 120)[1:]
    for m in MODES:
        if (v, m) not in eq:
            continue
        mats = []
        for s, r in eq[(v, m)].items():
            it, best, hard, t = load_trace(r["_dir"])
            if len(t) == 0:
                mats.append(np.full(len(grid), float(r["final_best_score"]))); continue
            idx = np.searchsorted(t, grid, side="right") - 1
            first = float(best[0])
            vals = np.where(idx >= 0, best[np.clip(idx, 0, len(best) - 1)], np.nan if first is None else first)
            vals = np.where(grid > t[-1], best[-1], vals)
            mats.append(vals)
        y = np.nanmean(np.array(mats), axis=0) + 1.0
        plots.append(f"\\addplot[{COLOR[m]}, thick, mark={MARK[m]}, mark repeat=12, mark size=1.3pt] coordinates {{{coords(grid, y)}}};")
        if legend:
            plots.append(f"\\addlegendentry{{{LABEL[m]}}}")
    opts = [f"width={width}", "height=4.8cm", "ymode=log", "xlabel={Wall-clock time (s)}",
            "ylabel={Mean best $f$ + 1}", f"title={{$V = {VSIZE[v]}$, equal budget {budget:.1f}\\,s}}",
            "xmin=0", f"xmax={budget:.3f}", "grid=major", "grid style={gray!30}",
            "legend style={font=\\scriptsize, at={(0.97,0.97)}, anchor=north east}",
            "tick label style={font=\\scriptsize}", "label style={font=\\small}", "title style={font=\\small}"]
    with open(path, "w") as fh:
        fh.write(pgf_axis(opts, plots))


def fig_iter_time(main, path, width):
    plots = []
    vs = [v for v in VARIANTS if (v, "no_ls") in main]
    for m in MODES:
        pts = [(VSIZE[v], 1000 * col(main[(v, m)], "avg_iter_time_sec").mean()) for v in vs if (v, m) in main]
        plots.append(f"\\addplot[{COLOR[m]}, thick, mark={MARK[m]}, mark size=1.6pt] coordinates {{"
                     + " ".join(f"({x},{y:.4f})" for x, y in pts) + "};")
        plots.append(f"\\addlegendentry{{{LABEL[m]}}}")
    opts = [f"width={width}", "height=4.8cm", "xmode=log", "ymode=log", "xlabel={Problem size $V$}",
            "ylabel={Time per iteration (ms)}", "xtick={30,60,120,250,500}", "xticklabels={30,60,120,250,500}",
            "grid=major", "grid style={gray!30}", "legend style={font=\\scriptsize, at={(0.03,0.97)}, anchor=north west}",
            "tick label style={font=\\scriptsize}", "label style={font=\\small}"]
    with open(path, "w") as fh:
        fh.write(pgf_axis(opts, plots))


# ── main ──────────────────────────────────────────────────────────────────────
def describe(runs):
    f = col(runs, "final_best_score")
    return {"n": len(f), "mean": float(f.mean()), "std": float(f.std(ddof=1)) if len(f) > 1 else 0.0,
            "median": float(np.median(f)), "min": float(f.min()), "max": float(f.max()),
            "feasible": int(col(runs, "feasible").sum()), "converged": int(col(runs, "converged").sum()),
            "hard_mean": float(col(runs, "hard_penalty").mean()), "soft_mean": float(col(runs, "soft_penalty").mean()),
            "iter_ms": float(1000 * col(runs, "avg_iter_time_sec").mean()),
            "time_mean": float(col(runs, "total_time_sec").mean()),
            "iters_mean": float(col(runs, "iterations_run").mean())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--figures", default="analysis_out/figures")
    ap.add_argument("--tables", default="analysis_out/tables")
    a = ap.parse_args()
    os.makedirs(a.figures, exist_ok=True); os.makedirs(a.tables, exist_ok=True)

    main_r = load_summaries(os.path.join(a.results, "main"))
    eq_r = load_summaries(os.path.join(a.results, "equal_time"))
    budgets = {}
    bf = os.path.join(a.results, "equal_time_budgets.txt")
    if os.path.exists(bf):
        for line in open(bf):
            parts = dict(p.split("=") for p in line.split())
            budgets[int(parts["variant"])] = float(parts["time_budget_sec"])

    numbers = {"main": {}, "equal_time": {}, "budgets": budgets}
    for (v, m), runs in main_r.items():
        numbers["main"][f"v{v}_{m}"] = describe(runs)
        r = runs
        n = col(r, "ls_sampled").sum()
        if n:
            numbers["main"][f"v{v}_{m}"]["conflicted_before"] = float(col(r, "ls_conflicted_before").sum() / n)
            numbers["main"][f"v{v}_{m}"]["conflicted_after"] = float(col(r, "ls_conflicted_after").sum() / n)
            numbers["main"][f"v{v}_{m}"]["ls_time_share"] = float(col(r, "ls_time_sec").sum() / col(r, "total_time_sec").sum())
        feas = [int(r[s]["first_feasible_iter"]) for s in sorted(r) if int(r[s]["feasible"])]
        numbers["main"][f"v{v}_{m}"]["first_feasible_iter_median"] = float(np.median(feas)) if feas else None
        numbers["main"][f"v{v}_{m}"]["ever_feasible"] = int(sum(int(r[s]["first_feasible_iter"]) >= 0 for s in r))
        ffs = [1000 * float(r[s]["first_feasible_sec"]) for s in sorted(r) if int(r[s]["feasible"])]
        numbers["main"][f"v{v}_{m}"]["first_feasible_ms_median"] = float(np.median(ffs)) if ffs else None
        conv = [int(r[s]["iterations_run"]) for s in sorted(r) if int(r[s]["converged"])]
        numbers["main"][f"v{v}_{m}"]["conv_iter_median"] = float(np.median(conv)) if conv else None
        numbers["main"][f"v{v}_{m}"]["conv_iter_mean"] = float(np.mean(conv)) if conv else None
    for (v, m), runs in eq_r.items():
        numbers["equal_time"][f"v{v}_{m}"] = describe(runs)
        numbers["equal_time"][f"v{v}_{m}"]["iters_median"] = float(np.median(col(runs, "iterations_run")))

    tests = stats_tests(main_r) if main_r else []
    numbers["tests"] = tests

    # tables
    if main_r:
        open(os.path.join(a.tables, "results_main.tex"), "w").write(table_main(main_r, 50) + "\n")
        open(os.path.join(a.tables, "stats.tex"), "w").write(table_stats(tests) + "\n")
        open(os.path.join(a.tables, "speed.tex"), "w").write(table_speed(main_r) + "\n")
        open(os.path.join(a.tables, "ls_effect.tex"), "w").write(table_ls_effect(main_r) + "\n")
    if eq_r:
        eq_tests = equal_time_tests(eq_r)
        numbers["equal_time_tests"] = eq_tests
        open(os.path.join(a.tables, "equal_time.tex"), "w").write(table_equal_time(eq_r, budgets, eq_tests) + "\n")

    # figures
    if main_r:
        for v, w in [(1, "\\linewidth"), (2, "\\linewidth"), (3, "\\linewidth"), (4, "\\linewidth")]:
            fig_convergence(main_r, v, os.path.join(a.figures, f"convergence_v{VSIZE[v]}.pgf"), w, legend=(v == 1))
        for v in [3, 4]:
            fig_feasible(main_r, v, os.path.join(a.figures, f"feasible_v{VSIZE[v]}.pgf"), "\\linewidth", legend=(v == 3))
        for v in [2, 3, 4]:
            fig_boxplot(main_r, v, os.path.join(a.figures, f"boxplot_v{VSIZE[v]}.pgf"), "\\linewidth")
        fig_iter_time(main_r, os.path.join(a.figures, "iteration_time.pgf"), "\\columnwidth")
    if eq_r:
        for v in [3, 4]:
            if v in budgets:
                fig_walltime(eq_r, v, os.path.join(a.figures, f"walltime_v{VSIZE[v]}.pgf"), "\\linewidth",
                             budgets[v], legend=(v == 3))

    with open(os.path.join(a.tables, "numbers.json"), "w") as fh:
        json.dump(numbers, fh, indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item())

    # report
    with open(os.path.join(a.tables, "report.md"), "w") as fh:
        for name, d in [("E1 main", numbers["main"]), ("E2 equal_time", numbers["equal_time"])]:
            fh.write(f"## {name}\n\n| run | n | feasible | f=0 | mean | std | median | min | max | H | W | ms/iter |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
            for k in sorted(d, key=lambda s: (int(s[1]), MODES.index(s.split('_', 1)[1]))):
                x = d[k]
                fh.write(f"| {k} | {x['n']} | {x['feasible']} | {x['converged']} | {x['mean']:.1f} | {x['std']:.1f} | "
                         f"{x['median']:.1f} | {x['min']:.1f} | {x['max']:.1f} | {x['hard_mean']:.1f} | {x['soft_mean']:.1f} | {x['iter_ms']:.3f} |\n")
            fh.write("\n")
        fh.write("## Tests (Holm-adjusted over all comparisons)\n\n")
        for o in tests:
            fh.write(f"- V={o['V']} {LABEL[o['x']]} vs {LABEL[o['y']]}: mean {o['mean_x']:.1f} vs {o['mean_y']:.1f}, "
                     f"Wilcoxon n={o['wilcoxon_n']} p={o['p']:.3g} p_holm={o['p_holm']:.3g}, A12={o['A12']:.2f}, "
                     f"W/T/L={o['wins']}/{o['ties']}/{o['losses']}, McNemar feas p={o['feas']['p']:.3g} "
                     f"({o['feas']['a_only']}:{o['feas']['b_only']}), opt p={o['opt']['p']:.3g} ({o['opt']['a_only']}:{o['opt']['b_only']})\n")
    print("Analysis written to", a.tables, "and", a.figures)


if __name__ == "__main__":
    main()
