#!/usr/bin/env python3
"""Reproducible simulations for the HDLSS PCA manuscript.

The code works entirely in dual/sample space.  It generates Gaussian factors,
which satisfy the independent-factor sub-Gaussian assumption, and records both
population and finite-score-oracle errors.  All random-number streams are
created from the master seed through NumPy SeedSequence.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MASTER_SEED = 20261109


def center_matrix(n: int) -> np.ndarray:
    return np.eye(n) - np.ones((n, n)) / n


def summaries(values: np.ndarray, experiment: str, d: int, n: int, m: int,
              alpha: str, estimator: str, metric: str, reps: int) -> dict:
    values = np.asarray(values, dtype=float)
    return {
        "experiment": experiment, "distribution": "Gaussian", "d": d,
        "n": n, "m": m, "alpha": alpha, "replicates": reps,
        "estimator": estimator, "metric": metric,
        "mean": float(np.mean(values)),
        "mc_se": float(np.std(values, ddof=1) / math.sqrt(values.size)),
        "median": float(np.median(values)),
        "p90": float(np.quantile(values, 0.90)),
    }


def simulate(d: int, n: int, lambdas: np.ndarray, reps: int, rng: np.random.Generator,
             batch: int = 40, need_cdm: bool = False, need_direction: bool = False):
    m = len(lambdas); r = n - 1
    P = center_matrix(n)
    out = {k: [] for k in ["pca", "bnrm", "nrm", "score", "cdm", "cdm_score", "sin_angle"]}
    for start in range(0, reps, batch):
        b = min(batch, reps - start)
        zs = rng.standard_normal((b, m, n))
        z0 = rng.standard_normal((b, d - m, n))
        sig_gram = np.einsum("bmn,m,bmk->bnk", zs, lambdas, zs, optimize=True)
        noise_gram = np.einsum("bdn,bdk->bnk", z0, z0, optimize=True)
        gram = np.einsum("ij,bjk,kl->bil", P, sig_gram + noise_gram, P, optimize=True) / r
        score_gram = np.einsum("ij,bjk,kl->bil", P, sig_gram, P, optimize=True) / r
        evals, evecs = np.linalg.eigh(gram)
        sevals = np.linalg.eigvalsh(score_gram)
        a = evals[:, ::-1][:, :r]
        V = evecs[:, :, ::-1]
        s = sevals[:, ::-1][:, :m]
        nu = np.mean(a[:, m:r], axis=1)
        bnrm = a[:, :m] - nu[:, None]
        nrm = np.empty((b, m))
        for i in range(m):
            nrm[:, i] = a[:, i] - np.mean(a[:, i+1:r], axis=1)
        out["pca"].append(a[:, :m]); out["bnrm"].append(bnrm)
        out["nrm"].append(nrm); out["score"].append(s)
        if need_direction:
            v = V[:, :, 0]
            sigcoef = np.sqrt(lambdas[0]) * np.einsum("bn,bn->b", zs[:, 0, :], v)
            z0v = np.einsum("bdn,bn->bd", z0, v, optimize=True)
            noise2 = np.einsum("bd,bd->b", z0v, z0v)
            out["sin_angle"].append(np.sqrt(noise2 / (sigcoef**2 + noise2)))
        if need_cdm:
            n1 = n // 2; n2 = n - n1; r1 = n1 - 1; r2 = n2 - 1
            P1, P2 = center_matrix(n1), center_matrix(n2)
            zsa, zsb = zs[:, :, :n1], zs[:, :, n1:]
            z0a, z0b = z0[:, :, :n1], z0[:, :, n1:]
            cross_sig = np.einsum("bmi,m,bmj->bij", zsa, lambdas, zsb, optimize=True)
            cross_noise = np.einsum("bdi,bdj->bij", z0a, z0b, optimize=True)
            scale = math.sqrt(r1 * r2)
            C = np.einsum("ij,bjk,kl->bil", P1, cross_sig + cross_noise, P2, optimize=True) / scale
            Cs = np.einsum("ij,bjk,kl->bil", P1, cross_sig, P2, optimize=True) / scale
            out["cdm"].append(np.linalg.svd(C, compute_uv=False)[:, :m])
            out["cdm_score"].append(np.linalg.svd(Cs, compute_uv=False)[:, :m])
    return {k: (np.concatenate(v, axis=0) if v else None) for k, v in out.items()}


def plot_results(df: pd.DataFrame, figdir: Path):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    colors = {"PCA": "#777777", "bNRM": "#006BA4", "CDM": "#FF800E", "sequential NRM": "#ABABAB", "theoretical target": "#C85200"}

    # Fixed-score population versus oracle errors.
    sub = df[(df.experiment == "fixed_rank1") & (df.metric.isin(["population relative error", "score-oracle relative error"]))]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), sharex=True)
    for j, metric in enumerate(["population relative error", "score-oracle relative error"]):
        for est in ["PCA", "bNRM", "CDM"]:
            q = sub[(sub.metric == metric) & (sub.estimator == est)]
            if len(q): ax[j].plot(q.d, q["median"], marker="o", label=est, color=colors[est])
        ax[j].set_xscale("log"); ax[j].set_xticks([250, 500, 1000, 2000]); ax[j].set_xticklabels(["250", "500", "1k", "2k"]); ax[j].minorticks_off()
        ax[j].set_yscale("log"); ax[j].grid(alpha=.25)
        ax[j].set_xlabel("dimension $d$"); ax[j].set_title(metric.capitalize(), fontsize=10)
    ax[0].set_ylabel("median relative error")
    ax[1].legend(frameon=False)
    fig.tight_layout(); fig.savefig(figdir / "fixed_score_errors.pdf", bbox_inches="tight"); fig.savefig(figdir / "fixed_score_errors.png", dpi=220, bbox_inches="tight"); plt.close(fig)

    # Multiple-spike sequential subtraction.
    sub = df[(df.experiment == "two_equal_spikes") & (df.metric == "estimate / score oracle")]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for est in ["sequential NRM", "theoretical target", "bNRM"]:
        q = sub[sub.estimator == est]
        ax.plot(q.d, q["median"], marker="o", label=est, color=colors[est])
    ax.axhline(1, lw=.8, color="black"); ax.set_xscale("log"); ax.set_xticks([250, 500, 1000, 2000]); ax.set_xticklabels(["250", "500", "1k", "2k"]); ax.minorticks_off(); ax.set_xlabel("dimension $d$")
    ax.set_ylabel("median ratio, component 1"); ax.grid(alpha=.25); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(figdir / "multispike_subtraction.pdf", bbox_inches="tight"); fig.savefig(figdir / "multispike_subtraction.png", dpi=220, bbox_inches="tight"); plt.close(fig)

    # Direction phase separation.
    sub = df[(df.experiment == "direction_phase") & (df.metric == "sine angle")]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for alpha, color in [(0.7, "#006BA4"), (1.2, "#FF800E")]:
        q = sub[np.isclose(sub.alpha.astype(float), alpha)]
        ax.plot(q.d, q["median"], marker="o", label=rf"$\alpha={alpha}$", color=color)
    ax.set_xscale("log"); ax.set_xticks([250, 500, 1000, 2000]); ax.set_xticklabels(["250", "500", "1k", "2k"]); ax.minorticks_off(); ax.set_xlabel("dimension $d$"); ax.set_ylabel("median sine angle")
    ax.grid(alpha=.25); ax.legend(frameon=False, loc="center right")
    fig.tight_layout(); fig.savefig(figdir / "direction_phase.pdf", bbox_inches="tight"); fig.savefig(figdir / "direction_phase.png", dpi=220, bbox_inches="tight"); plt.close(fig)

    # Growing-sample population transition.
    sub = df[df.experiment == "jointhdlss"]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for metric, ls in [("population relative error", "-"), ("score-oracle relative error", "--")]:
        q = sub[(sub.estimator == "bNRM") & (sub.metric == metric)]
        ax.plot(q.d, q["median"], marker="o", ls=ls, label=metric)
    ax.set_xscale("log"); ax.set_xticks([400, 900, 1600, 2500]); ax.set_xticklabels(["400", "900", "1.6k", "2.5k"]); ax.minorticks_off(); ax.set_yscale("log"); ax.set_xlabel("dimension $d$")
    ax.set_ylabel("median relative error"); ax.grid(alpha=.25); ax.legend(frameon=False, loc="center right")
    fig.tight_layout(); fig.savefig(figdir / "joint_population_transition.pdf", bbox_inches="tight"); fig.savefig(figdir / "joint_population_transition.png", dpi=220, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results"); ap.add_argument("--reps", type=int, default=1000)
    args = ap.parse_args(); outdir = Path(args.out); figdir = outdir / "figures"; outdir.mkdir(parents=True, exist_ok=True); figdir.mkdir(exist_ok=True)
    reps = args.reps; rows = []
    ss = np.random.SeedSequence(MASTER_SEED); children = iter(ss.spawn(100))

    # Experiment 1: fixed n, rank one; distinguish population and score-oracle errors.
    for d in [250, 500, 1000, 2000]:
        n = 20; lam = d ** 0.7
        res = simulate(d, n, np.array([lam]), reps, np.random.default_rng(next(children)), need_cdm=True)
        for est, arr, oracle in [("PCA", res["pca"], res["score"]), ("bNRM", res["bnrm"], res["score"]), ("CDM", res["cdm"], res["cdm_score"])]:
            poperr = np.abs(arr[:, 0] - lam) / lam
            orerr = np.abs(arr[:, 0] - oracle[:, 0]) / np.maximum(oracle[:, 0], 1e-12)
            rows += [summaries(poperr, "fixed_rank1", d, n, 1, "0.7", est, "population relative error", reps),
                     summaries(orerr, "fixed_rank1", d, n, 1, "0.7", est, "score-oracle relative error", reps)]

    # Experiment 2: equal spikes; isolate the exact sequential-NRM target.
    for d in [250, 500, 1000, 2000]:
        n = 20; r = n - 1; lam = d ** 0.9
        res = simulate(d, n, np.array([lam, lam]), reps, np.random.default_rng(next(children)))
        seq_ratio = res["nrm"][:, 0] / res["score"][:, 0]
        target_ratio = 1.0 - res["score"][:, 1] / ((r - 1) * res["score"][:, 0])
        block_ratio = res["bnrm"][:, 0] / res["score"][:, 0]
        rows += [summaries(seq_ratio, "two_equal_spikes", d, n, 2, "0.9", "sequential NRM", "estimate / score oracle", reps),
                 summaries(target_ratio, "two_equal_spikes", d, n, 2, "0.9", "theoretical target", "estimate / score oracle", reps),
                 summaries(block_ratio, "two_equal_spikes", d, n, 2, "0.9", "bNRM", "estimate / score oracle", reps),
                 summaries(np.abs(seq_ratio-target_ratio), "two_equal_spikes", d, n, 2, "0.9", "sequential NRM", "absolute target discrepancy", reps),
                 summaries(np.abs(block_ratio-1), "two_equal_spikes", d, n, 2, "0.9", "bNRM", "absolute target discrepancy", reps)]

    # Experiment 3: eigenvectors require the trace scale, unlike eigenvalues.
    for alpha in [0.7, 1.2]:
        for d in [250, 500, 1000, 2000]:
            n = 20; lam = d ** alpha
            res = simulate(d, n, np.array([lam]), reps, np.random.default_rng(next(children)), need_direction=True)
            rows.append(summaries(res["sin_angle"], "direction_phase", d, n, 1, str(alpha), "PCA direction", "sine angle", reps))

    # Experiment 4: J-HDLSS sequence n ~= sqrt(d), rank one.
    for d in [400, 900, 1600, 2500]:
        n = int(round(math.sqrt(d))); lam = d ** 0.7
        res = simulate(d, n, np.array([lam]), reps, np.random.default_rng(next(children)))
        poperr = np.abs(res["bnrm"][:, 0] - lam) / lam
        orerr = np.abs(res["bnrm"][:, 0] - res["score"][:, 0]) / np.maximum(res["score"][:, 0], 1e-12)
        rows += [summaries(poperr, "jointhdlss", d, n, 1, "0.7", "bNRM", "population relative error", reps),
                 summaries(orerr, "jointhdlss", d, n, 1, "0.7", "bNRM", "score-oracle relative error", reps)]

    df = pd.DataFrame(rows).sort_values(["experiment", "metric", "estimator", "d"])
    df.to_csv(outdir / "simulation_summary.csv", index=False)
    plot_results(df, figdir)
    manifest = {
        "master_seed": MASTER_SEED, "replicates_per_configuration": reps,
        "numpy_version": np.__version__, "pandas_version": pd.__version__,
        "models": "Gaussian independent factors with isotropic unit bulk",
        "dimensions": {"fixed": [250,500,1000,2000], "joint": [400,900,1600,2500]},
        "notes": "Monte Carlo SE is sample SD divided by sqrt(replicates). No post-selection or discarded draws."
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(df.to_string(index=False, max_rows=200))

if __name__ == "__main__":
    main()
