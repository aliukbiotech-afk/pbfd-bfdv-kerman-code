"""Epidemiological statistics for the Kerman PBFD survey (136 birds).

Pure python (numpy / pandas / scipy / statsmodels). Used by the Colab notebook
and runnable locally:  python analysis/scripts/pbfd_stats.py
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
from pbfd_names import species_label, load_species_groups

SEED = 20261103
BAND_ORDER = ["None", "Low", "Moderate", "High"]


# ------------------------------------------------------------------ helpers
def wilson_ci(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, centre - half), min(1.0, centre + half)


def chi2_stat(table: np.ndarray) -> float:
    r, c = table.sum(1, keepdims=True), table.sum(0, keepdims=True)
    e = r @ c / table.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(e > 0, (table - e) ** 2 / e, 0.0)
    return float(x.sum())


def monte_carlo_independence(x: pd.Series, y: pd.Series, n_iter: int = 100_000, seed: int = SEED):
    """Permutation (Monte Carlo) test of independence for an r x c table.

    Conditional on both margins (as Fisher-Freeman-Halton); statistic = Pearson X^2.
    Returns observed X^2, Monte Carlo p (with +1 correction), Cramer's V, table.
    """
    rng = np.random.default_rng(seed)
    xa, ya = pd.Categorical(x), pd.Categorical(y)
    xi, yi = xa.codes, ya.codes
    nr, nc = len(xa.categories), len(ya.categories)

    def tab(a, b):
        t = np.zeros((nr, nc)); np.add.at(t, (a, b), 1); return t

    t_obs = tab(xi, yi)
    s_obs = chi2_stat(t_obs)
    ge = 0
    yy = yi.copy()
    for _ in range(n_iter):
        rng.shuffle(yy)
        if chi2_stat(tab(xi, yy)) >= s_obs - 1e-12:
            ge += 1
    p = (ge + 1) / (n_iter + 1)
    n = t_obs.sum(); k = min(nr, nc)
    v = math.sqrt(s_obs / (n * (k - 1))) if k > 1 else float("nan")
    table = pd.DataFrame(t_obs.astype(int), index=xa.categories, columns=ya.categories)
    return dict(chi2=s_obs, p_mc=p, cramers_v=v, n_iter=n_iter, table=table)


def expected_cells_report(table: pd.DataFrame):
    t = table.values.astype(float)
    e = t.sum(1, keepdims=True) @ t.sum(0, keepdims=True) / t.sum()
    return dict(min_expected=float(e.min()), pct_cells_expected_lt5=float((e < 5).mean() * 100))


# ------------------------------------------------------------------ main analysis
def load_birds(meta_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(meta_dir / "species_mapping.csv")
    groups = load_species_groups(meta_dir)
    df["species_group"] = df["host_excel"].map(dict(zip(groups.host_excel, groups.group_id)))
    df["species_label"] = df["species_group"].map(lambda g: species_label(g, groups))
    df["band_label"] = pd.Categorical(df["band_label"], BAND_ORDER, ordered=True)
    df["sex"] = df["sex"].str.lower().str.strip()
    return df


def prevalence_tables(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    out = {}
    k, n = int(df.positive.sum()), len(df)
    lo, hi = wilson_ci(k, n)
    out["overall"] = pd.DataFrame([dict(group="All birds", n=n, positive=k, prevalence=k / n, ci_low=lo, ci_high=hi)])
    for col, name in [("species_group", "by_species"), ("sex", "by_sex")]:
        rows = []
        for g, sub in df.groupby(col, observed=True):
            k, n = int(sub.positive.sum()), len(sub)
            lo, hi = wilson_ci(k, n)
            rows.append(dict(group=g, label=sub["species_label"].iloc[0] if col == "species_group" else g,
                             n=n, positive=k, prevalence=k / n, ci_low=lo, ci_high=hi,
                             low=int((sub.band_label == "Low").sum()), moderate=int((sub.band_label == "Moderate").sum()),
                             high=int((sub.band_label == "High").sum())))
        out[name] = pd.DataFrame(rows).sort_values("n", ascending=False)
    return out


def association_tests(df: pd.DataFrame, n_iter: int = 100_000) -> dict:
    from scipy import stats
    res = {}
    # sex x positivity (2x2 Fisher exact)
    t = pd.crosstab(df.sex, df.positive)
    odds, p = stats.fisher_exact(t.values)
    res["sex_vs_positivity"] = dict(test="Fisher exact (two-sided)", table=t.to_dict(), odds_ratio=odds, p=p)
    # species x positivity
    mc = monte_carlo_independence(df.species_group, df.positive, n_iter)
    res["species_vs_positivity"] = dict(test="Monte Carlo exact test of independence (X^2 statistic, margins fixed)",
                                        chi2=mc["chi2"], p=mc["p_mc"], cramers_v=mc["cramers_v"],
                                        **expected_cells_report(mc["table"]), table=mc["table"].to_dict())
    # collapsed: groups with n < 10 merged into 'Other'
    sizes = df.species_group.value_counts()
    coll = df.species_group.where(df.species_group.map(sizes) >= 10, "Other(<10)")
    mc = monte_carlo_independence(coll, df.positive, n_iter)
    res["species_collapsed_vs_positivity"] = dict(test="as above, groups n<10 collapsed", chi2=mc["chi2"], p=mc["p_mc"],
                                                  cramers_v=mc["cramers_v"], **expected_cells_report(mc["table"]),
                                                  table=mc["table"].to_dict())
    # replicate the thesis analysis: species x 4-level band (incl. None)
    mc = monte_carlo_independence(df.species_group, df.band_label.astype(str), n_iter)
    asym = chi2_asymptotic(mc["table"])
    res["thesis_replication_species_vs_band4"] = dict(
        test="species x band (None/Low/Moderate/High) - replicates thesis Table 4-1", chi2=mc["chi2"],
        p_asymptotic=asym, p_mc=mc["p_mc"], cramers_v=mc["cramers_v"], **expected_cells_report(mc["table"]))
    # intensity among positives only
    pos = df[df.positive == 1]
    for col in ["species_group", "sex"]:
        mc = monte_carlo_independence(pos[col], pos.band_label.astype(str), n_iter)
        res[f"intensity_among_positives_vs_{col}"] = dict(
            test="Monte Carlo exact test, positives only (Low/Moderate/High)", n=len(pos), chi2=mc["chi2"],
            p=mc["p_mc"], cramers_v=mc["cramers_v"], **expected_cells_report(mc["table"]), table=mc["table"].to_dict())
    kw = stats.kruskal(*[g.pcr_band.values for _, g in pos.groupby("species_group") if len(g)])
    res["intensity_among_positives_kruskal_species"] = dict(H=kw.statistic, p=kw.pvalue)
    mw = stats.mannwhitneyu(pos[pos.sex == "male"].pcr_band, pos[pos.sex == "female"].pcr_band)
    res["intensity_among_positives_mannwhitney_sex"] = dict(U=mw.statistic, p=mw.pvalue)
    return res


def chi2_asymptotic(table: pd.DataFrame) -> float:
    from scipy import stats
    return float(stats.chi2_contingency(table.values, correction=False)[1])


def firth_logit(X: np.ndarray, y: np.ndarray, max_iter: int = 200, tol: float = 1e-8):
    """Firth penalised logistic regression (Heinze & Schemper 2002); returns beta, se."""
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(max_iter):
        eta = X @ beta; p = 1 / (1 + np.exp(-eta)); w = p * (1 - p)
        XtW = X.T * w
        info = XtW @ X
        info_inv = np.linalg.pinv(info)
        h = np.einsum("ij,jk,ik->i", X * w[:, None] ** 0.5, info_inv, X * w[:, None] ** 0.5)
        score = X.T @ (y - p + h * (0.5 - p))
        step = info_inv @ score
        beta_new = beta + step
        if np.max(np.abs(step)) < tol:
            beta = beta_new; break
        beta = beta_new
    eta = X @ beta; p = 1 / (1 + np.exp(-eta)); info = (X.T * (p * (1 - p))) @ X
    se = np.sqrt(np.diag(np.linalg.pinv(info)))
    return beta, se


def regression(df: pd.DataFrame) -> pd.DataFrame:
    """Firth logistic regression: positive ~ species_group (ref = largest) + sex (ref = female)."""
    from scipy import stats
    ref = df.species_group.value_counts().idxmax()
    sp = pd.get_dummies(df.species_group).drop(columns=ref).astype(float)
    X = pd.concat([pd.Series(1.0, index=df.index, name="Intercept"), sp,
                   (df.sex == "male").astype(float).rename("sex[male]")], axis=1)
    beta, se = firth_logit(X.values, df.positive.values.astype(float))
    z = beta / se
    out = pd.DataFrame(dict(term=X.columns, beta=beta, se=se, OR=np.exp(beta),
                            OR_low=np.exp(beta - 1.96 * se), OR_high=np.exp(beta + 1.96 * se),
                            p_wald=2 * stats.norm.sf(np.abs(z))))
    out.attrs["reference_species_group"] = ref
    return out


def run_all(meta_dir: Path, out_dir: Path, n_iter: int = 100_000) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = load_birds(meta_dir)
    prev = prevalence_tables(df)
    for k, t in prev.items():
        t.to_csv(out_dir / f"prevalence_{k}.csv", index=False)
    tests = association_tests(df, n_iter)
    reg = regression(df); reg.to_csv(out_dir / "firth_logistic_regression.csv", index=False)
    tests["firth_regression_reference_group"] = reg.attrs["reference_species_group"]

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        return o
    json.dump(clean(tests), open(out_dir / "association_tests.json", "w"), indent=2)
    df.to_csv(out_dir / "birds_analysis_table.csv", index=False)
    return dict(prevalence=prev, tests=tests, regression=reg, birds=df)


if __name__ == "__main__":
    import sys
    root = Path(__file__).resolve().parents[1]
    n_iter = int(sys.argv[1]) if len(sys.argv) > 1 else 20_000
    r = run_all(root / "metadata", root / "results" / "local_stats", n_iter)
    print(r["prevalence"]["overall"].to_string(index=False))
    print(r["prevalence"]["by_species"][["group", "label", "n", "positive", "prevalence", "ci_low", "ci_high"]].to_string(index=False))
    print(r["prevalence"]["by_sex"][["group", "n", "positive", "prevalence"]].to_string(index=False))
    for k, v in r["tests"].items():
        if isinstance(v, dict):
            print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "table"})
    print(r["regression"].round(3).to_string(index=False))
