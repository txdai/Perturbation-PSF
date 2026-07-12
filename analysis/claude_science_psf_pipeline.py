"""
Perturbation PSF — reusable kernel-estimation pipeline.

Estimates cell-type-resolved radial response kernels K(r) for barcode-negative
recipient cells around sparse barcode-positive perturbed source cells, versus a
matched safe-harbour (msafe) control, on Spatial Perturb-seq (Stereo-seq) data.

Inputs (checkpoints/):
  typed_cells.parquet      per-cell coords(µm), perturbation, cell_type, niche, depth
  program_scores.parquet   per-cell pathway program scores + control genes
Outputs:
  all_kernels.parquet      K(r) ± CI for every (perturbation, program)
  all_kernel_params.csv    peak_gain, peak_r, integrated, nearfield per (pert, program)

Usage:
  python psf_pipeline.py --rmax 250 --out checkpoints
  # or import: from psf_pipeline import build_neighborhoods, fit_kernel, run_all
"""
import numpy as np, pandas as pd, json, argparse, warnings
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf
from patsy import dmatrix
warnings.filterwarnings("ignore")
trapz = np.trapezoid

KNOTS = [30, 60, 100, 160]
RGRID = np.linspace(2, 250, 80)
SPLINE = f"bs(dist_um, knots={KNOTS}, degree=3, include_intercept=False)"
CONTROL = "msafe"
_Bg = np.asarray(dmatrix(SPLINE, {"dist_um": RGRID}, return_type="dataframe"))


def load(indir="checkpoints"):
    tc = pd.read_parquet(f"{indir}/typed_cells.parquet").reset_index()
    ps = pd.read_parquet(f"{indir}/program_scores.parquet")
    ps.index = ps.index.astype(str)
    return tc, ps


def build_neighborhoods(tc, sources, rmax=250.0):
    """All barcode-negative recipients within rmax µm of each source cell."""
    recip = tc[tc.perturbation == "Neg"]
    rtrees, ridx = {}, {}
    for sec, sub in recip.groupby("section", observed=True):
        rtrees[sec] = cKDTree(sub[["x_um", "y_um"]].values)
        ridx[sec] = sub.index.values
    rows = []
    for s in sources.itertuples():
        if s.section not in rtrees:
            continue
        tr, idxmap = rtrees[s.section], ridx[s.section]
        for j in tr.query_ball_point([s.x_um, s.y_um], rmax):
            ci = idxmap[j]
            dx = tc.at[ci, "x_um"] - s.x_um; dy = tc.at[ci, "y_um"] - s.y_um
            dist = np.hypot(dx, dy)
            if dist < 1e-6:
                continue
            rows.append((s.global_id, s.perturbation, s.section, ci, dist,
                         np.arctan2(dy, dx), tc.at[ci, "cell_type"], tc.at[ci, "niche"]))
    return pd.DataFrame(rows, columns=["source_id", "perturbation", "section",
                        "recipient_id", "dist_um", "angle", "recip_type", "recip_niche"])


def make_frame(nb, tc, ps, programs):
    """Join recipient program scores + depth onto the neighborhood table."""
    gid = tc["global_id"].astype(str).values
    nb = nb.copy()
    nb["recip_gid"] = gid[nb.recipient_id.values]
    S = ps.reindex(nb.recip_gid.values)[programs].reset_index(drop=True)
    af = pd.concat([nb.reset_index(drop=True), S], axis=1)
    af["recip_nCount"] = ps.reindex(nb.recip_gid.values)["nCount"].values
    af["log_depth"] = np.log1p(af.recip_nCount)
    for c in ("recip_type", "recip_niche", "section"):
        af[c] = af[c].astype(str)
    return af.dropna(subset=programs).reset_index(drop=True)


def _contrast(m):
    names = list(m.params.index)
    L = np.zeros((len(RGRID), len(names))); L[:, names.index("treat")] = 1.0
    for k, n in enumerate([x for x in names if x.startswith("treat:bs")]):
        L[:, names.index(n)] = _Bg[:, k]
    beta = m.params.values; cov = m.cov_params().values
    return L @ beta, np.sqrt(np.einsum("ij,jk,ik->i", L, cov, L))


def fit_kernel(af, pert, program):
    """OLS spline kernel with treat×distance interaction, source-clustered SE.
       Returns K(r) DataFrame with 95% CI."""
    d = af[af.perturbation.isin([pert, CONTROL])].copy()
    d["treat"] = (d.perturbation == pert).astype(int)
    f = (f"{program} ~ {SPLINE} + treat + treat:{SPLINE} "
         f"+ C(recip_niche) + C(section) + C(recip_type) + log_depth")
    m = smf.ols(f, data=d).fit(cov_type="cluster", cov_kwds={"groups": d.source_id})
    K, se = _contrast(m)
    return pd.DataFrame({"perturbation": pert, "program": program, "r": RGRID,
                         "K": K, "se": se, "lo": K - 1.96 * se, "hi": K + 1.96 * se})


def summarize(res):
    K, r, se = res.K.values, res.r.values, res.se.values
    ipk = int(np.argmax(np.abs(K)))
    return dict(perturbation=res.perturbation.iloc[0], program=res.program.iloc[0],
                peak_gain=K[ipk], peak_r=r[ipk], peak_sig=bool(abs(K[ipk]) > 1.96 * se[ipk]),
                integrated=trapz(K, r), nearfield_mean=K[r < 25].mean(), se_at_peak=se[ipk])


def run_all(indir="checkpoints", outdir="checkpoints", rmax=250.0,
            perts=None, programs=None, min_sources=5):
    tc, ps = load(indir)
    programs = programs or list(pd.read_csv(f"{indir}/programs_used.csv").program)
    # single-guide sources per perturbation (exclude Neg/Multiple/control from the pert list)
    counts = tc.perturbation.value_counts()
    guide_perts = [p for p in counts.index if p not in ("Neg", "Multiple", CONTROL)]
    perts = perts or [p for p in guide_perts if counts[p] >= min_sources]
    allsrc = tc[tc.perturbation.isin(perts + [CONTROL])].copy()
    nb = build_neighborhoods(tc, allsrc, rmax)
    af = make_frame(nb, tc, ps, programs)
    allK, summ = [], []
    for pert in perts:
        nsrc = af[af.perturbation == pert].source_id.nunique()
        if nsrc < min_sources:
            continue
        for prog in programs:
            try:
                res = fit_kernel(af, pert, prog)
            except Exception:
                continue
            allK.append(res); summ.append({**summarize(res), "n_sources": nsrc})
    Kdf = pd.concat(allK, ignore_index=True)
    sdf = pd.DataFrame(summ)
    Kdf.to_parquet(f"{outdir}/all_kernels.parquet")
    sdf.to_csv(f"{outdir}/all_kernel_params.csv", index=False)
    return Kdf, sdf


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="checkpoints")
    ap.add_argument("--out", default="checkpoints")
    ap.add_argument("--rmax", type=float, default=250.0)
    ap.add_argument("--min-sources", type=int, default=5)
    a = ap.parse_args()
    Kdf, sdf = run_all(a.indir, a.out, a.rmax, min_sources=a.min_sources)
    print(f"Fit kernels for {sdf.perturbation.nunique()} perturbations × "
          f"{sdf.program.nunique()} programs = {len(sdf)} kernels")
    print(f"Perturbations: {sorted(sdf.perturbation.unique())}")
    print(f"Significant peaks (source-clustered): {int(sdf.peak_sig.sum())}")
