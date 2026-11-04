"""BAO / RSD tool family: fetch, fit, plot spectroscopic BAO data.

Datasets: DESI DR2/DR1, eBOSS DR16, BOSS DR12, 6dFGS, SDSS MGS.
Data from CobayaSampler/bao_data (bundled in data/bao/).
"""
from __future__ import annotations

import json
from typing import Annotated

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, ok, resolve_outdir, mpl, save_fig

__all__ = [
    "fetch_bao_measurements",
    "plot_bao_distances",
    "fit_bao_cosmology",
    "plot_bao_contours",
    "plot_fsigma8",
    "compare_bao_surveys",
]

DISPATCH_KERNELS = {
    "fit_bao": {
        "function": "bao.kernels.fit_bao",
        "args_shape": "z_arr, obs_arr, data_vec, cov_mat (lists); fit_w0wa, Om_range, H0rd_range, cmb_prior, sn_prior",
        "pip_deps": ["numpy", "scipy"],
        "duration_hint_s": 600,
        "returns": "best_fit {Om, H0rd, [w0, wa]}, chi2_min, chi2_grid, marginalised",
    },
}

BAO_DATASETS = [
    "desi_dr2", "desi_dr1", "eboss_dr16", "eboss_dr16_qso",
    "boss_dr12", "boss_dr12_bao", "sixdf", "sdss_mgs",
]


def _load(dataset: str):
    from .io import load_gaussian
    return load_gaussian(dataset)


def _planck_lcdm(z, quantity):
    """Planck 2018 LCDM prediction for BAO observables."""
    from .kernels import _Ez, _comoving_integral, C_KM_S
    Om = 0.3153
    h_rd = 0.6736 * 147.09  # h*rd in Mpc = 99.04
    c100 = C_KM_S / 100.0
    if quantity == "DM_over_rs":
        return c100 / h_rd * _comoving_integral(z, Om)
    elif quantity == "DH_over_rs":
        return c100 / (h_rd * _Ez(z, Om))
    elif quantity == "DV_over_rs":
        dm_rd = c100 / h_rd * _comoving_integral(z, Om)
        dh_rd = c100 / (h_rd * _Ez(z, Om))
        return (z * dm_rd ** 2 * dh_rd) ** (1.0 / 3.0)
    elif quantity == "rs_over_DV":
        dm_rd = c100 / h_rd * _comoving_integral(z, Om)
        dh_rd = c100 / (h_rd * _Ez(z, Om))
        dv_rd = (z * dm_rd ** 2 * dh_rd) ** (1.0 / 3.0)
        return 1.0 / dv_rd
    return None


@validate_call
def fetch_bao_measurements(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    dataset: Annotated[str, Field(description="BAO dataset key: desi_dr2, desi_dr1, eboss_dr16, eboss_dr16_qso, boss_dr12, boss_dr12_bao, sixdf, sdss_mgs")] = "desi_dr2",
) -> ArtifactResult:
    """Fetch BAO distance measurements for a survey and write to CSV.

    Use this first to get the data. Returns redshifts, DM/rd, DH/rd, DV/rd
    values, errors (diagonal of covariance), and the full covariance as a
    separate file.
    """
    if dataset not in BAO_DATASETS:
        raise ValueError(f"Unknown BAO dataset {dataset!r}. Choose from: {BAO_DATASETS}")
    out = resolve_outdir(output_dir)
    d = _load(dataset)
    import pandas as pd
    rows = []
    diag = np.sqrt(np.diag(d["cov"]))
    for i, (z, v, o) in enumerate(zip(d["z"], d["value"], d["observable"], strict=False)):
        rows.append({"z": z, "value": v, "error": diag[i], "observable": o})
    df = pd.DataFrame(rows)
    csv_path = out / f"bao_{dataset}.csv"
    df.to_csv(csv_path, index=False)
    cov_path = out / f"bao_{dataset}_cov.csv"
    np.savetxt(cov_path, d["cov"], delimiter=",")
    return ok(
        f"{d['label']}: {len(rows)} measurements written.",
        [csv_path, cov_path],
        dataset=dataset, label=d["label"], n_points=len(rows),
        columns=list(df.columns),
    )


@validate_call
def plot_bao_distances(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    datasets: Annotated[str, Field(description="Comma-separated BAO dataset keys (default: all)")] = "desi_dr2,desi_dr1,boss_dr12_bao,sixdf,sdss_mgs",
) -> ArtifactResult:
    """Plot Hubble-diagram-style D/rd/sqrt(z) vs z for BAO surveys with Planck LCDM curve.

    Use after fetch_bao_measurements to visualise how BAO distances compare
    across surveys and against the standard cosmological model.
    """
    plt = mpl()
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ds_list = [s.strip() for s in datasets.split(",")]

    colors = {"desi_dr2": "C0", "desi_dr1": "C1", "eboss_dr16": "C2",
              "eboss_dr16_qso": "C3", "boss_dr12": "C4", "boss_dr12_bao": "C5",
              "sixdf": "C6", "sdss_mgs": "C7"}
    markers = {"DM_over_rs": "o", "DH_over_rs": "s", "DV_over_rs": "D", "rs_over_DV": "^"}

    for ds in ds_list:
        d = _load(ds)
        diag = np.sqrt(np.diag(d["cov"]))
        for i, (z, v, o) in enumerate(zip(d["z"], d["value"], d["observable"], strict=False)):
            if o == "f_sigma8":
                continue
            if o == "rs_over_DV":
                v_plot = 1.0 / v
                e_plot = diag[i] / v ** 2
                o_label = "DV_over_rs"
            else:
                v_plot = v
                e_plot = diag[i]
                o_label = o
            norm = v_plot / np.sqrt(z) if z > 0 else v_plot
            enorm = e_plot / np.sqrt(z) if z > 0 else e_plot
            label = f"{d['label']} ({o_label.replace('_over_rs', '/rd')})" if i == 0 or True else None
            ax.errorbar(z, norm, yerr=enorm, fmt=markers.get(o_label, "o"),
                        color=colors.get(ds, "gray"), ms=5, capsize=2,
                        label=label if i < 3 else None)

    # Planck LCDM curve
    zz = np.linspace(0.05, 2.5, 200)
    for obs, ls, lbl in [("DM_over_rs", "-", "DM/rd"), ("DH_over_rs", "--", "DH/rd"),
                          ("DV_over_rs", ":", "DV/rd")]:
        yy = np.array([_planck_lcdm(z, obs) for z in zz])
        ax.plot(zz, yy / np.sqrt(zz), ls, color="k", lw=1, alpha=0.6,
                label=f"Planck LCDM {lbl}" if obs == "DM_over_rs" else None)

    ax.set_xlabel("Redshift z")
    ax.set_ylabel("D / rd / sqrt(z)")
    ax.set_title("BAO Distance Measurements")
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    ax.set_xlim(0, 2.6)
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "bao_distances")
    return ok("BAO distance plot.", [path])


@validate_call
def fit_bao_cosmology(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    dataset: Annotated[str, Field(description="BAO dataset key")] = "desi_dr2",
    fit_w0wa: Annotated[bool, Field(description="If True, fit w0waCDM instead of LCDM")] = False,
    cmb_prior: Annotated[str | None, Field(description='JSON: {"Om":[mean,sigma],"H0rd":[mean,sigma]} or null')] = None,
    sn_prior: Annotated[str | None, Field(description='JSON: {"w0":[mean,sigma],"wa":[mean,sigma]} or null')] = None,
    n_grid: Annotated[int, Field(ge=20, le=200, description="Grid resolution per parameter")] = 80,
) -> ArtifactResult:
    """Fit BAO data for Om, H0*rd (LCDM) or Om, H0*rd, w0, wa (w0waCDM).

    This is a HEAVY tool (chi-squared grid scan). Can be dispatched to HPC.
    Optional Gaussian priors from CMB or supernovae can be applied.
    Sanity: DESI DR2 BAO-only LCDM should give Om ~ 0.297, H0*rd ~ 101.5.
    """
    from mcp_server.dispatch import remote_site, run_kernel

    d = _load(dataset)
    # filter to BAO-only observables for fitting
    bao_mask = [o in ("DM_over_rs", "DH_over_rs", "DV_over_rs") for o in d["observable"]]
    z_fit = d["z"][bao_mask].tolist()
    obs_fit = [o for o, m in zip(d["observable"], bao_mask, strict=False) if m]
    val_fit = d["value"][bao_mask].tolist()
    cov_fit = d["cov"][np.ix_(bao_mask, bao_mask)].tolist()

    cmb_p = json.loads(cmb_prior) if cmb_prior else None
    sn_p = json.loads(sn_prior) if sn_prior else None

    kwargs = dict(
        z_arr=z_fit, obs_arr=obs_fit, data_vec=val_fit, cov_mat=cov_fit,
        fit_w0wa=fit_w0wa, n_grid_Om=n_grid, n_grid_H0rd=n_grid,
        cmb_prior=cmb_p, sn_prior=sn_p,
    )
    if fit_w0wa:
        kwargs["n_grid_w"] = max(20, n_grid // 3)

    site = remote_site()
    if site:
        res = run_kernel("bao.kernels.fit_bao", kwargs,
                         pip_deps=["numpy", "scipy"], duration=600)
        payload = res["result"]
    else:
        from .kernels import fit_bao
        payload = fit_bao(**kwargs)

    out = resolve_outdir(output_dir)
    result_path = out / f"bao_fit_{dataset}.json"
    # save without the full grid to keep JSON small
    summary = {k: v for k, v in payload.items() if k != "chi2_grid"}
    summary.pop("Om_grid", None)
    summary.pop("H0rd_grid", None)
    result_path.write_text(json.dumps(summary, indent=2))

    # save full grid as npz
    grid_path = out / f"bao_fit_{dataset}_grid.npz"
    np.savez_compressed(grid_path,
                        Om_grid=np.array(payload["Om_grid"]),
                        H0rd_grid=np.array(payload["H0rd_grid"]),
                        chi2_grid=np.array(payload["chi2_grid"]))

    bf = payload.get("best_fit", {})
    marg = payload.get("marginalised", {})
    msg = (f"{payload['model']} fit: Om={bf.get('Om', '?')}, H0*rd={bf.get('H0rd', '?')} "
           f"(chi2/dof={payload['chi2_min']:.1f}/{payload['ndof']})")
    if marg:
        msg += f"; marginalised Om={marg.get('Om','?')} +/- {marg.get('Om_sigma','?')}"

    return ok(msg, [result_path, grid_path], **summary)

fit_bao_cosmology.weight = "heavy"


@validate_call
def plot_bao_contours(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    grid_file: Annotated[str, Field(description="Path to bao_fit_*_grid.npz from fit_bao_cosmology")],
    dataset_label: Annotated[str, Field(description="Label for the plot title")] = "BAO",
) -> ArtifactResult:
    """Plot Om-H0*rd contour from a BAO fit grid.

    Use after fit_bao_cosmology. Reads the npz grid file and shows
    1-sigma and 2-sigma contours.
    """
    plt = mpl()
    from pathlib import Path
    data = np.load(Path(grid_file))
    Om = data["Om_grid"]
    H0rd = data["H0rd_grid"]
    chi2 = data["chi2_grid"]
    dchi2 = chi2 - chi2.min()

    fig, ax = plt.subplots(figsize=(7, 5))
    levels = [2.30, 6.18]  # 1-sigma, 2-sigma for 2 params
    cs = ax.contour(H0rd, Om, dchi2, levels=levels, colors=["C0", "C1"])
    ax.clabel(cs, fmt={2.30: r"$1\sigma$", 6.18: r"$2\sigma$"})
    ax.contourf(H0rd, Om, dchi2, levels=[0, 2.30, 6.18, 100],
                colors=["C0", "C1", "white"], alpha=0.2)
    imin = np.unravel_index(np.argmin(chi2), chi2.shape)
    ax.plot(H0rd[imin[1]], Om[imin[0]], "k+", ms=12, mew=2)
    ax.set_xlabel(r"$H_0 r_d$ [km/s]")
    ax.set_ylabel(r"$\Omega_m$")
    ax.set_title(f"{dataset_label} BAO Constraints")
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "bao_contours")
    return ok("BAO contour plot.", [path])


@validate_call
def plot_fsigma8(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    datasets: Annotated[str, Field(description="Comma-separated dataset keys with f*sigma8 data")] = "eboss_dr16,boss_dr12",
) -> ArtifactResult:
    """Plot f*sigma8(z) measurements vs Planck LCDM prediction.

    Use to visualise the growth-rate measurements from RSD analyses.
    Datasets must contain f_sigma8 observables (eBOSS DR16, BOSS DR12).
    """
    plt = mpl()
    fig, ax = plt.subplots(figsize=(8, 5))
    ds_list = [s.strip() for s in datasets.split(",")]
    colors = {"eboss_dr16": "C0", "eboss_dr16_qso": "C1", "boss_dr12": "C2"}

    for ds in ds_list:
        d = _load(ds)
        diag = np.sqrt(np.diag(d["cov"]))
        for i, (z, v, o) in enumerate(zip(d["z"], d["value"], d["observable"], strict=False)):
            if o != "f_sigma8":
                continue
            ax.errorbar(z, v, yerr=diag[i], fmt="o", color=colors.get(ds, "gray"),
                        ms=6, capsize=3, label=d["label"] if i == 2 or i == 0 else None)

    # Planck LCDM f*sigma8 approximation
    zz = np.linspace(0.1, 2.0, 100)
    Om = 0.3153
    s8 = 0.8111
    fs8 = s8 * Om ** 0.55 * (1 + zz) ** (-0.5) * np.exp(-0.15 * zz)  # rough approx
    ax.plot(zz, fs8, "k-", lw=1.5, alpha=0.6, label="Planck LCDM (approx)")
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title("Growth Rate Measurements")
    ax.legend()
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "fsigma8")
    return ok("f*sigma8 plot.", [path])


@validate_call
def compare_bao_surveys(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    datasets: Annotated[str, Field(description="Comma-separated BAO dataset keys")] = "desi_dr2,desi_dr1,boss_dr12_bao",
) -> ArtifactResult:
    """Compare BAO measurements across surveys: residuals relative to Planck LCDM.

    Plots (data - Planck) / error for each data point, grouped by survey.
    Useful for spotting tensions between surveys.
    """
    plt = mpl()
    fig, ax = plt.subplots(figsize=(10, 5))
    ds_list = [s.strip() for s in datasets.split(",")]
    colors = {"desi_dr2": "C0", "desi_dr1": "C1", "eboss_dr16": "C2",
              "boss_dr12_bao": "C3", "boss_dr12": "C4", "sixdf": "C5", "sdss_mgs": "C6"}
    offset = 0
    tick_labels = []
    tick_pos = []

    for ds in ds_list:
        d = _load(ds)
        diag = np.sqrt(np.diag(d["cov"]))
        for i, (z, v, o) in enumerate(zip(d["z"], d["value"], d["observable"], strict=False)):
            if o == "f_sigma8":
                continue
            if o == "rs_over_DV":
                planck = _planck_lcdm(z, o) if _planck_lcdm(z, o) else 1.0 / _planck_lcdm(z, "DV_over_rs")
            else:
                planck = _planck_lcdm(z, o)
            if planck is None:
                continue
            residual = (v - planck) / diag[i]
            ax.barh(offset, residual, color=colors.get(ds, "gray"), height=0.6, alpha=0.7)
            tick_labels.append(f"{d['label']} z={z:.2f} {o.replace('_over_rs','/rd')}")
            tick_pos.append(offset)
            offset += 1

    ax.axvline(0, color="k", lw=0.8)
    ax.axvline(-1, color="gray", ls="--", lw=0.5)
    ax.axvline(1, color="gray", ls="--", lw=0.5)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(tick_labels, fontsize=6)
    ax.set_xlabel("(Data - Planck LCDM) / sigma")
    ax.set_title("BAO Survey Comparison vs Planck LCDM")
    fig.tight_layout()
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "bao_comparison")
    return ok("BAO comparison plot.", [path])
