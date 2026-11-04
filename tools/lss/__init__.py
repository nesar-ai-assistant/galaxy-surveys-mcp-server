"""LSS clustering catalog tool family: fetch catalogs, measure multipoles, plot.

Tier 1: fetch_galaxy_lss_catalog (raw galaxy + random catalogs)
Tier 2: measure_correlation_multipoles, compute_nz_histogram,
        fetch_published_multipoles, plot_clustering_multipoles
"""
from __future__ import annotations

from typing import Annotated

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, ok, resolve_outdir, mpl, save_fig, tier

__all__ = [
    "fetch_galaxy_lss_catalog",
    "measure_correlation_multipoles",
    "compute_nz_histogram",
    "fetch_published_multipoles",
    "plot_clustering_multipoles",
]

DISPATCH_KERNELS = {
    "measure_xi_multipoles": {
        "function": "lss.kernels.measure_xi_multipoles",
        "args_shape": "data_ra/dec/z/weight, rand_ra/dec/z/weight (lists); s_bins, cosmology",
        "pip_deps": ["numpy", "scipy", "astropy", "Corrfunc"],
        "duration_hint_s": 900,
        "returns": "s_mid, xi0, xi2, n_data, n_rand, used_corrfunc",
    },
}

_CATALOG_KEYS = [
    "eboss_dr16_lrg_ngc", "eboss_dr16_qso_ngc", "eboss_dr16_elg_ngc",
    "boss_dr12_cmass_north", "boss_dr12_lowz_north",
    "desi_dr1_lrg_ngc", "desi_dr1_elg_ngc", "desi_dr1_qso_ngc",
]

_PUBLISHED_KEYS = [
    "boss_dr12_pk_z1_ngc", "boss_dr12_pk_z3_ngc",
    "eboss_dr16_lrg_xi", "eboss_dr16_qso_xi",
]


@tier(1)
@validate_call
def fetch_galaxy_lss_catalog(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    catalog: Annotated[str, Field(description=(
        "LSS catalog key: eboss_dr16_lrg_ngc, eboss_dr16_qso_ngc, "
        "eboss_dr16_elg_ngc, boss_dr12_cmass_north, boss_dr12_lowz_north, "
        "desi_dr1_lrg_ngc, desi_dr1_elg_ngc, desi_dr1_qso_ngc"
    ))] = "eboss_dr16_lrg_ngc",
    kind: Annotated[str, Field(description="'data' or 'random'")] = "data",
    max_rows: Annotated[int, Field(ge=100, le=500000, description="Max rows to download")] = 50000,
    z_min: Annotated[float | None, Field(description="Min redshift filter (null = catalog default)")] = None,
    z_max: Annotated[float | None, Field(description="Max redshift filter (null = catalog default)")] = None,
    random_fraction: Annotated[float, Field(ge=0.001, le=1.0, description="Fraction of rows to keep (subsampling)")] = 1.0,
) -> ArtifactResult:
    """Fetch a subset of a galaxy clustering catalog (data or randoms) from BOSS/eBOSS/DESI.

    Downloads on demand via HTTP Range (partial) for uncompressed FITS or streaming for
    gzipped FITS. Returns RA, DEC, Z, weights as CSV. Use this to get raw tier-1 data
    for correlation function measurement.

    Catalogs include eBOSS DR16 LRG/QSO/ELG (NGC), BOSS DR12 CMASS/LOWZ (North),
    and DESI DR1 LRG/ELG/QSO (NGC). Each comes with systematic + FKP weights.
    """
    if catalog not in _CATALOG_KEYS:
        raise ValueError(f"Unknown catalog {catalog!r}. Choose from: {_CATALOG_KEYS}")
    if kind not in ("data", "random"):
        raise ValueError("kind must be 'data' or 'random'")

    from .io import download_lss_subset
    out = resolve_outdir(output_dir)
    csv_path, meta = download_lss_subset(
        catalog, kind=kind, max_rows=max_rows,
        z_min=z_min, z_max=z_max, random_fraction=random_fraction,
        output_dir=out,
    )
    return ok(
        f"{meta['survey']} {kind}: {meta['rows_after_filter']} rows "
        f"(from {meta['total_rows_remote']} remote, z=[{meta['z_range'][0]:.2f},{meta['z_range'][1]:.2f}]).",
        [csv_path],
        **meta,
    )


@tier(2)
@validate_call
def measure_correlation_multipoles(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    data_file: Annotated[str, Field(description="Path to data catalog CSV (from fetch_galaxy_lss_catalog)")],
    random_file: Annotated[str, Field(description="Path to random catalog CSV")],
    s_max: Annotated[float, Field(ge=50, le=300, description="Max separation in Mpc/h")] = 200.0,
    n_bins: Annotated[int, Field(ge=10, le=100, description="Number of separation bins")] = 40,
) -> ArtifactResult:
    """Measure correlation function multipoles xi_0(s), xi_2(s) from data+random catalogs.

    This is a HEAVY tier-1 to tier-2 reduction. Uses Landy-Szalay estimator with FKP +
    systematic weights. Uses Corrfunc if available, else scipy cKDTree pair-counting.

    Sanity check: for BOSS CMASS, expect BAO bump near s ~ 100-110 Mpc/h in xi_0,
    and xi_0(s=20 Mpc/h) of order 0.1-0.3. Can be dispatched to HPC.
    """
    import pandas as pd
    from pathlib import Path
    from mcp_server.dispatch import remote_site, run_kernel

    data = pd.read_csv(Path(data_file))
    rand = pd.read_csv(Path(random_file))

    # Use WEIGHT_TOTAL if available, else 1.0
    dw = data["WEIGHT_TOTAL"].values if "WEIGHT_TOTAL" in data.columns else np.ones(len(data))
    rw = rand["WEIGHT_TOTAL"].values if "WEIGHT_TOTAL" in rand.columns else np.ones(len(rand))

    kwargs = dict(
        data_ra=data["RA"].tolist(), data_dec=data["DEC"].tolist(),
        data_z=data["Z"].tolist(), data_weight=dw.tolist(),
        rand_ra=rand["RA"].tolist(), rand_dec=rand["DEC"].tolist(),
        rand_z=rand["Z"].tolist(), rand_weight=rw.tolist(),
        s_bins_min=0.0, s_bins_max=s_max, n_s_bins=n_bins,
    )

    site = remote_site()
    if site:
        res = run_kernel("lss.kernels.measure_xi_multipoles", kwargs,
                         pip_deps=["numpy", "scipy", "astropy", "Corrfunc"],
                         duration=900)
        payload = res["result"]
    else:
        from .kernels import measure_xi_multipoles
        payload = measure_xi_multipoles(**kwargs)

    out = resolve_outdir(output_dir)
    xi_data = pd.DataFrame({
        "s_mid": payload["s_mid"],
        "xi0": payload["xi0"],
        "xi2": payload["xi2"],
    })
    csv_path = out / "xi_multipoles.csv"
    xi_data.to_csv(csv_path, index=False)

    return ok(
        f"Measured xi_0, xi_2 in {n_bins} bins to s_max={s_max} Mpc/h. "
        f"N_data={payload['n_data']}, N_rand={payload['n_rand']}. "
        f"Corrfunc={'yes' if payload['used_corrfunc'] else 'no (scipy fallback)'}.",
        [csv_path],
        n_data=payload["n_data"], n_rand=payload["n_rand"],
        used_corrfunc=payload["used_corrfunc"],
        s_range=[0.0, s_max], n_bins=n_bins,
    )

measure_correlation_multipoles.weight = "heavy"


@tier(2)
@validate_call
def compute_nz_histogram(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    catalog_file: Annotated[str, Field(description="Path to catalog CSV (from fetch_galaxy_lss_catalog)")],
    n_bins: Annotated[int, Field(ge=10, le=200, description="Number of redshift bins")] = 50,
    catalog_key: Annotated[str | None, Field(description="If a DESI catalog key, also fetch the official n(z)")] = None,
) -> ArtifactResult:
    """Compute a redshift histogram n(z) from a galaxy catalog.

    Tier-2 summary statistic. Weights are applied if WEIGHT_TOTAL column exists.
    If catalog_key is a DESI catalog, also fetches the official n(z) for comparison.
    """
    import pandas as pd
    from pathlib import Path

    df = pd.read_csv(Path(catalog_file))
    z = df["Z"].values
    w = df["WEIGHT_TOTAL"].values if "WEIGHT_TOTAL" in df.columns else np.ones(len(df))

    counts, edges = np.histogram(z, bins=n_bins, weights=w)
    z_mid = 0.5 * (edges[:-1] + edges[1:])
    dz = edges[1] - edges[0]

    nz_df = pd.DataFrame({"z_mid": z_mid, "z_low": edges[:-1], "z_high": edges[1:],
                           "n_weighted": counts, "n_density": counts / dz})

    out = resolve_outdir(output_dir)
    csv_path = out / "nz_histogram.csv"
    nz_df.to_csv(csv_path, index=False)

    files = [csv_path]
    # Optionally fetch official DESI n(z)
    official_meta = {}
    if catalog_key and catalog_key.startswith("desi_"):
        from .io import fetch_desi_nz
        nz_path, nz_meta = fetch_desi_nz(catalog_key, out)
        if nz_path:
            files.append(nz_path)
            official_meta = {"official_nz": str(nz_path), **nz_meta}

    return ok(
        f"n(z) histogram: {n_bins} bins, z=[{z.min():.3f}, {z.max():.3f}], "
        f"{len(df)} galaxies.",
        files, n_bins=n_bins, z_range=[float(z.min()), float(z.max())],
        total_galaxies=len(df), **official_meta,
    )


@tier(2)
@validate_call
def fetch_published_multipoles(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    dataset: Annotated[str, Field(description=(
        "Published multipole key: boss_dr12_pk_z1_ngc, boss_dr12_pk_z3_ngc, "
        "eboss_dr16_lrg_xi, eboss_dr16_qso_xi"
    ))] = "boss_dr12_pk_z1_ngc",
) -> ArtifactResult:
    """Fetch released P(k) or xi(s) multipoles and covariances from BOSS DR12 or eBOSS DR16.

    These are tier-2 summary statistics measured by the survey teams and published
    alongside the BAO/RSD papers. Use for comparison with your own measurements
    from tier-1 catalogs, or as direct inputs to cosmological fits.
    """
    if dataset not in _PUBLISHED_KEYS:
        raise ValueError(f"Unknown dataset {dataset!r}. Choose from: {_PUBLISHED_KEYS}")

    from .io import PUBLISHED_MULTIPOLES
    out = resolve_outdir(output_dir)
    entry = PUBLISHED_MULTIPOLES[dataset]

    if entry["type"] == "pk":
        from .io import fetch_published_pk_beutler
        csv_path, meta = fetch_published_pk_beutler(dataset, out)
    else:
        from .io import fetch_published_xi
        csv_path, meta = fetch_published_xi(dataset, out)

    return ok(
        f"Published {entry['type']} multipoles: {entry['survey']}. "
        f"{meta.get('n_bins', meta.get('n_k_bins', '?'))} bins. "
        f"Citation: {entry['citation']}.",
        [csv_path], **meta,
    )


@tier(2)
@validate_call
def plot_clustering_multipoles(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    data_file: Annotated[str, Field(description="Path to multipole CSV (from measure_correlation_multipoles or fetch_published_multipoles)")],
    label: Annotated[str, Field(description="Label for the plot legend")] = "measured",
    published_file: Annotated[str | None, Field(description="Optional: path to published multipoles for comparison")] = None,
    published_label: Annotated[str | None, Field(description="Label for published data")] = None,
    plot_type: Annotated[str, Field(description="'xi' for xi(s) or 'pk' for P(k)")] = "xi",
) -> ArtifactResult:
    """Plot clustering multipoles: xi_0, xi_2 vs s or P_0, P_2 vs k.

    Optionally overlay published multipoles for comparison. Creates a two-panel
    plot showing monopole (top) and quadrupole (bottom).
    """
    import pandas as pd
    from pathlib import Path

    plt = mpl()
    df = pd.read_csv(Path(data_file))

    fig, (ax0, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    if plot_type == "xi":
        x_col = "s_mid" if "s_mid" in df.columns else "s"
        x = df[x_col].values
        xlabel = r"$s$ [Mpc/$h$]"

        if "xi0" in df.columns:
            ax0.plot(x, x**2 * df["xi0"].values, "o-", ms=4, label=label)
            ax0.set_ylabel(r"$s^2 \xi_0(s)$")
        if "xi2" in df.columns:
            ax2.plot(x, x**2 * df["xi2"].values, "s-", ms=4, label=label)
            ax2.set_ylabel(r"$s^2 \xi_2(s)$")

        if published_file:
            pf = pd.read_csv(Path(published_file))
            px_col = "s_mid" if "s_mid" in pf.columns else "s"
            px = pf[px_col].values
            plbl = published_label or "published"
            if "xi0" in pf.columns:
                err0 = pf["xi0_err"].values if "xi0_err" in pf.columns else None
                ax0.errorbar(px, px**2 * pf["xi0"].values,
                            yerr=px**2 * err0 if err0 is not None else None,
                            fmt="^", ms=4, alpha=0.7, label=plbl)
            if "xi2" in pf.columns:
                err2 = pf["xi2_err"].values if "xi2_err" in pf.columns else None
                ax2.errorbar(px, px**2 * pf["xi2"].values,
                            yerr=px**2 * err2 if err2 is not None else None,
                            fmt="^", ms=4, alpha=0.7, label=plbl)

    elif plot_type == "pk":
        x_col = "k"
        x = df[x_col].values
        xlabel = r"$k$ [$h$/Mpc]"

        if "P0" in df.columns:
            err0 = df["P0_err"].values if "P0_err" in df.columns else None
            ax0.errorbar(x, df["P0"].values, yerr=err0, fmt="o-", ms=4, label=label)
            ax0.set_ylabel(r"$P_0(k)$ [(Mpc/$h$)$^3$]")
        if "P2" in df.columns:
            err2 = df["P2_err"].values if "P2_err" in df.columns else None
            ax2.errorbar(x, df["P2"].values, yerr=err2, fmt="s-", ms=4, label=label)
            ax2.set_ylabel(r"$P_2(k)$ [(Mpc/$h$)$^3$]")

        if published_file:
            pf = pd.read_csv(Path(published_file))
            plbl = published_label or "published"
            if "P0" in pf.columns:
                ax0.errorbar(pf["k"].values, pf["P0"].values,
                            yerr=pf["P0_err"].values if "P0_err" in pf.columns else None,
                            fmt="^", ms=4, alpha=0.7, label=plbl)
            if "P2" in pf.columns:
                ax2.errorbar(pf["k"].values, pf["P2"].values,
                            yerr=pf["P2_err"].values if "P2_err" in pf.columns else None,
                            fmt="^", ms=4, alpha=0.7, label=plbl)

    ax0.set_title("Clustering Multipoles")
    ax0.legend(fontsize=9)
    ax2.legend(fontsize=9)
    ax2.set_xlabel(xlabel)
    ax0.axhline(0, color="gray", lw=0.5, ls="--")
    ax2.axhline(0, color="gray", lw=0.5, ls="--")
    fig.tight_layout()

    out = resolve_outdir(output_dir)
    stem = "clustering_multipoles_pk" if plot_type == "pk" else "clustering_multipoles_xi"
    path = save_fig(fig, out, stem)
    return ok(f"Clustering multipoles plot ({plot_type}).", [path])
