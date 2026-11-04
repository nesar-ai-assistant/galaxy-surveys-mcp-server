"""SDSS SkyServer tools: SQL queries, galaxy catalogs, photometry, colours.

Ported from the user's sdss_query_v2.py reference code; drops smolagents,
uses the ArtifactResult + file-path contract.
"""
from __future__ import annotations

import hashlib
from typing import Annotated

import numpy as np
from pydantic import Field, validate_call

from ..common import (ArtifactResult, ok, resolve_outdir, mpl, save_fig,
                      data_cache_dir, data_kind)

__all__ = [
    "sdss_sql_query",
    "fetch_sdss_galaxies",
    "compute_sdss_absolute_magnitudes",
    "plot_color_magnitude",
    "plot_sdss_redshift_distribution",
]

DISPATCH_KERNELS: dict = {}


def _cache_dir():
    return data_cache_dir() / "sdss"


def _query_sdss_sql(sql: str, timeout: int = 300, data_release: int = 18):
    """Execute an SDSS SQL query and return a pandas DataFrame (with caching)."""
    import pandas as pd
    cache = _cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    h = hashlib.md5(sql.encode()).hexdigest()
    cache_path = cache / f"query_{h}.csv"
    if cache_path.is_file():
        return pd.read_csv(cache_path)

    from astroquery.sdss import SDSS
    result = SDSS.query_sql(sql, timeout=timeout, data_release=data_release)
    if result is None or len(result) == 0:
        return pd.DataFrame()
    df = result.to_pandas()
    df.to_csv(cache_path, index=False)
    return df


def _filter_invalid(df):
    """Filter SDSS invalid sentinel values: -9999 for mags, -1000 for errors."""
    if df.empty:
        return df.copy()
    out = df.copy()
    for col in df.columns:
        if "Mag" in col and "err" not in col.lower():
            out = out[out[col] > -9999]
        elif "err" in col.lower():
            out = out[out[col] > -1000]
    return out


@data_kind("catalog")
@validate_call
def sdss_sql_query(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    sql: Annotated[str, Field(min_length=10, description="Read-only SELECT SQL query for SDSS SkyServer DR18.")],
    timeout: Annotated[int, Field(ge=30, le=600, description="Query timeout in seconds")] = 300,
    filter_invalid: Annotated[bool, Field(description="Filter SDSS sentinel values")] = True,
) -> ArtifactResult:
    """Execute an arbitrary read-only SQL SELECT on the SDSS SkyServer (DR18).

    Use this for custom SDSS queries. Only SELECT is allowed; the query runs
    via astroquery and results are cached on disk. Results written to CSV.
    Network access required.
    """
    sql_up = sql.strip().upper()
    if not sql_up.startswith("SELECT"):
        raise ValueError("Only SELECT queries are allowed.")
    for forbidden in ("INSERT", "UPDATE", "DELETE", "DROP", "ALTER", "CREATE"):
        if forbidden in sql_up:
            raise ValueError(f"Forbidden keyword {forbidden} in query.")

    df = _query_sdss_sql(sql, timeout=timeout)
    if df.empty:
        return ok("Query returned 0 rows.", [], n_rows=0)
    if filter_invalid:
        df = _filter_invalid(df)
    out = resolve_outdir(output_dir)
    path = out / "sdss_query_result.csv"
    df.to_csv(path, index=False)
    return ok(f"Retrieved {len(df)} rows with {len(df.columns)} columns.",
              [path], n_rows=len(df), columns=list(df.columns))


@data_kind("catalog")
@validate_call
def fetch_sdss_galaxies(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    limit: Annotated[int, Field(ge=10, le=100000, description="Max galaxies to return")] = 1000,
    min_z: Annotated[float, Field(ge=0.0, le=1.0)] = 0.01,
    max_z: Annotated[float, Field(ge=0.001, le=1.0)] = 0.1,
    max_r_mag: Annotated[float, Field(ge=10, le=25)] = 22.0,
    clean_photometry: Annotated[bool, Field(description="Require clean=1")] = True,
) -> ArtifactResult:
    """Fetch a galaxy sample from SDSS with spectroscopic redshifts and ugriz photometry.

    Use this to get a galaxy catalog for colour-magnitude diagrams, redshift
    distributions, or cross-matching with GSWLC. Requires network access.
    Returns CSV with magnitudes, extinction, redshift, clean flag.
    """
    clean_clause = "AND p.clean = 1" if clean_photometry else ""
    sql = f"""
    SELECT TOP {limit}
        s.specObjID, s.bestObjID, s.ra, s.dec, s.z AS redshift, s.class,
        p.objID,
        p.petroMag_u, p.petroMag_g, p.petroMag_r, p.petroMag_i, p.petroMag_z,
        p.modelMag_u, p.modelMag_g, p.modelMag_r, p.modelMag_i, p.modelMag_z,
        p.extinction_u, p.extinction_g, p.extinction_r, p.extinction_i, p.extinction_z,
        p.clean
    FROM SpecObj AS s
    JOIN PhotoObj AS p ON s.bestObjID = p.objID
    WHERE s.class = 'GALAXY'
        AND s.z BETWEEN {min_z} AND {max_z}
        AND p.petroMag_r < {max_r_mag}
        {clean_clause}
    ORDER BY s.z
    """
    df = _query_sdss_sql(sql)
    if df.empty:
        return ok("No galaxies found with these parameters.", [], n_rows=0)
    df = _filter_invalid(df)
    out = resolve_outdir(output_dir)
    path = out / "sdss_galaxies.csv"
    df.to_csv(path, index=False)
    return ok(f"Fetched {len(df)} SDSS galaxies.", [path],
              n_rows=len(df), z_range=[float(df["redshift"].min()), float(df["redshift"].max())],
              columns=list(df.columns))


@data_kind("catalog")
@validate_call
def compute_sdss_absolute_magnitudes(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    input_file: Annotated[str, Field(description="Path to CSV from fetch_sdss_galaxies")],
    correct_extinction: Annotated[bool, Field(description="Apply extinction correction first")] = True,
) -> ArtifactResult:
    """Compute absolute magnitudes and colours from SDSS apparent magnitudes.

    Reads the CSV from fetch_sdss_galaxies, applies extinction correction
    (subtracting extinction_X from each band), computes distance modulus
    using FlatLambdaCDM(H0=70, Om0=0.3), and derives absolute mags and
    colours (u-g, g-r, r-i, i-z, u-r).
    """
    import pandas as pd
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u

    df = pd.read_csv(input_file)
    if df.empty:
        return ok("Input file is empty.", [])

    bands = ["u", "g", "r", "i", "z"]
    # Extinction correction
    if correct_extinction:
        for b in bands:
            mag_col = f"modelMag_{b}"
            ext_col = f"extinction_{b}"
            if mag_col in df.columns and ext_col in df.columns:
                df[f"modelMag_{b}_corr"] = df[mag_col] - df[ext_col]
        suffix = "_corr"
    else:
        suffix = ""

    # Distance modulus
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    valid_z = df["redshift"].values
    valid_z = np.clip(valid_z, 1e-4, None)
    dl = cosmo.luminosity_distance(valid_z).to(u.pc).value
    dm = 5 * np.log10(dl) - 5

    for b in bands:
        col = f"modelMag_{b}{suffix}"
        if col in df.columns:
            df[f"M_{b}"] = df[col].values - dm

    # Colours
    colour_pairs = [("u", "g"), ("g", "r"), ("r", "i"), ("i", "z"), ("u", "r")]
    for b1, b2 in colour_pairs:
        c1 = f"modelMag_{b1}{suffix}"
        c2 = f"modelMag_{b2}{suffix}"
        if c1 in df.columns and c2 in df.columns:
            df[f"{b1}_{b2}"] = df[c1] - df[c2]

    out = resolve_outdir(output_dir)
    path = out / "sdss_absolute_mags.csv"
    df.to_csv(path, index=False)
    return ok(f"Computed absolute magnitudes and colours for {len(df)} galaxies.",
              [path], n_rows=len(df), columns=list(df.columns))


@data_kind("catalog")
@validate_call
def plot_color_magnitude(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    input_file: Annotated[str, Field(description="Path to CSV with absolute mags (from compute_sdss_absolute_magnitudes)")],
) -> ArtifactResult:
    """Plot colour-magnitude diagram (g-r vs M_r) showing the red sequence and blue cloud.

    Use after compute_sdss_absolute_magnitudes. Classic bimodal galaxy
    population diagnostic.
    """
    import pandas as pd
    plt = mpl()

    df = pd.read_csv(input_file)
    if "g_r" not in df.columns or "M_r" not in df.columns:
        raise ValueError("Input file must have 'g_r' and 'M_r' columns (from compute_sdss_absolute_magnitudes).")

    mask = np.isfinite(df["g_r"]) & np.isfinite(df["M_r"])
    gr = df.loc[mask, "g_r"].values
    Mr = df.loc[mask, "M_r"].values

    fig, ax = plt.subplots(figsize=(8, 6))
    h = ax.hexbin(Mr, gr, gridsize=60, cmap="viridis", mincnt=1)
    fig.colorbar(h, ax=ax, label="count")
    ax.set_xlabel(r"$M_r$ (absolute magnitude)")
    ax.set_ylabel(r"$g - r$ (colour)")
    ax.set_title("SDSS Colour-Magnitude Diagram")
    ax.invert_xaxis()
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "color_magnitude")
    return ok(f"CMD plotted for {mask.sum()} galaxies.", [path])


@data_kind("catalog")
@validate_call
def plot_sdss_redshift_distribution(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    input_file: Annotated[str, Field(description="Path to CSV with redshift column (from fetch_sdss_galaxies)")],
    n_bins: Annotated[int, Field(ge=10, le=200)] = 50,
) -> ArtifactResult:
    """Plot the redshift distribution N(z) of an SDSS galaxy sample.

    Use after fetch_sdss_galaxies. Histogram of spectroscopic redshifts.
    """
    import pandas as pd
    plt = mpl()

    df = pd.read_csv(input_file)
    z = df["redshift"].dropna().values

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(z, bins=n_bins, edgecolor="k", alpha=0.7, color="steelblue")
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("N(z)")
    ax.set_title(f"SDSS Galaxy Redshift Distribution (N={len(z):,})")
    ax.axvline(np.median(z), color="red", ls="--", label=f"median z = {np.median(z):.3f}")
    ax.legend()
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "sdss_redshift_dist")
    return ok(f"Redshift distribution for {len(z)} galaxies.", [path],
              n_galaxies=len(z), z_median=round(float(np.median(z)), 4))
