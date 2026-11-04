"""GSWLC tools: query, describe, plot, cross-match the GALEX-SDSS-WISE Legacy Catalog.

Ported from the user's extract_gswlc_data.py; drops smolagents.
The GSWLC-A2 catalog (~100 MB) is downloaded on demand from
https://salims.pages.iu.edu/gswlc/GSWLC-A2.dat.gz.
"""
from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Annotated

import numpy as np
from pydantic import Field, validate_call

from ..common import (ArtifactResult, ok, resolve_outdir, mpl, save_fig, data_kind,
                      data_cache_dir, fetch_file, get_cached)

__all__ = [
    "query_gswlc_catalog",
    "describe_gswlc_catalog",
    "plot_star_forming_main_sequence",
    "crossmatch_sdss_gswlc",
]

DISPATCH_KERNELS: dict = {}

GSWLC_URL = "https://salims.pages.iu.edu/gswlc/GSWLC-A2.dat.gz"

COLUMN_NAMES = [
    "ObjID", "GLXID", "plate", "mjd", "fiberID", "ra", "dec", "z",
    "chi2", "logMstar", "sigma_logMstar", "logSFR", "sigma_logSFR",
    "A_FUV", "sigma_A_FUV", "A_B", "sigma_A_B", "A_V", "sigma_A_V",
    "flag_sed", "UV_survey", "flag_uv", "flag_midir", "flag_mgs",
]

COLUMN_DESCRIPTIONS = {
    "ObjID": "SDSS DR10 ObjID",
    "GLXID": "GALEX unique object ID (-99 if no UV detection)",
    "plate": "SDSS spectroscopic plate",
    "mjd": "MJD of spectroscopic observation",
    "fiberID": "SDSS fiber ID",
    "ra": "Right ascension (degrees, J2000)",
    "dec": "Declination (degrees, J2000)",
    "z": "Spectroscopic redshift",
    "chi2": "Reduced chi-squared of SED fit",
    "logMstar": "Log stellar mass (Msun)",
    "sigma_logMstar": "Uncertainty in log stellar mass",
    "logSFR": "Log star formation rate (Msun/yr)",
    "sigma_logSFR": "Uncertainty in log SFR",
    "A_FUV": "FUV attenuation (mag)",
    "sigma_A_FUV": "Uncertainty in A_FUV",
    "A_B": "B-band attenuation (mag)",
    "sigma_A_B": "Uncertainty in A_B",
    "A_V": "V-band attenuation (mag)",
    "sigma_A_V": "Uncertainty in A_V",
    "flag_sed": "SED fit quality flag (0=good)",
    "UV_survey": "GALEX survey depth (1=AIS, 2=MIS, 3=DIS)",
    "flag_uv": "UV detection flag",
    "flag_midir": "Mid-IR detection flag (0=22um, 7=12um only, etc.)",
    "flag_mgs": "SDSS Main Galaxy Sample flag (1=MGS member)",
}


def _ensure_catalog() -> Path:
    """Download GSWLC-A2 if not cached; return path to the .dat file."""
    cache = data_cache_dir() / "gswlc"
    cache.mkdir(parents=True, exist_ok=True)
    dat_path = cache / "GSWLC-A2.dat"
    if dat_path.is_file() and dat_path.stat().st_size > 1_000_000:
        return dat_path
    gz_path = cache / "GSWLC-A2.dat.gz"
    if not gz_path.is_file():
        fetch_file(GSWLC_URL, "GSWLC-A2.dat.gz", subdir="gswlc", timeout=300)
        gz_path = data_cache_dir() / "gswlc" / "GSWLC-A2.dat.gz"
    with gzip.open(gz_path, "rb") as f_in, open(dat_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return dat_path


def _load_catalog():
    """Load the full GSWLC-A2 catalog into a pandas DataFrame (cached)."""
    def _factory():
        import pandas as pd
        path = _ensure_catalog()
        return pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=COLUMN_NAMES)
    return get_cached("gswlc_a2", _factory)


def _apply_filters(df, filters: dict):
    """Apply dict-based filters: min/max, eq, in, gt/lt/gte/lte."""
    for col, crit in filters.items():
        if col not in df.columns:
            raise ValueError(f"Unknown column {col!r}. Valid: {list(df.columns)}")
        if not isinstance(crit, dict):
            raise ValueError(f"Filter for {col} must be a dict, got {type(crit)}")
        if "min" in crit:
            df = df[df[col] >= crit["min"]]
        if "max" in crit:
            df = df[df[col] <= crit["max"]]
        if "eq" in crit:
            df = df[df[col] == crit["eq"]]
        if "in" in crit:
            df = df[df[col].isin(crit["in"])]
        if "gt" in crit:
            df = df[df[col] > crit["gt"]]
        if "lt" in crit:
            df = df[df[col] < crit["lt"]]
        if "gte" in crit:
            df = df[df[col] >= crit["gte"]]
        if "lte" in crit:
            df = df[df[col] <= crit["lte"]]
    return df


@data_kind("catalog")
@validate_call
def query_gswlc_catalog(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    columns: Annotated[str | None, Field(description="Comma-separated column names (all if null)")] = None,
    filters: Annotated[str | None, Field(description='JSON filter dict: {"logMstar":{"min":10},"z":{"max":0.1}}')] = None,
    limit: Annotated[int, Field(ge=1, le=1000000, description="Max rows to return")] = 10000,
    sort_by: Annotated[str | None, Field(description="Column to sort by")] = None,
    sort_ascending: Annotated[bool, Field(description="Sort ascending")] = True,
) -> ArtifactResult:
    """Query the GSWLC-A2 catalog with flexible column selection and filtering.

    Use this to extract galaxy subsets from the 640k-galaxy GSWLC-A2 catalog.
    First call downloads the catalog (~35 MB compressed) if not cached.
    Filters use the same syntax as the reference code: min/max, eq, in, gt/lt/gte/lte.
    """
    import json

    df = _load_catalog()

    if filters:
        filt = json.loads(filters)
        df = _apply_filters(df, filt)

    if sort_by:
        if sort_by not in df.columns:
            raise ValueError(f"Unknown sort column {sort_by!r}")
        df = df.sort_values(by=sort_by, ascending=sort_ascending)

    df = df.head(limit)

    if columns:
        col_list = [c.strip() for c in columns.split(",")]
        missing = [c for c in col_list if c not in df.columns]
        if missing:
            raise ValueError(f"Unknown columns: {missing}. Valid: {list(COLUMN_NAMES)}")
        df = df[col_list]

    out = resolve_outdir(output_dir)
    path = out / "gswlc_query.csv"
    df.to_csv(path, index=False)
    return ok(f"GSWLC-A2: {len(df)} galaxies returned.",
              [path], n_rows=len(df), columns=list(df.columns))


@data_kind("catalog")
@validate_call
def describe_gswlc_catalog() -> ArtifactResult:
    """Describe the GSWLC-A2 catalog: columns, types, and filter syntax.

    Use this first to understand what columns are available and how to
    construct filters for query_gswlc_catalog.
    """
    return ok(
        "GSWLC-A2: GALEX-SDSS-WISE Legacy Catalog (Salim et al. 2016, 2018). "
        "640,659 galaxies with z < 0.3, stellar masses, SFRs, dust from SED fitting.",
        [],
        columns=COLUMN_DESCRIPTIONS,
        filter_syntax={
            "range": '{"logMstar": {"min": 10, "max": 11}}',
            "equality": '{"flag_mgs": {"eq": 1}}',
            "list": '{"plate": {"in": [1234, 5678]}}',
            "comparison": '{"z": {"gt": 0.05, "lt": 0.1}}',
        },
        catalog_size=640659,
        citation="Salim et al., ApJS 227, 2 (2016); Salim, Boquien & Lee, ApJ 859, 11 (2018)",
    )


@data_kind("catalog")
@validate_call
def plot_star_forming_main_sequence(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    input_file: Annotated[str | None, Field(description="Path to CSV from query_gswlc_catalog (null = use full catalog)")] = None,
    max_galaxies: Annotated[int, Field(ge=1000, le=700000)] = 100000,
) -> ArtifactResult:
    """Plot the star-forming main sequence: log SFR vs log M* coloured by density.

    Use after query_gswlc_catalog or directly (will load the full catalog).
    Reproduces the classic SFR-M* relation showing the main sequence, green
    valley, and quenched population.
    """
    import pandas as pd
    plt = mpl()

    if input_file:
        df = pd.read_csv(input_file)
    else:
        df = _load_catalog()

    # Require logMstar and logSFR
    needed = {"logMstar", "logSFR"}
    if not needed.issubset(df.columns):
        raise ValueError(f"Input must have columns: {needed}")

    # Filter out bad SED fits and non-detections
    mask = (df["logMstar"] > 7) & (df["logMstar"] < 13) & (df["logSFR"] > -5) & (df["logSFR"] < 3)
    if "flag_sed" in df.columns:
        mask &= df["flag_sed"] == 0
    df = df[mask]

    if len(df) > max_galaxies:
        df = df.sample(n=max_galaxies, random_state=42)

    fig, ax = plt.subplots(figsize=(8, 6))
    h = ax.hexbin(df["logMstar"].values, df["logSFR"].values,
                  gridsize=80, cmap="viridis", mincnt=1, bins="log")
    fig.colorbar(h, ax=ax, label="log(count)")

    # Main sequence reference line: logSFR = 0.76*logM - 7.64 (Salim+07)
    x_ms = np.linspace(8.5, 11.5, 100)
    y_ms = 0.76 * x_ms - 7.64
    ax.plot(x_ms, y_ms, "r--", lw=2, label="MS (Salim+07)")

    ax.set_xlabel(r"$\log\, M_\star / M_\odot$")
    ax.set_ylabel(r"$\log\, \mathrm{SFR}\, [M_\odot/\mathrm{yr}]$")
    ax.set_title("Star-Forming Main Sequence (GSWLC-A2)")
    ax.set_xlim(8.5, 12)
    ax.set_ylim(-3, 2)
    ax.legend()
    out = resolve_outdir(output_dir)
    path = save_fig(fig, out, "sfms")
    return ok(f"Star-forming main sequence plotted for {len(df)} galaxies.", [path])


@data_kind("catalog")
@validate_call
def crossmatch_sdss_gswlc(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    sdss_file: Annotated[str, Field(description="Path to SDSS CSV (from fetch_sdss_galaxies)")],
    match_radius_arcsec: Annotated[float, Field(ge=0.5, le=30.0)] = 2.0,
) -> ArtifactResult:
    """Cross-match an SDSS galaxy catalog with GSWLC-A2 on sky coordinates.

    Matches within match_radius_arcsec using astropy SkyCoord. Returns a
    merged CSV with columns from both catalogs. Use ObjID matching when both
    have it (faster), else falls back to ra/dec.
    """
    import pandas as pd
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    sdss = pd.read_csv(sdss_file)
    gswlc = _load_catalog()

    # Try ObjID matching first
    if "objID" in sdss.columns and "ObjID" in gswlc.columns:
        merged = pd.merge(sdss, gswlc, left_on="objID", right_on="ObjID",
                          how="inner", suffixes=("_sdss", "_gswlc"))
        if len(merged) > 0:
            out = resolve_outdir(output_dir)
            path = out / "sdss_gswlc_matched.csv"
            merged.to_csv(path, index=False)
            return ok(f"ObjID match: {len(merged)} galaxies.",
                      [path], n_matched=len(merged), method="ObjID",
                      columns=list(merged.columns))

    # Fallback: coordinate matching
    sdss_coords = SkyCoord(ra=sdss["ra"].values, dec=sdss["dec"].values, unit="deg")
    gswlc_coords = SkyCoord(ra=gswlc["ra"].values, dec=gswlc["dec"].values, unit="deg")
    idx, d2d, _ = sdss_coords.match_to_catalog_sky(gswlc_coords)
    mask = d2d < match_radius_arcsec * u.arcsec

    matched_sdss = sdss[mask].reset_index(drop=True)
    matched_gswlc = gswlc.iloc[idx[mask]].reset_index(drop=True)
    merged = pd.concat([matched_sdss.add_suffix("_sdss"),
                        matched_gswlc.add_suffix("_gswlc")], axis=1)

    out = resolve_outdir(output_dir)
    path = out / "sdss_gswlc_matched.csv"
    merged.to_csv(path, index=False)
    return ok(f"Coordinate match: {len(merged)} galaxies within {match_radius_arcsec}\".",
              [path], n_matched=len(merged), method="sky_coord",
              match_radius_arcsec=match_radius_arcsec)
