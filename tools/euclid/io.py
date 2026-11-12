"""Euclid DR1 I/O: ADQL builders + TAP execution (astroquery.esa.euclid, HTTP fallback)
and the ECLIPSE strong-lens catalogue (Zenodo).

Euclid DR1 ('DR1-Foundation', public 2026-11-12) is served by the ESA Euclid
Science Archive TAP service (anonymous queries allowed):
  https://eas.esac.esa.int/tap-server/tap   (astroquery EuclidClass(environment="PDR"))
Tables discovered from tap_schema (schema 'dr1'):
  dr1.mer_catalogue_wide_survey                    MER photometric object catalogue (Wide Survey)
  dr1.phz_photo_z_wide                             photometric redshifts (phz_median, phz_70_int, ...)
  dr1.phz_classification_wide                      galaxy / star / QSO probabilities
  dr1.spectro_z_spe_galaxy_candidates_wide_survey  NISP spectroscopic z candidates (spe_rank 0..4)
  dr1.spectro_z_spe_classification_wide_survey     spectroscopic class (galaxy/qso/star)
  dr1.spectro_z_spe_quality_wide_survey            spectroscopic quality / error flags
(deep-field equivalents: *_deep_survey / *_deep). All join on object_id.
"""
from __future__ import annotations

import io
import json
import urllib.error
import urllib.parse
import urllib.request

TAP_URL = "https://eas.esac.esa.int/tap-server/tap"
ESA_ENVIRONMENT = "PDR"   # public data release environment of astroquery.esa.euclid

TABLES = {
    "mer": "dr1.mer_catalogue_wide_survey",
    "phz": "dr1.phz_photo_z_wide",
    "phz_class": "dr1.phz_classification_wide",
    "spe_z": "dr1.spectro_z_spe_galaxy_candidates_wide_survey",
    "spe_class": "dr1.spectro_z_spe_classification_wide_survey",
    "spe_quality": "dr1.spectro_z_spe_quality_wide_survey",
}

MAX_RADIUS_DEG = 3.0
MAX_AREA_DEG2 = 30.0

# MER fluxes are in micro-Jansky: AB mag = -2.5 log10(f/uJy) + 23.9.
# Magnitude cuts are applied as FLUX cuts and mag_vis is computed client-side: the
# archive's distributed SQL engine evaluates LOG() on negative fluxes even behind a
# "flux > 0" guard and aborts the whole query.
AB_ZP_UJY = 23.9


def vis_flux_limit(mag_max: float) -> float:
    return 10 ** ((AB_ZP_UJY - mag_max) / 2.5)


def region_clause(ra: float | None, dec: float | None, radius_deg: float | None,
                  bounds: list[float] | None) -> tuple[str, float]:
    """ADQL WHERE fragment for a cone (ra, dec, radius_deg) or box bounds
    [ra_min, ra_max, dec_min, dec_max]; returns (clause, area_deg2)."""
    import math
    if bounds is not None:
        if len(bounds) != 4:
            raise ValueError("bounds must be [ra_min, ra_max, dec_min, dec_max] in degrees")
        ra0, ra1, d0, d1 = (float(b) for b in bounds)
        if not (-90 <= d0 < d1 <= 90):
            raise ValueError("need -90 <= dec_min < dec_max <= 90")
        if not (0 <= ra0 <= 360 and 0 <= ra1 <= 360) or ra0 == ra1:
            raise ValueError("ra_min / ra_max must lie in [0, 360] and differ")
        width = (ra1 - ra0) % 360.0
        area = width * (math.sin(math.radians(d1)) - math.sin(math.radians(d0))) * 180.0 / math.pi
        if ra0 < ra1:
            ra_part = f"m.right_ascension BETWEEN {ra0} AND {ra1}"
        else:   # wraps through RA = 0
            ra_part = f"(m.right_ascension >= {ra0} OR m.right_ascension <= {ra1})"
        return f"{ra_part} AND m.declination BETWEEN {d0} AND {d1}", area
    if ra is None or dec is None or radius_deg is None:
        raise ValueError("give either (ra, dec, radius_deg) for a cone or bounds=[ra_min, ra_max, dec_min, dec_max]")
    if not (0 < radius_deg <= MAX_RADIUS_DEG):
        raise ValueError(f"radius_deg must be in (0, {MAX_RADIUS_DEG}]")
    if not (-90 <= dec <= 90):
        raise ValueError("dec must be in [-90, 90]")
    area = 2 * math.pi * (1 - math.cos(math.radians(radius_deg))) * (180 / math.pi) ** 2
    return (f"1=CONTAINS(POINT('ICRS', m.right_ascension, m.declination), "
            f"CIRCLE('ICRS', {float(ra) % 360.0}, {float(dec)}, {float(radius_deg)}))"), area


def build_query(sample: str, region: str, *, max_rows: int, sample_mod: int = 1, sample_seed: int = 0,
                z_min: float | None = None, z_max: float | None = None, mag_vis_max: float | None = 24.5,
                min_spec_z_prob: float = 0.99, require_galaxy_class: bool = True,
                drop_spurious: bool = True, require_phz_flags_zero: bool = False) -> str:
    """ADQL for a photometric (MER x PHZ) or spectroscopic (MER x SPE) galaxy sample."""
    where = [region]
    if drop_spurious:
        where.append("m.spurious_flag = 0")
    if mag_vis_max is not None:
        where.append(f"m.flux_vis_sersic > {vis_flux_limit(float(mag_vis_max)):.6g}")
    if sample_mod > 1:
        # deterministic pseudo-random subsample, evaluated server-side
        where.append(f"MOD(ABS(m.object_id), {int(sample_mod)}) = {int(sample_seed) % int(sample_mod)}")
    cols = ["m.object_id", "m.right_ascension AS ra", "m.declination AS dec",
            "m.flux_vis_sersic", "m.point_like_prob", "m.det_quality_flag"]
    if sample == "photometric":
        zc = "p.phz_median"
        cols += ["p.phz_median", "p.phz_mode_1", "p.phz_70_int", "p.phz_90_int", "p.phz_flags",
                 "p.phz_classification", "p.tom_bin_id"]
        frm = f"{TABLES['mer']} AS m JOIN {TABLES['phz']} AS p ON m.object_id = p.object_id"
        if require_galaxy_class:
            where.append("p.phz_classification = 2")
        if require_phz_flags_zero:
            where.append("p.phz_flags = 0")
        where.append("p.phz_median IS NOT NULL")
    elif sample == "spectroscopic":
        zc = "s.spe_z"
        cols += ["s.spe_z", "s.spe_z_err", "s.spe_z_prob", "s.spe_z_rel", "c.spe_class", "c.spe_gal_prob"]
        frm = (f"{TABLES['mer']} AS m JOIN {TABLES['spe_z']} AS s ON m.object_id = s.object_id "
               f"JOIN {TABLES['spe_class']} AS c ON m.object_id = c.object_id")
        # spe_rank = 0 is the highest-probability redshift solution in DR1
        where.append("s.spe_rank = 0")
        where.append(f"s.spe_z_prob >= {float(min_spec_z_prob)}")
        if require_galaxy_class:
            where.append("c.spe_class = 'galaxy'")
    else:
        raise ValueError("sample must be 'photometric' or 'spectroscopic'")
    if z_min is not None:
        where.append(f"{zc} >= {float(z_min)}")
    if z_max is not None:
        where.append(f"{zc} <= {float(z_max)}")
    return f"SELECT TOP {int(max_rows)} " + ", ".join(cols) + f" FROM {frm} WHERE " + " AND ".join(where)


def run_tap(adql: str, timeout: float = 300.0):
    """Run ADQL on the ESA Euclid archive -> pandas DataFrame.

    Uses astroquery.esa.euclid (async job); on failure falls back to a plain
    HTTP synchronous TAP request. Returns (df, backend).
    """
    errors = []
    try:
        from astroquery.esa.euclid import EuclidClass
        euclid = EuclidClass(environment=ESA_ENVIRONMENT, show_server_messages=False)
        job = euclid.launch_job_async(adql, verbose=False)
        tab = job.get_results()
        return _postprocess(tab.to_pandas()), "astroquery.esa.euclid (PDR)"
    except Exception as exc:   # noqa: BLE001 — fall back to raw TAP
        errors.append(f"astroquery: {type(exc).__name__}: {exc}")
    try:
        import pandas as pd
        data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                       "QUERY": adql}).encode()
        req = urllib.request.Request(f"{TAP_URL}/sync", data=data,
                                     headers={"User-Agent": "cosmo-mcp-server/0.1"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode()
        except urllib.error.HTTPError as http_exc:
            body = http_exc.read().decode(errors="replace")
            raise RuntimeError(f"HTTP {http_exc.code}: {_votable_error(body)}") from http_exc
        if text.lstrip().startswith("<"):
            raise RuntimeError(_votable_error(text))
        return _postprocess(pd.read_csv(io.StringIO(text))), "TAP sync (HTTP)"
    except Exception as exc:   # noqa: BLE001
        errors.append(f"HTTP TAP: {type(exc).__name__}: {exc}")
    raise RuntimeError("Euclid archive query failed (" + " | ".join(errors) + f"). Service: {TAP_URL}. "
                       "Try a smaller region / max_rows, or retry later.")


def _postprocess(df):
    """Add client-side AB magnitudes (NaN for non-positive flux)."""
    import numpy as np
    if "flux_vis_sersic" in df:
        f = df["flux_vis_sersic"].to_numpy(dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            df["mag_vis"] = np.where(f > 0, -2.5 * np.log10(np.where(f > 0, f, 1.0)) + AB_ZP_UJY, np.nan)
    return df


def _votable_error(text: str) -> str:
    import re
    because = re.search(r"Because:\s*(.*?)\s*</INFO>", text, re.S)
    if because:
        return because.group(1)[:500].strip()
    m = re.search(r'name="QUERY_STATUS" value="ERROR">\s*(.*?)\s*</INFO>', text, re.S)
    return (m.group(1)[:500] if m else text[:500]).strip()


# ------------------------------------------------------------------ ECLIPSE strong lenses
ZENODO_RECORD = "19654904"
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
ECLIPSE_FILE = "ECLIPSE_object_catalogue.zip"
ECLIPSE_CSV = "ECLIPSE_object_catalogue.csv"
ECLIPSE_CITATION = ("Lines, Holloway, Nightingale, Nersesian et al. (Euclid Collaboration), Euclid DR1 Strong Lens "
                    "Catalogue (ECLIPSE), Zenodo, doi:10.5281/zenodo.19654904 (2026)")


def eclipse_file_info() -> dict:
    """List the Zenodo record and return metadata for the object-catalogue file."""
    req = urllib.request.Request(ZENODO_API, headers={"User-Agent": "cosmo-mcp-server/0.1",
                                                      "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            rec = json.load(resp)
    except Exception as exc:
        raise RuntimeError(f"Zenodo API request failed for {ZENODO_API}: {exc}") from exc
    for f in rec.get("files", []):
        if f["key"] == ECLIPSE_FILE:
            return {"key": f["key"], "size": f["size"], "checksum": f.get("checksum"),
                    "url": f["links"]["self"], "record_title": rec["metadata"]["title"],
                    "doi": rec.get("doi"), "all_files": [g["key"] for g in rec["files"]]}
    raise RuntimeError(f"{ECLIPSE_FILE} not found in Zenodo record {ZENODO_RECORD}")


def load_eclipse():
    """Download (once, cached) + read the ECLIPSE object catalogue -> (DataFrame, path, file_info)."""
    import hashlib
    import zipfile

    import pandas as pd

    from ..common import data_cache_dir, fetch_file
    sub = data_cache_dir() / "euclid_dr1" / "eclipse"
    csv = sub / ECLIPSE_CSV
    info = {"key": ECLIPSE_FILE, "url": f"{ZENODO_API}/files/{ECLIPSE_FILE}/content",
            "doi": f"10.5281/zenodo.{ZENODO_RECORD}"}
    if not csv.is_file():
        try:
            info = eclipse_file_info()
        except RuntimeError:
            pass   # fall back to the canonical content URL
        zpath = fetch_file(info["url"], ECLIPSE_FILE, subdir="euclid_dr1/eclipse", timeout=300)
        md5 = (info.get("checksum") or "").removeprefix("md5:")
        if md5 and hashlib.md5(zpath.read_bytes()).hexdigest() != md5:
            zpath.unlink(missing_ok=True)
            raise RuntimeError(f"Checksum mismatch for {ECLIPSE_FILE}; deleted, please retry")
        with zipfile.ZipFile(zpath) as zf:
            member = next(n for n in zf.namelist() if n.endswith(ECLIPSE_CSV))
            csv.write_bytes(zf.read(member))
    return pd.read_csv(csv), csv, info
