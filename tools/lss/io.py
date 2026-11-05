"""I/O helpers for Large-Scale Structure clustering catalogs.

Handles partial / subset downloads of BOSS DR12, eBOSS DR16 and DESI DR1
clustering catalogs (FITS, possibly gzipped).  Downloads are always subsets
(row-limited, z-range filtered, random fraction) so the user never has to
pull hundreds of GB.
"""
from __future__ import annotations

import gzip
import re
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits

from ..common import data_cache_dir, fetch_file

# ---------- catalog registry --------------------------------------------------

_TF = {"L": "i1", "B": "u1", "I": ">i2", "J": ">i4", "K": ">i8", "E": ">f4", "D": ">f8"}

_SAS16 = "https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/"
_SAS12 = "https://data.sdss.org/sas/dr12/boss/lss/"
_DESI1 = "https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/"

_BOSS12_PAPERS = "https://data.sdss.org/sas/dr12/boss/papers/clustering/"
_EBOSS_DVC = "https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/dataveccov/lrg_elg_qso/"

# key -> (data_url, random_url, nz_url_or_None, is_gzipped, weight_cols, z_col, survey)
LSS_CATALOGS: dict[str, dict] = {
    "eboss_dr16_lrg_ngc": {
        "data_url": _SAS16 + "eBOSS_LRG_clustering_data-NGC-vDR16.fits",
        "random_url": _SAS16 + "eBOSS_LRG_clustering_random-NGC-vDR16.fits",
        "nz_url": None,
        "gzipped": False,
        "weight_cols": ["WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.6, 1.0), "survey": "eBOSS DR16 LRG NGC",
        "citation": "Alam et al., Phys. Rev. D 103, 083533 (2021)",
    },
    "eboss_dr16_qso_ngc": {
        "data_url": _SAS16 + "eBOSS_QSO_clustering_data-NGC-vDR16.fits",
        "random_url": _SAS16 + "eBOSS_QSO_clustering_random-NGC-vDR16.fits",
        "nz_url": None,
        "gzipped": False,
        "weight_cols": ["WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.8, 2.2), "survey": "eBOSS DR16 QSO NGC",
        "citation": "Alam et al., Phys. Rev. D 103, 083533 (2021)",
    },
    "eboss_dr16_elg_ngc": {
        "data_url": _SAS16 + "eBOSS_ELG_clustering_data-NGC-vDR16.fits",
        "random_url": _SAS16 + "eBOSS_ELG_clustering_random-NGC-vDR16.fits",
        "nz_url": None,
        "gzipped": False,
        "weight_cols": ["WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.6, 1.1), "survey": "eBOSS DR16 ELG NGC",
        "citation": "Alam et al., Phys. Rev. D 103, 083533 (2021)",
    },
    "boss_dr12_cmass_north": {
        "data_url": _SAS12 + "galaxy_DR12v5_CMASS_North.fits.gz",
        "random_url": _SAS12 + "random0_DR12v5_CMASS_North.fits.gz",
        "nz_url": None,
        "gzipped": True,
        "weight_cols": ["WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.43, 0.7), "survey": "BOSS DR12 CMASS North",
        "citation": "Alam et al., MNRAS 470, 2617 (2017)",
    },
    "boss_dr12_lowz_north": {
        "data_url": _SAS12 + "galaxy_DR12v5_LOWZ_North.fits.gz",
        "random_url": _SAS12 + "random0_DR12v5_LOWZ_North.fits.gz",
        "nz_url": None,
        "gzipped": True,
        "weight_cols": ["WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.15, 0.43), "survey": "BOSS DR12 LOWZ North",
        "citation": "Alam et al., MNRAS 470, 2617 (2017)",
    },
    "desi_dr1_lrg_ngc": {
        "data_url": _DESI1 + "LRG_NGC_clustering.dat.fits",
        "random_url": _DESI1 + "LRG_NGC_0_clustering.ran.fits",
        "nz_url": _DESI1 + "LRG_NGC_nz.txt",
        "gzipped": False,
        "weight_cols": ["WEIGHT", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.4, 1.1), "survey": "DESI DR1 LRG NGC",
        "citation": "DESI Collaboration, arXiv:2306.06308 (2023)",
    },
    "desi_dr1_elg_ngc": {
        "data_url": _DESI1 + "ELG_LOPnotqso_NGC_clustering.dat.fits",
        "random_url": _DESI1 + "ELG_LOPnotqso_NGC_0_clustering.ran.fits",
        "nz_url": _DESI1 + "ELG_LOPnotqso_NGC_nz.txt",
        "gzipped": False,
        "weight_cols": ["WEIGHT", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.8, 1.6), "survey": "DESI DR1 ELG NGC",
        "citation": "DESI Collaboration, arXiv:2306.06308 (2023)",
    },
    "desi_dr1_qso_ngc": {
        "data_url": _DESI1 + "QSO_NGC_clustering.dat.fits",
        "random_url": _DESI1 + "QSO_NGC_0_clustering.ran.fits",
        "nz_url": _DESI1 + "QSO_NGC_nz.txt",
        "gzipped": False,
        "weight_cols": ["WEIGHT", "WEIGHT_FKP"],
        "z_col": "Z", "z_range": (0.8, 2.1), "survey": "DESI DR1 QSO NGC",
        "citation": "DESI Collaboration, arXiv:2306.06308 (2023)",
    },
}

# ---------- released multipole files (summary statistics) ---------------------

# BOSS DR12 Beutler+17 P(k) multipoles (fullshape prerecon)
PUBLISHED_MULTIPOLES: dict[str, dict] = {
    "boss_dr12_pk_z1_ngc": {
        "type": "pk", "survey": "BOSS DR12 NGC z1 (zeff=0.38)",
        "data_urls": {
            "P0": _BOSS12_PAPERS + "Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz",
            "P2": _BOSS12_PAPERS + "Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz",
            "P4": _BOSS12_PAPERS + "Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz",
        },
        "tar_members": {
            "P0": "public_material_RSD/Beutleretal_pk_monopole_DR12_NGC_z1_prerecon_120.dat",
            "P2": "public_material_RSD/Beutleretal_pk_quadrupole_DR12_NGC_z1_prerecon_120.dat",
            "P4": "public_material_RSD/Beutleretal_pk_hexadecapole_DR12_NGC_z1_prerecon_120.dat",
            "cov": "public_material_RSD/Beutleretal_cov_patchy_z1_NGC_1_15_1_15_1_10_2045_60.dat",
        },
        "citation": "Beutler et al., MNRAS 466, 2242 (2017)",
    },
    "boss_dr12_pk_z3_ngc": {
        "type": "pk", "survey": "BOSS DR12 NGC z3 (zeff=0.61)",
        "data_urls": {
            "P0": _BOSS12_PAPERS + "Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz",
        },
        "tar_members": {
            "P0": "public_material_RSD/Beutleretal_pk_monopole_DR12_NGC_z3_prerecon_120.dat",
            "P2": "public_material_RSD/Beutleretal_pk_quadrupole_DR12_NGC_z3_prerecon_120.dat",
            "P4": "public_material_RSD/Beutleretal_pk_hexadecapole_DR12_NGC_z3_prerecon_120.dat",
            "cov": "public_material_RSD/Beutleretal_cov_patchy_z3_NGC_1_15_1_15_1_10_2045_60.dat",
        },
        "citation": "Beutler et al., MNRAS 466, 2242 (2017)",
    },
    "eboss_dr16_lrg_xi": {
        "type": "xi", "survey": "eBOSS DR16 LRG (NGC+SGC, 0.6<z<1.0)",
        "data_urls": {
            "xi": _EBOSS_DVC + "LRG_xi/Data_LRGxi_NGCSGC_0.6z1.0_prerecon.txt",
            "cov": _EBOSS_DVC + "LRG_xi/Covariance_LRGxi_NGCSGC_0.6z1.0_prerecon.txt",
        },
        "citation": "Bautista et al., MNRAS 500, 736 (2021); arXiv:2007.08993",
    },
    "eboss_dr16_qso_xi": {
        "type": "xi", "survey": "eBOSS DR16 QSO (NGC+SGC, 0.8<z<2.2)",
        "data_urls": {
            "xi": _EBOSS_DVC + "QSO_xi/Data_QSOxi_NGCSGC_0.8z2.2_prerecon.txt",
            "cov": _EBOSS_DVC + "QSO_xi/Covariance_QSOxi_NGCSGC_0.8z2.2_prerecon.txt",
        },
        "citation": "Hou et al., MNRAS 500, 1201 (2021); arXiv:2007.08998",
    },
}


# ---------- FITS partial-read helpers -----------------------------------------

def _fits_dtype(header):
    """Build a numpy structured dtype from a FITS binary table header."""
    fields = []
    for i in range(1, header["TFIELDS"] + 1):
        form = header[f"TFORM{i}"].strip()
        m = re.match(r"(\d*)([A-Z])", form)
        if not m:
            raise ValueError(f"Cannot parse TFORM={form}")
        rep = int(m.group(1) or 1)
        code = m.group(2)
        name = header[f"TTYPE{i}"]
        if code == "A":
            fields.append((name, f"S{rep}"))
        elif rep == 1:
            fields.append((name, _TF[code]))
        else:
            fields.append((name, _TF[code], (rep,)))
    dt = np.dtype(fields)
    assert dt.itemsize == header["NAXIS1"], (dt.itemsize, header["NAXIS1"])
    return dt


def _range_bytes(url: str, start: int, end: int) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}",
                                               "User-Agent": "cosmo-mcp-server/0.1"})
    return urllib.request.urlopen(req, timeout=120).read()


def _parse_header_ends(buf: bytes):
    txt = buf.decode("latin1")
    cards = [txt[i:i + 80] for i in range(0, len(txt), 80)]
    return [i for i, c in enumerate(cards) if c.rstrip() == "END"]


def _fits_header_range(url: str):
    """Read a FITS header via HTTP Range, return (header, data_offset)."""
    buf = _range_bytes(url, 0, 57599)
    ends = _parse_header_ends(buf)
    if len(ends) < 2:
        buf += _range_bytes(url, len(buf), len(buf) + 57600 - 1)
        ends = _parse_header_ends(buf)
    p_end = ((ends[0] + 1) * 80 + 2879) // 2880 * 2880
    d0 = ((ends[1] + 1) * 80 + 2879) // 2880 * 2880
    header = fits.Header.fromstring(buf[p_end:d0].decode("latin1"))
    return header, d0


def _fits_header_gzip(url: str):
    """Read a gzipped FITS header by streaming, return (header, gzip_stream)."""
    resp = urllib.request.urlopen(urllib.request.Request(
        url, headers={"User-Agent": "cosmo-mcp-server/0.1"}), timeout=120)
    g = gzip.GzipFile(fileobj=resp)
    buf = b""
    while True:
        buf += g.read(2880)
        ends = _parse_header_ends(buf)
        if len(ends) >= 2:
            break
    p_end = ((ends[0] + 1) * 80 + 2879) // 2880 * 2880
    header = fits.Header.fromstring(buf[p_end:].decode("latin1"))
    return header, g


def read_fits_rows_range(url: str, row_start: int, n_rows: int,
                         header, data_offset: int, dt) -> np.ndarray:
    """Read n_rows starting at row_start from an uncompressed remote FITS."""
    start = data_offset + row_start * dt.itemsize
    end = start + n_rows * dt.itemsize - 1
    raw = _range_bytes(url, start, end)
    return np.frombuffer(raw, dtype=dt)


def read_fits_rows_gzip(url: str, n_rows: int, header, dt) -> np.ndarray:
    """Read the first n_rows from a gzipped remote FITS (streaming)."""
    hdr, g = _fits_header_gzip(url)
    raw = g.read(n_rows * dt.itemsize)
    return np.frombuffer(raw, dtype=dt)


def download_lss_subset(catalog_key: str,
                        kind: str = "data",
                        max_rows: int = 50000,
                        z_min: float | None = None,
                        z_max: float | None = None,
                        random_fraction: float = 1.0,
                        output_dir: Path | None = None,
                        ) -> tuple[Path, dict]:
    """Download a subset of an LSS clustering catalog.

    Returns (csv_path, metadata_dict).
    """
    import pandas as pd

    cat = LSS_CATALOGS[catalog_key]
    url = cat["data_url"] if kind == "data" else cat["random_url"]
    z_col = cat["z_col"]
    is_gz = cat["gzipped"]
    wt_cols = cat["weight_cols"]

    # Determine columns to keep
    keep_cols = ["RA", "DEC", z_col] + wt_cols
    if "NZ" not in keep_cols:
        keep_cols.append("NZ")

    if is_gz:
        # For gzipped: stream rows from the start
        header, g = _fits_header_gzip(url)
        dt = _fits_dtype(header)
        total_rows = header["NAXIS2"]
        # Read up to max_rows * 1/random_fraction to have enough after sampling
        read_n = min(total_rows, int(max_rows / max(random_fraction, 0.01) * 1.5))
        raw = g.read(read_n * dt.itemsize)
        arr = np.frombuffer(raw, dtype=dt)
    else:
        # For uncompressed: use HTTP Range for strided reads
        header, d0 = _fits_header_range(url)
        dt = _fits_dtype(header)
        total_rows = header["NAXIS2"]
        # Catalog files are ordered (by sky position / chunk), so a contiguous head block
        # covers a different footprint for data and randoms and biases xi(s) (constant
        # offset at large s). Read evenly spaced blocks spanning the whole file instead.
        read_n = min(total_rows, int(max_rows / max(random_fraction, 0.01) * 1.5))
        if read_n >= total_rows:
            arr = read_fits_rows_range(url, 0, total_rows, header, d0, dt)
        else:
            n_blocks = int(min(128, max(1, read_n // 500)))
            per = read_n // n_blocks
            starts = np.linspace(0, total_rows - per, n_blocks).astype(int)
            arr = np.concatenate([read_fits_rows_range(url, int(s0), per, header, d0, dt) for s0 in starts])

    # Build DataFrame with only the columns we want
    available = set(dt.names)
    cols = [c for c in keep_cols if c in available]
    data = {c: np.array(arr[c], dtype=float) for c in cols}
    df = pd.DataFrame(data)

    # z-range filter
    zlo = z_min if z_min is not None else (cat["z_range"][0] if "z_range" in cat else -np.inf)
    zhi = z_max if z_max is not None else (cat["z_range"][1] if "z_range" in cat else np.inf)
    df = df[(df[z_col] >= zlo) & (df[z_col] <= zhi)]

    # Random subsample
    if random_fraction < 1.0 and len(df) > 0:
        n_keep = max(1, int(len(df) * random_fraction))
        df = df.sample(n=min(n_keep, len(df)), random_state=42)

    # Row limit (random subsample, not head, to keep the footprint uniform)
    if len(df) > max_rows:
        df = df.sample(n=max_rows, random_state=0)

    df = df.reset_index(drop=True)

    # Compute combined weight column (product of all weight cols present)
    wt_present = [c for c in wt_cols if c in df.columns]
    if wt_present:
        w = np.ones(len(df))
        for c in wt_present:
            w *= df[c].values
        df["WEIGHT_TOTAL"] = w

    outdir = output_dir or data_cache_dir() / "lss"
    outdir.mkdir(parents=True, exist_ok=True)
    tag = f"{catalog_key}_{kind}"
    csv_path = outdir / f"{tag}.csv"
    df.to_csv(csv_path, index=False)

    meta = {
        "catalog": catalog_key, "kind": kind, "survey": cat["survey"],
        "url": url, "total_rows_remote": total_rows,
        "rows_read": int(len(arr)), "rows_after_filter": len(df),
        "z_range": [float(zlo), float(zhi)],
        "columns": list(df.columns),
        "citation": cat["citation"],
    }
    return csv_path, meta


# ---------- published multipole fetchers --------------------------------------

def fetch_published_xi(key: str, output_dir: Path) -> tuple[Path, dict]:
    """Fetch released xi(s) multipoles (eBOSS DR16 format)."""
    entry = PUBLISHED_MULTIPOLES[key]
    data_url = entry["data_urls"]["xi"]
    cov_url = entry["data_urls"]["cov"]

    output_dir.mkdir(parents=True, exist_ok=True)
    data_path = fetch_file(data_url, subdir="lss_multipoles")
    cov_path = fetch_file(cov_url, subdir="lss_multipoles")

    # Parse data file: columns are s xi0 xi0_err xi2 xi2_err [xi4 xi4_err]
    lines = [ln for ln in data_path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    rows = [list(map(float, ln.split())) for ln in lines]
    arr = np.array(rows)
    cols = ["s", "xi0", "xi0_err", "xi2", "xi2_err"]
    if arr.shape[1] >= 7:
        cols += ["xi4", "xi4_err"]

    import pandas as pd
    df = pd.DataFrame(arr[:, :len(cols)], columns=cols)
    out_csv = output_dir / f"published_{key}_data.csv"
    df.to_csv(out_csv, index=False)

    # Parse covariance
    cov_lines = [ln for ln in cov_path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    cov_rows = [list(map(float, ln.split())) for ln in cov_lines]
    cov_arr = np.array(cov_rows)
    cov_csv = output_dir / f"published_{key}_cov.csv"
    np.savetxt(cov_csv, cov_arr, delimiter=",")

    return out_csv, {
        "key": key, "type": "xi", "survey": entry["survey"],
        "n_bins": len(df), "s_range": [float(df["s"].min()), float(df["s"].max())],
        "columns": cols, "citation": entry["citation"],
        "cov_shape": list(cov_arr.shape),
    }


def fetch_published_pk_beutler(key: str, output_dir: Path) -> tuple[Path, dict]:
    """Fetch Beutler+17 P(k) multipoles from the BOSS DR12 tarball."""
    import tarfile

    entry = PUBLISHED_MULTIPOLES[key]
    tar_url = list(entry["data_urls"].values())[0]
    tar_path = fetch_file(tar_url, subdir="lss_multipoles")

    output_dir.mkdir(parents=True, exist_ok=True)
    members = entry["tar_members"]

    results = {}
    with tarfile.open(tar_path, "r:gz") as tf:
        for label, member_name in members.items():
            if label == "cov":
                continue
            f = tf.extractfile(member_name)
            if f is None:
                continue
            lines = [ln.decode() for ln in f.readlines()
                     if not ln.decode().strip().startswith("#") and not ln.decode().strip().startswith("###")
                     and ln.decode().strip()]
            data_lines = [ln for ln in lines if len(ln.split()) >= 4]
            rows = [list(map(float, ln.split())) for ln in data_lines]
            arr = np.array(rows)
            results[label] = arr  # cols: k_central k_mean P_ell(k) sig_ii

    # Combine into one CSV
    import pandas as pd
    frames = []
    for label in ["P0", "P2", "P4"]:
        if label in results:
            a = results[label]
            df = pd.DataFrame({"k": a[:, 1], label: a[:, 2], f"{label}_err": a[:, 3]})
            frames.append(df)

    if not frames:
        raise RuntimeError(f"No multipole data found for {key}")

    merged = frames[0]
    for f in frames[1:]:
        merged = pd.merge(merged, f, on="k", how="outer")
    merged = merged.sort_values("k").reset_index(drop=True)
    out_csv = output_dir / f"published_{key}_data.csv"
    merged.to_csv(out_csv, index=False)

    # Covariance
    if "cov" in members:
        cov_member = members["cov"]
        with tarfile.open(tar_path, "r:gz") as tf:
            f = tf.extractfile(cov_member)
            cov_lines = [ln.decode() for ln in f.readlines()
                         if not ln.decode().strip().startswith("#") and not ln.decode().strip().startswith("###")
                         and ln.decode().strip()]
            cov_data = [list(map(float, ln.split())) for ln in cov_lines]
            cov_arr = np.array(cov_data)
        cov_csv = output_dir / f"published_{key}_cov.csv"
        np.savetxt(cov_csv, cov_arr, delimiter=",")
    else:
        cov_csv = None

    return out_csv, {
        "key": key, "type": "pk", "survey": entry["survey"],
        "n_k_bins": len(merged),
        "k_range": [float(merged["k"].min()), float(merged["k"].max())],
        "columns": list(merged.columns), "citation": entry["citation"],
    }


def fetch_desi_nz(catalog_key: str, output_dir: Path) -> tuple[Path | None, dict]:
    """Fetch the n(z) text file for a DESI DR1 catalog."""
    cat = LSS_CATALOGS.get(catalog_key)
    if not cat or not cat.get("nz_url"):
        return None, {"error": f"No n(z) URL for {catalog_key}"}
    nz_path = fetch_file(cat["nz_url"], subdir="lss_nz")
    # Parse: zmid zlow zhigh n(z) Nbin Vol_bin
    lines = [ln for ln in nz_path.read_text().splitlines()
             if ln.strip() and not ln.startswith("#")]
    rows = [list(map(float, ln.split())) for ln in lines]
    arr = np.array(rows)
    import pandas as pd
    cols = ["zmid", "zlow", "zhigh", "nz", "Nbin", "Vol_bin"]
    df = pd.DataFrame(arr[:, :len(cols)], columns=cols[:arr.shape[1]])
    out_csv = output_dir / f"nz_{catalog_key}.csv"
    df.to_csv(out_csv, index=False)
    return out_csv, {
        "catalog": catalog_key, "survey": cat["survey"],
        "n_bins": len(df), "z_range": [float(df["zmid"].min()), float(df["zmid"].max())],
        "citation": cat["citation"],
    }
