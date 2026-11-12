"""DES Y6 angular BAO: chi2(alpha) likelihood profiles + released clustering data vectors.

Source: DES Y6 BAO public release (no login)
  https://desdr-server.ncsa.illinois.edu/despublic/y6a2_files/y6_bao/
  (manifest y6_bao_manifest.20250124.json; release page
  https://des.ncsa.illinois.edu/releases/y6a2/Y6bao)

DES Y6 BAO measures alpha = (D_M/r_d) / (D_M/r_d)_Planck at z_eff = 0.851
from photometric red galaxies (6 tomographic bins, 0.6 < z_ph < 1.2) using
three statistics: angular correlation function w(theta) (ACF), angular power
spectrum C_ell (APS) and projected correlation function xi_p(s_perp) (PCF).
The likelihood is NON-Gaussian; it is released as Delta chi2(alpha).

Bundled in data/bao/des_y6/ (small text files, ~0.5 MB):
  likelihood_{fid,acf,aps,pcf}.txt   -- DES release Likelihood/DESY6BAO_likelihood.zip
                                        (fid = ACF+APS+PCF consensus)
  likelihood_DESY6_{DR1tiles,deccut}_noDESI.csv
                                     -- DESI-overlap-removed profiles
                                        (Mena-Fernandez et al., arXiv:2601.14864;
                                        cosmosis-standard-library likelihood/bao/des/y6-no-desi)
Downloaded on demand (cached): Data/Clustering/DESY6BAO_datavectors.zip and
Covariances/DESY6BAO_covariances.zip.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np

from ..common import data_cache_dir, fetch_file

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "bao" / "des_y6"
BASE_URL = "https://desdr-server.ncsa.illinois.edu/despublic/y6a2_files/y6_bao"
DATAVEC_URL = f"{BASE_URL}/Data/Clustering/DESY6BAO_datavectors.zip"
COV_URL = f"{BASE_URL}/Covariances/DESY6BAO_covariances.zip"
DATAVEC_SHA256 = "17541b7dacd8b8fc13123865d9553629ee7c9524db16d94ee1d3bfdab12db977"
COV_SHA256 = "1cb1512be1c953d5a5c2ec6c17993235068d62a4c3805aa7239138e1adfd3f4e"

Z_EFF = 0.851
# Fiducial (Planck 2018) D_M/r_d used to define alpha. The release README rounds
# it to 20.38; the official CosmoSIS module uses D_A = 1627.70 Mpc, r_d = 147.8 Mpc.
DM_RD_FID = 1627.70 * (1.0 + Z_EFF) / 147.8   # = 20.3848

CITATION_Y6 = "DES Collaboration, arXiv:2402.10696 (2024)"
CITATION_NODESI = "Mena-Fernandez et al. (DES), Phys. Rev. D 113, 123521 (2026), arXiv:2601.14864"

# variant -> (file, description, published alpha, published sigma, citation)
VARIANTS: dict[str, dict] = {
    "fid": {"file": "likelihood_fid.txt", "label": "DES Y6 BAO (ACF+APS+PCF consensus)",
            "published_alpha": 0.9571, "published_sigma": 0.0201,
            "published_DM_over_rd": 19.51, "published_DM_over_rd_sigma": 0.41,
            "citation": CITATION_Y6},
    "acf": {"file": "likelihood_acf.txt", "label": "DES Y6 BAO w(theta) only", "citation": CITATION_Y6},
    "aps": {"file": "likelihood_aps.txt", "label": "DES Y6 BAO C_ell only", "citation": CITATION_Y6},
    "pcf": {"file": "likelihood_pcf.txt", "label": "DES Y6 BAO xi_p(s_perp) only", "citation": CITATION_Y6},
    "nodesi_dr1tiles": {"file": "likelihood_DESY6_DR1tiles_noDESI.csv",
                        "label": "DES Y6 BAO no-DESI (DESI DR1 tiles removed; for DESI DR1/DR2 joint fits)",
                        "citation": CITATION_NODESI},
    "nodesi_deccut": {"file": "likelihood_DESY6_deccut_noDESI.csv",
                      "label": "DES Y6 BAO no-DESI (Dec cut; for future DESI DR3/DR4 joint fits)",
                      "citation": CITATION_NODESI},
}

# BAO dataset keys (used by tools.bao) -> variant
DATASET_VARIANT = {"des_y6": "fid", "des_y6_nodesi": "nodesi_dr1tiles"}


def load_profile(variant: str = "fid") -> dict:
    """Return the chi2(alpha) profile (min-subtracted) + derived constraint."""
    if variant not in VARIANTS:
        raise ValueError(f"Unknown DES Y6 BAO variant {variant!r}; choose from {list(VARIANTS)}")
    info = VARIANTS[variant]
    path = DATA_DIR / info["file"]
    if path.suffix == ".csv":
        arr = np.loadtxt(path, delimiter=",", skiprows=1)
    else:
        arr = np.loadtxt(path)
    alpha, chi2 = arr[:, 0], arr[:, 1] - arr[:, 1].min()
    order = np.argsort(alpha)
    alpha, chi2 = alpha[order], chi2[order]
    imin = int(np.argmin(chi2))
    lo, hi = _crossing(alpha, chi2, imin, 1.0)
    return {
        "variant": variant, "label": info["label"], "citation": info["citation"],
        "z": Z_EFF, "observable": "DM_over_rs", "fid": DM_RD_FID,
        "alpha": alpha, "chi2": chi2,
        "alpha_best": float(alpha[imin]), "alpha_lo": lo, "alpha_hi": hi,
        "alpha_sigma": 0.5 * (hi - lo),
        "DM_over_rd": float(alpha[imin] * DM_RD_FID),
        "DM_over_rd_sigma": float(0.5 * (hi - lo) * DM_RD_FID),
        "published": {k: v for k, v in info.items() if k.startswith("published")},
        "file": str(path),
    }


def _crossing(x, y, imin, level):
    """Linear-interpolated points where y crosses ``level`` either side of imin."""
    def side(idx):
        for a, b in zip(idx[:-1], idx[1:], strict=False):
            if y[a] < level <= y[b]:
                return float(x[a] + (level - y[a]) * (x[b] - x[a]) / (y[b] - y[a]))
        return float(x[idx[-1]])
    return side(list(range(imin, -1, -1))), side(list(range(imin, len(x))))


def profile_payload(variant: str) -> dict:
    """JSON-safe profile for the dispatchable fit kernel."""
    p = load_profile(variant)
    return {"z": p["z"], "observable": p["observable"], "fid": p["fid"],
            "alpha": p["alpha"].tolist(), "chi2": p["chi2"].tolist(), "label": p["label"]}


# ------------------------------------------------------------------ clustering data vectors
STATISTICS = {
    "acf": {"x": "theta_rad", "f": "w_theta", "cov": "acf_cov_{c}.txt", "binning": "acf_binning.txt",
            "desc": "angular correlation function w(theta); theta in radians"},
    "aps": {"x": "ell", "f": "C_ell", "cov": "aps_cov_{c}.txt", "binning": "aps_binning.txt",
            "desc": "angular power spectrum C_ell (bandpowers)"},
    "pcf": {"x": "s_perp_Mpc_h", "f": "xi_p", "cov": "xip_cov_{c}.txt", "binning": "pcf_binning.txt",
            "desc": "projected correlation function xi_p(s_perp); s_perp in Mpc/h"},
}
COVARIANCES = {"planck_datalike": "Planck_datalike", "mice_mocklike": "MICE_mocklike"}
N_ZBINS = 6


def _unzipped(url: str, sha256: str, name: str) -> Path:
    zpath = fetch_file(url, subdir="des_y6_bao", sha256=sha256)
    dest = data_cache_dir() / "des_y6_bao" / name
    if not (dest.is_dir() and (any(dest.rglob("*.txt")) or any(dest.rglob("*.dat")))):
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zpath) as zf:
            zf.extractall(dest)
    return dest


def load_clustering(statistic: str = "acf", covariance: str = "planck_datalike") -> dict:
    """Load a DES Y6 BAO clustering data vector (6 z-bins) with matching covariance.

    The released covariances span each statistic's full binning; the released
    data vectors may use a sub-range (PCF: 20 of 34 s_perp bins). The data vector
    is matched to its covariance block by the error column (released errors come
    from the Planck data-like covariance) and the covariance is cut accordingly.
    """
    if statistic not in STATISTICS:
        raise ValueError(f"statistic must be one of {list(STATISTICS)}")
    if covariance not in COVARIANCES:
        raise ValueError(f"covariance must be one of {list(COVARIANCES)}")
    st = STATISTICS[statistic]
    dv_dir = _unzipped(DATAVEC_URL, DATAVEC_SHA256, "datavectors")
    cov_dir = _unzipped(COV_URL, COV_SHA256, "covariances")
    cov_root = next(p.parent for p in cov_dir.rglob(st["binning"]))
    cov_full = np.loadtxt(cov_root / st["cov"].format(c=COVARIANCES[covariance]))
    cov_planck = np.loadtxt(cov_root / st["cov"].format(c="Planck_datalike"))
    binning = np.loadtxt(cov_root / st["binning"])
    nper = cov_full.shape[0] // N_ZBINS
    sig_full = np.sqrt(np.diag(cov_planck))

    zbin, xs, fs, es, idx = [], [], [], [], []
    for b in range(1, N_ZBINS + 1):
        f = next(dv_dir.rglob(f"{statistic}_data_bin{b}.dat"))
        a = np.loadtxt(f)
        block = sig_full[(b - 1) * nper: b * nper]
        for row in a:
            k = int(np.argmin(np.abs(block / row[2] - 1.0)))
            if abs(block[k] / row[2] - 1.0) > 1e-4:
                raise RuntimeError(f"Could not match {f.name} x={row[0]} to a covariance bin")
            idx.append((b - 1) * nper + k)
            zbin.append(b)
            xs.append(row[0])
            fs.append(row[1])
            es.append(row[2])
    idx = np.array(idx)
    cov = cov_full[np.ix_(idx, idx)]
    return {"statistic": statistic, "covariance_kind": covariance, "description": st["desc"],
            "x_name": st["x"], "f_name": st["f"],
            "zbin": np.array(zbin), "x": np.array(xs), "value": np.array(fs), "error": np.array(es),
            "cov": cov, "cov_index": idx, "binning": binning, "n_per_bin_full": nper}
