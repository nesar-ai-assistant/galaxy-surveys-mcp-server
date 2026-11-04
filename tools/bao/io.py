"""BAO data I/O: load bundled CobayaSampler/bao_data tables."""
from __future__ import annotations

import numpy as np
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "bao"

# ---- dataset definitions ----

# DESI DR2 combined (7 tracers → 13 data points)
DESI_DR2 = {
    "label": "DESI DR2",
    "mean_file": "desi_dr2_ALL_GCcomb_mean.txt",
    "cov_file": "desi_dr2_ALL_GCcomb_cov.txt",
}

# DESI DR1 combined
DESI_DR1 = {
    "label": "DESI DR1",
    "mean_file": "desi_2024_gaussian_bao_ALL_GCcomb_mean.txt",
    "cov_file": "desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
}

# eBOSS DR16 BAO+FS LRG (3 z-bins × 3 observables = 9 points)
EBOSS_DR16_LRG = {
    "label": "eBOSS DR16 LRG (BAO+FS)",
    "mean_file": "sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8.dat",
    "cov_file": "sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8_covtot.txt",
}

EBOSS_DR16_QSO = {
    "label": "eBOSS DR16 QSO (BAO+FS)",
    "mean_file": "sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8.dat",
    "cov_file": "sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8_covtot.txt",
}

BOSS_DR12 = {
    "label": "BOSS DR12 consensus (BAO+FS)",
    "mean_file": "sdss_DR12Consensus_final.dat",
    "cov_file": "final_consensus_covtot_dM_Hz_fsig.txt",
    "rs_fid": 147.78,  # Mpc — stored DM in Mpc, H*rs in km/s
}

BOSS_DR12_BAO = {
    "label": "BOSS DR12 consensus (BAO-only)",
    "mean_file": "sdss_DR12Consensus_bao.dat",
    "cov_file": "BAO_consensus_covtot_dM_Hz.txt",
    "rs_fid": 147.78,
}

# 6dFGS single point (Gaussian)
SIXDF = {
    "label": "6dFGS",
    "z": 0.106,
    "value": 0.336,
    "error": 0.015,
    "observable": "rs_over_DV",
    "citation": "Beutler et al. (2011)",
}

# SDSS DR7 MGS (1D grid likelihood)
MGS = {
    "label": "SDSS MGS",
    "z": 0.15,
    "prob_file": "sdss_MGS_prob.txt",
    "prob_bounds": (0.8005, 1.1985),
    "rs_rescale": 4.29720761315,  # DV_fid / rs_fid
    "observable": "DV_over_rs",
    "citation": "Ross et al. (2015)",
}


SURVEY_MAP = {
    "desi_dr2": DESI_DR2,
    "desi_dr1": DESI_DR1,
    "eboss_dr16": EBOSS_DR16_LRG,
    "eboss_dr16_qso": EBOSS_DR16_QSO,
    "boss_dr12": BOSS_DR12,
    "boss_dr12_bao": BOSS_DR12_BAO,
    "sixdf": SIXDF,
    "sdss_mgs": MGS,
}


def load_gaussian(dataset_key: str) -> dict:
    """Load a Gaussian BAO dataset.  Returns dict with z, value, observable, cov."""
    info = SURVEY_MAP[dataset_key]

    # special cases
    if dataset_key == "sixdf":
        return {
            "label": info["label"],
            "z": np.array([info["z"]]),
            "value": np.array([info["value"]]),
            "observable": [info["observable"]],
            "cov": np.array([[info["error"] ** 2]]),
        }
    if dataset_key == "sdss_mgs":
        return {
            "label": info["label"],
            "z": np.array([info["z"]]),
            "value": np.array([4.4657]),  # peak of the probability distribution
            "observable": [info["observable"]],
            "cov": np.array([[0.1681 ** 2]]),
            "prob_file": str(DATA_DIR / info["prob_file"]),
        }

    # Standard CobayaSampler format: z value observable
    mean_path = DATA_DIR / info["mean_file"]
    rows = []
    with open(mean_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            z = float(parts[0])
            val = float(parts[1])
            obs = parts[2]
            # Normalise BOSS DR12: convert DM (Mpc) -> DM/rd, bao_Hz_rs -> DH/rd
            if "rs_fid" in info:
                rd = info["rs_fid"]
                if obs == "DM_over_rs":
                    val = val / rd
                    obs = "DM_over_rs"
                elif obs == "bao_Hz_rs":
                    # BOSS DR12 consensus stores H(z) in km/s/Mpc
                    # DH/rd = c / (H(z) * rd_fid)
                    val = 299792.458 / (val * rd)
                    obs = "DH_over_rs"
            rows.append((z, val, obs))

    z_arr = np.array([r[0] for r in rows])
    val_arr = np.array([r[1] for r in rows])
    obs_arr = [r[2] for r in rows]

    cov_path = DATA_DIR / info["cov_file"]
    cov = np.loadtxt(cov_path)
    # BOSS DR12 cov is in original units; need to convert
    if "rs_fid" in info:
        rd = info["rs_fid"]
        npts = len(rows)
        scale = np.ones(npts)
        for k, obs in enumerate(obs_arr):
            if obs == "DM_over_rs":
                scale[k] = 1.0 / rd
            elif obs == "DH_over_rs":
                # DH/rd = c/val_orig -> dDH/dval = -c/val_orig^2
                # Use finite-diff scale; easier: just recompute from the raw cov
                pass
            elif obs == "f_sigma8":
                scale[k] = 1.0  # already dimensionless
        # For DR12, cov is in raw units (DM in Mpc^2, Hz_rs in (km/s)^2);
        # easier to use the cov as-is and match the data vector to original units
        # Actually let's just keep the cov transformation simple
        # Re-read original data and use matching units
        pass  # cov already loaded, DR12 cov is pre-scaled to match data_file units

    # For DR12 the cov is already in the units of the data file
    # Since we transformed data -> DM/rd and DH/rd, we need to transform cov too
    if "rs_fid" in info:
        rd = info["rs_fid"]
        # Re-read raw data to get original values for Jacobian
        raw_rows = []
        with open(mean_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                raw_rows.append((float(parts[0]), float(parts[1]), parts[2]))

        npts = len(raw_rows)
        J = np.eye(npts)
        for k, (_, raw_val, raw_obs) in enumerate(raw_rows):
            if raw_obs == "DM_over_rs":
                J[k, k] = 1.0 / rd  # DM_new = DM_old / rd
            elif raw_obs == "bao_Hz_rs":
                J[k, k] = -299792.458 / (raw_val ** 2 * rd)  # DH/rd = c/(H*rd)
            elif raw_obs == "f_sigma8":
                J[k, k] = 1.0
        cov = J @ cov @ J.T

    return {
        "label": info["label"],
        "z": z_arr,
        "value": val_arr,
        "observable": obs_arr,
        "cov": cov,
    }
