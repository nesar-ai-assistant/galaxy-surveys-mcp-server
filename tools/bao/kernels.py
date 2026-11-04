"""Dispatchable BAO fitting kernel.

Pure function: JSON-safe args -> JSON-safe dict.
Imports only numpy/scipy/astropy.  Runs on facility nodes or locally.
"""
from __future__ import annotations

import numpy as np

# numpy 2.x renamed trapz -> trapezoid
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)

# --------------- cosmology helpers (flat wCDM / w0waCDM) ---------------

C_KM_S = 299792.458  # km/s


def _Ez(z, Om, w0=-1.0, wa=0.0):
    """E(z) = H(z)/H0 for flat w0wa CDM."""
    a = 1.0 / (1.0 + z)
    Ode = 1.0 - Om
    de = Ode * a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
    return np.sqrt(Om * (1.0 + z) ** 3 + de)


def _comoving_distance(z, Om, H0, w0=-1.0, wa=0.0, npts=500):
    """DC(z) in Mpc via simple trapezoid integration."""
    zz = np.linspace(0, z, npts)
    integrand = 1.0 / _Ez(zz, Om, w0, wa)
    return _trapz(integrand, zz) * C_KM_S / H0


def _DM(z, Om, H0, w0=-1.0, wa=0.0):
    """Comoving transverse distance DM = DC for flat cosmology (Mpc)."""
    return _comoving_distance(z, Om, H0, w0, wa)


def _DH(z, Om, H0, w0=-1.0, wa=0.0):
    """Hubble distance c/H(z) in Mpc."""
    return C_KM_S / (H0 * _Ez(z, Om, w0, wa))


def _DV(z, Om, H0, w0=-1.0, wa=0.0):
    """Volume-averaged distance in Mpc."""
    dm = _DM(z, Om, H0, w0, wa)
    dh = _DH(z, Om, H0, w0, wa)
    return (z * dm ** 2 * dh) ** (1.0 / 3.0)


def _theory_vec(z_arr, obs_arr, Om, H0rd, w0=-1.0, wa=0.0):
    """Compute theory predictions for a vector of (z, observable) pairs.

    H0rd = h * rd in Mpc (h = H0/100).  Range typically 95-105 Mpc.
    The BAO observables DM/rd, DH/rd, DV/rd depend on Om and h*rd only
    (not H0 and rd separately) for flat models.

    DM/rd = (c/100) * int_0^z dz'/E(z') / H0rd
    DH/rd = (c/100) / (E(z) * H0rd)
    DV/rd = (z * DM^2 * DH)^{1/3} / rd
    """
    c100 = C_KM_S / 100.0   # = 2997.92458 Mpc
    theory = np.empty(len(z_arr))
    for i, (z, obs) in enumerate(zip(z_arr, obs_arr, strict=False)):
        if obs == "DM_over_rs":
            # DM/rd = (c/100)/H0rd * int_0^z dz/E(z')
            theory[i] = c100 / H0rd * _comoving_integral(z, Om, w0, wa)
        elif obs == "DH_over_rs":
            # DH/rd = (c/100) / (H0rd * E(z))
            theory[i] = c100 / (H0rd * _Ez(z, Om, w0, wa))
        elif obs == "DV_over_rs":
            dm_rd = c100 / H0rd * _comoving_integral(z, Om, w0, wa)
            dh_rd = c100 / (H0rd * _Ez(z, Om, w0, wa))
            theory[i] = (z * dm_rd ** 2 * dh_rd) ** (1.0 / 3.0)
        else:
            theory[i] = 0.0
    return theory


def _comoving_integral(z, Om, w0=-1.0, wa=0.0, npts=500):
    """Dimensionless comoving integral int_0^z dz'/E(z')."""
    zz = np.linspace(0, z, npts)
    integrand = 1.0 / _Ez(zz, Om, w0, wa)
    return _trapz(integrand, zz)


# --------------- main kernel ---------------

def fit_bao(z_arr, obs_arr, data_vec, cov_mat,
            fit_w0wa=False,
            Om_range=(0.1, 0.6), H0rd_range=(90, 115),
            w0_range=(-2.0, 0.0), wa_range=(-3.0, 2.0),
            n_grid_Om=80, n_grid_H0rd=80, n_grid_w=30,
            cmb_prior=None, sn_prior=None):
    """Chi-squared grid fit of BAO data.

    Parameters
    ----------
    z_arr, obs_arr, data_vec : lists (JSON arrays)
    cov_mat : 2-D list (JSON array of arrays)
    fit_w0wa : bool — if True, scan (Om, H0*rd, w0, wa) with optional priors
    cmb_prior : dict {"Om": [mean,sigma], "H0rd": [mean,sigma]} or None
    sn_prior : dict {"w0": [mean,sigma], "wa": [mean,sigma]} or None
    Returns dict with best_fit, grid arrays, chi2_min, marginalised 1D.
    """
    z_arr = np.asarray(z_arr, dtype=float)
    data_vec = np.asarray(data_vec, dtype=float)
    cov_mat = np.asarray(cov_mat, dtype=float)
    icov = np.linalg.inv(cov_mat)

    Om_grid = np.linspace(*Om_range, n_grid_Om)
    H0rd_grid = np.linspace(*H0rd_range, n_grid_H0rd)

    if not fit_w0wa:
        chi2 = np.full((n_grid_Om, n_grid_H0rd), np.inf)
        for i, Om in enumerate(Om_grid):
            for j, H0rd in enumerate(H0rd_grid):
                th = _theory_vec(z_arr, obs_arr, Om, H0rd)
                diff = data_vec - th
                chi2[i, j] = diff @ icov @ diff
                if cmb_prior:
                    if "Om" in cmb_prior:
                        m, s = cmb_prior["Om"]
                        chi2[i, j] += ((Om - m) / s) ** 2
                    if "H0rd" in cmb_prior:
                        m, s = cmb_prior["H0rd"]
                        chi2[i, j] += ((H0rd - m) / s) ** 2

        imin = np.unravel_index(np.argmin(chi2), chi2.shape)
        Om_best, H0rd_best = Om_grid[imin[0]], H0rd_grid[imin[1]]
        chi2_min = float(chi2[imin])

        # marginalise
        lnL = -0.5 * (chi2 - chi2_min)
        L = np.exp(lnL)
        L_Om = L.sum(axis=1)
        L_Om /= _trapz(L_Om, Om_grid)
        L_H0rd = L.sum(axis=0)
        L_H0rd /= _trapz(L_H0rd, H0rd_grid)

        # 1-sigma from marginalised
        Om_mean = _trapz(Om_grid * L_Om, Om_grid)
        Om_var = _trapz((Om_grid - Om_mean) ** 2 * L_Om, Om_grid)
        H0rd_mean = _trapz(H0rd_grid * L_H0rd, H0rd_grid)
        H0rd_var = _trapz((H0rd_grid - H0rd_mean) ** 2 * L_H0rd, H0rd_grid)

        return {
            "model": "LCDM",
            "best_fit": {"Om": round(float(Om_best), 5),
                         "H0rd": round(float(H0rd_best), 3)},
            "marginalised": {
                "Om": round(float(Om_mean), 5),
                "Om_sigma": round(float(np.sqrt(Om_var)), 5),
                "H0rd": round(float(H0rd_mean), 3),
                "H0rd_sigma": round(float(np.sqrt(H0rd_var)), 3),
            },
            "chi2_min": round(chi2_min, 3),
            "ndof": len(data_vec) - 2,
            "Om_grid": Om_grid.tolist(),
            "H0rd_grid": H0rd_grid.tolist(),
            "chi2_grid": chi2.tolist(),
        }
    else:
        # 4-parameter scan: marginalise w0, wa; keep Om, H0rd grid
        w0_grid = np.linspace(*w0_range, n_grid_w)
        wa_grid = np.linspace(*wa_range, n_grid_w)
        chi2_marg = np.full((n_grid_Om, n_grid_H0rd), np.inf)
        best_overall = {"chi2": np.inf}

        # For each (Om, H0rd), minimise over (w0, wa)
        for i, Om in enumerate(Om_grid):
            for j, H0rd in enumerate(H0rd_grid):
                best_chi2 = np.inf
                for w0 in w0_grid:
                    for wa_v in wa_grid:
                        th = _theory_vec(z_arr, obs_arr, Om, H0rd, w0, wa_v)
                        diff = data_vec - th
                        c2 = float(diff @ icov @ diff)
                        if cmb_prior:
                            if "Om" in cmb_prior:
                                m, s = cmb_prior["Om"]
                                c2 += ((Om - m) / s) ** 2
                            if "H0rd" in cmb_prior:
                                m, s = cmb_prior["H0rd"]
                                c2 += ((H0rd - m) / s) ** 2
                        if sn_prior:
                            if "w0" in sn_prior:
                                m, s = sn_prior["w0"]
                                c2 += ((w0 - m) / s) ** 2
                            if "wa" in sn_prior:
                                m, s = sn_prior["wa"]
                                c2 += ((wa_v - m) / s) ** 2
                        if c2 < best_chi2:
                            best_chi2 = c2
                        if c2 < best_overall["chi2"]:
                            best_overall = {"chi2": c2, "Om": Om, "H0rd": H0rd,
                                            "w0": w0, "wa": wa_v}
                chi2_marg[i, j] = best_chi2

        chi2_min = best_overall["chi2"]
        return {
            "model": "w0waCDM",
            "best_fit": {k: round(float(v), 5) for k, v in best_overall.items()},
            "chi2_min": round(float(chi2_min), 3),
            "ndof": len(data_vec) - 4,
            "Om_grid": Om_grid.tolist(),
            "H0rd_grid": H0rd_grid.tolist(),
            "chi2_grid": chi2_marg.tolist(),
        }
