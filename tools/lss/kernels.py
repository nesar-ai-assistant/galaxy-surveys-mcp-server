"""Dispatchable correlation-function measurement kernel.

Pure function: JSON-safe args -> JSON-safe dict + files to cwd.
Imports only numpy/scipy/astropy + Corrfunc.
Runs on facility nodes or locally.
"""
from __future__ import annotations

import numpy as np


def measure_xi_multipoles(
    data_ra, data_dec, data_z, data_weight,
    rand_ra, rand_dec, rand_z, rand_weight,
    s_bins_min: float = 0.0,
    s_bins_max: float = 200.0,
    n_s_bins: int = 40,
    cosmology_Om: float = 0.31,
    cosmology_h: float = 0.676,
    use_corrfunc: bool = True,
) -> dict:
    """Measure correlation function multipoles xi_0(s), xi_2(s) from data+randoms.

    Uses Landy-Szalay estimator with FKP/systematic weights.
    Returns dict with s_mid, xi0, xi2, DD, DR, RR counts.

    Parameters
    ----------
    data/rand_ra/dec/z/weight : lists or arrays of float
    s_bins_min/max, n_s_bins : binning in comoving separation s [Mpc/h]
    cosmology_Om, cosmology_h : flat LCDM for comoving distances
    use_corrfunc : if True and available, use Corrfunc; else scipy cKDTree
    """
    data_ra = np.asarray(data_ra, dtype=np.float64)
    data_dec = np.asarray(data_dec, dtype=np.float64)
    data_z = np.asarray(data_z, dtype=np.float64)
    data_weight = np.asarray(data_weight, dtype=np.float64)
    rand_ra = np.asarray(rand_ra, dtype=np.float64)
    rand_dec = np.asarray(rand_dec, dtype=np.float64)
    rand_z = np.asarray(rand_z, dtype=np.float64)
    rand_weight = np.asarray(rand_weight, dtype=np.float64)

    s_edges = np.linspace(s_bins_min, s_bins_max, n_s_bins + 1)
    s_mid = 0.5 * (s_edges[:-1] + s_edges[1:])

    # Comoving distances
    from scipy.integrate import cumulative_trapezoid
    C_KM_S = 299792.458
    z_grid = np.linspace(0, max(data_z.max(), rand_z.max()) * 1.05, 2000)
    Ode = 1.0 - cosmology_Om
    Ez = np.sqrt(cosmology_Om * (1 + z_grid)**3 + Ode)
    chi_grid = cumulative_trapezoid(1.0 / Ez, z_grid, initial=0) * C_KM_S / (100.0 * cosmology_h)  # Mpc/h

    data_chi = np.interp(data_z, z_grid, chi_grid)
    rand_chi = np.interp(rand_z, z_grid, chi_grid)

    # Convert to Cartesian [Mpc/h]
    def radec_to_xyz(ra, dec, chi):
        ra_r = np.deg2rad(ra)
        dec_r = np.deg2rad(dec)
        x = chi * np.cos(dec_r) * np.cos(ra_r)
        y = chi * np.cos(dec_r) * np.sin(ra_r)
        z = chi * np.sin(dec_r)
        return x, y, z

    dx, dy, dz = radec_to_xyz(data_ra, data_dec, data_chi)
    rx, ry, rz = radec_to_xyz(rand_ra, rand_dec, rand_chi)

    # Pair counting with multipole decomposition
    n_mu_bins = 20  # bins in |mu| from 0 to 1

    # Try Corrfunc first
    corrfunc_ok = False
    if use_corrfunc:
        try:
            from Corrfunc.theory.DDsmu import DDsmu
            corrfunc_ok = True
        except ImportError:
            pass

    if corrfunc_ok:
        from Corrfunc.theory.DDsmu import DDsmu
        nthreads = 2

        dd = DDsmu(1, nthreads, s_edges, 1.0, n_mu_bins,
                   dx, dy, dz, weights1=data_weight, weight_type="pair_product",
                   periodic=False, output_savg=True)
        dr = DDsmu(0, nthreads, s_edges, 1.0, n_mu_bins,
                   dx, dy, dz, weights1=data_weight,
                   X2=rx, Y2=ry, Z2=rz, weights2=rand_weight,
                   weight_type="pair_product",
                   periodic=False, output_savg=True)
        rr = DDsmu(1, nthreads, s_edges, 1.0, n_mu_bins,
                   rx, ry, rz, weights1=rand_weight, weight_type="pair_product",
                   periodic=False, output_savg=True)

        # Reshape to (n_s_bins, n_mu_bins) and compute Landy-Szalay per (s, mu)
        dd_w = np.array([d["weightavg"] * d["npairs"] for d in dd]).reshape(n_s_bins, n_mu_bins)
        dr_w = np.array([d["weightavg"] * d["npairs"] for d in dr]).reshape(n_s_bins, n_mu_bins)
        rr_w = np.array([d["weightavg"] * d["npairs"] for d in rr]).reshape(n_s_bins, n_mu_bins)

    else:
        # scipy cKDTree fallback for modest subsets
        from scipy.spatial import cKDTree

        def weighted_pair_counts_smu(x1, y1, z1, w1, x2, y2, z2, w2, s_edges, n_mu, auto=False):
            """Count weighted pairs in (s, mu) bins."""
            pos1 = np.column_stack([x1, y1, z1])
            pos2 = np.column_stack([x2, y2, z2])
            tree2 = cKDTree(pos2)
            counts = np.zeros((len(s_edges) - 1, n_mu), dtype=np.float64)

            # Process in chunks to manage memory
            chunk = 5000
            for i0 in range(0, len(pos1), chunk):
                i1 = min(i0 + chunk, len(pos1))
                pairs = tree2.query_ball_point(pos1[i0:i1], r=s_edges[-1])
                for ii, nbrs in enumerate(pairs):
                    if not nbrs:
                        continue
                    idx = i0 + ii
                    nbrs = np.array(nbrs)
                    if auto:
                        nbrs = nbrs[nbrs > idx]
                        if len(nbrs) == 0:
                            continue
                    dp = pos2[nbrs] - pos1[idx]
                    sep = np.sqrt(np.sum(dp**2, axis=1))
                    # mu = |cos(angle to LOS)|, LOS = midpoint
                    midpoint = 0.5 * (pos1[idx] + pos2[nbrs])
                    los_norm = np.sqrt(np.sum(midpoint**2, axis=1))
                    los_norm = np.maximum(los_norm, 1e-10)
                    mu = np.abs(np.sum(dp * midpoint, axis=1) / (sep * los_norm))
                    mu = np.clip(mu, 0, 1)

                    s_idx = np.searchsorted(s_edges, sep) - 1
                    mu_idx = np.minimum((mu * n_mu).astype(int), n_mu - 1)
                    valid = (s_idx >= 0) & (s_idx < len(s_edges) - 1)
                    ww = w1[idx] * w2[nbrs[valid]]
                    np.add.at(counts, (s_idx[valid], mu_idx[valid]), ww)
            return counts

        dd_w = weighted_pair_counts_smu(dx, dy, dz, data_weight,
                                        dx, dy, dz, data_weight,
                                        s_edges, n_mu_bins, auto=True)
        dr_w = weighted_pair_counts_smu(dx, dy, dz, data_weight,
                                        rx, ry, rz, rand_weight,
                                        s_edges, n_mu_bins, auto=False)
        rr_w = weighted_pair_counts_smu(rx, ry, rz, rand_weight,
                                        rx, ry, rz, rand_weight,
                                        s_edges, n_mu_bins, auto=True)

    # Normalise
    Nd = data_weight.sum()
    Nr = rand_weight.sum()
    f = Nd / Nr

    # Landy-Szalay per (s, mu) bin
    # xi(s, mu) = (DD - 2*f*DR + f^2*RR) / (f^2 * RR)
    rr_safe = np.maximum(rr_w, 1e-30)
    xi_smu = (dd_w - 2 * f * dr_w + f**2 * rr_w) / (f**2 * rr_safe)
    # Zero out bins with no RR
    xi_smu[rr_w == 0] = 0.0

    # Integrate multipoles: xi_ell = (2*ell+1)/2 * integral xi(s,mu) * L_ell(mu) dmu
    mu_edges = np.linspace(0, 1, n_mu_bins + 1)
    mu_mid = 0.5 * (mu_edges[:-1] + mu_edges[1:])
    dmu = mu_edges[1] - mu_edges[0]

    # Legendre polynomials at mu_mid
    from numpy.polynomial.legendre import legval
    L0 = np.ones_like(mu_mid)
    L2 = legval(mu_mid, [0, 0, 1])  # P_2(mu)

    xi0 = np.sum(xi_smu * L0[None, :] * dmu, axis=1) * (2 * 0 + 1)  # ell=0
    xi2 = np.sum(xi_smu * L2[None, :] * dmu, axis=1) * (2 * 2 + 1)  # ell=2

    # Factor of 2 because we integrate from 0 to 1 (half range for even Legendre)
    # Actually for even ell, integral from -1 to 1 of L_ell = 2 * integral from 0 to 1
    # So xi_ell = (2l+1) * integral_0^1 xi(s,mu) L_ell(mu) dmu  (already correct for symmetric)

    return {
        "s_mid": s_mid.tolist(),
        "xi0": xi0.tolist(),
        "xi2": xi2.tolist(),
        "n_data": int(len(data_ra)),
        "n_rand": int(len(rand_ra)),
        "sum_w_data": float(Nd),
        "sum_w_rand": float(Nr),
        "s_bins": s_edges.tolist(),
        "used_corrfunc": corrfunc_ok,
    }
