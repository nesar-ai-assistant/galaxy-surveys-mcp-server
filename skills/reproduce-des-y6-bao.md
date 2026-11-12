---
name: reproduce-des-y6-bao
description: Reproduce DES Y6 angular BAO (D_M/r_d = 19.51 +/- 0.41 at z=0.851) and combine it with DESI DR2 BAO for Om, H0*rd
---

# Reproduce DES Y6 BAO and combine with DESI DR2

DES Y6 measures α = (D_M/r_d)/(D_M/r_d)_Planck at z_eff = 0.851 from photometric red
galaxies (DES Collaboration, arXiv:2402.10696). The likelihood is **non-Gaussian** and is
released as a Δχ²(α) profile.

## Steps

1. `fetch_des_y6_bao_likelihood(out, variant="fid")`
   Expected: α = 0.9571 (+0.0202/−0.0200) and D_M/r_d = 19.51 ± 0.41.
   The fiducial value is (D_M/r_d)_Planck = 20.385.
2. Check consistency across statistics: `variant="acf"` gives α ≈ 0.951, `"aps"` gives ≈ 0.962
   and `"pcf"` gives ≈ 0.955.
3. Look at the inputs with `fetch_des_y6_clustering(out, statistic="acf")`. It returns
   144 w(θ) points in 6 z-bins plus a 144×144 covariance. Use `"aps"` for 210 points and
   `"pcf"` for 120 points.
4. Joint fit: `fit_bao_cosmology(out, datasets=["desi_dr2", "des_y6_nodesi"], n_grid=100)`
   Expected: Om ≈ 0.297 ± 0.009 and H0*rd ≈ 101.6. The DES-only no-DESI value is
   α = 0.969 ± 0.030 (Mena-Fernandez et al., arXiv:2601.14864).
5. For comparison, `datasets=["desi_dr2", "des_y6"]` gives Om ≈ 0.2966. It comes with a
   warning, because the full DES footprint overlaps DESI and the two are not independent.

## Pitfalls

- Do not add DES Y6 as a Gaussian point. The fit uses the tabulated profile. The
  `fetch_bao_measurements(dataset="des_y6")` row is only a summary.
- DES Y6 alone constrains a single distance ratio, so on its own it is degenerate in
  (Om, H0*rd).
- DES Y6 is ~2σ below Planck LCDM (α < 1). That is the published result; it is not a bug.

## Report

- α and D_M/r_d against the published values.
- Om and H0*rd for DESI DR2 alone, for DESI DR2 + DES Y6 (no-DESI), and the Δχ²
  contributed by DES.
