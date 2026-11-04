---
name: reproduce-desi-dr2-bao-om
description: Reproduce the DESI DR2 BAO-only LCDM constraint on Omega_m and H0*rd — expected Om = 0.295-0.300, H0*rd ~ 101.5
---

# Reproduce DESI DR2 BAO-only Om constraint

Goal: fit flat LCDM (Om, H0*rd) to DESI DR2 BAO data and verify the
published result of Om ≈ 0.297 ± 0.009.

## Steps

1. `fetch_bao_measurements(output_dir, dataset="desi_dr2")`
   Sanity: 13 data points.

2. `fit_bao_cosmology(output_dir, dataset="desi_dr2", n_grid=100)`
   Expected: Om ≈ 0.297 ± 0.009, H0*rd ≈ 101.5 Mpc.
   Chi2/dof should be near 1.

3. `plot_bao_contours(output_dir, grid_file="<output_dir>/bao_fit_desi_dr2_grid.npz", dataset_label="DESI DR2")`

4. Comparison: `fit_bao_cosmology(output_dir, dataset="desi_dr1", n_grid=100)`
   DR1 expected: Om ≈ 0.295 ± 0.015. The larger error confirms DR2's
   improved precision.

5. `compare_bao_surveys(output_dir, datasets="desi_dr2,desi_dr1")`

## Report
- Om and H0*rd best-fit and marginalised values for DR2 and DR1
- Chi2/dof for both
- Whether the values are consistent with the published DESI DR2 results
- Any > 2-sigma residuals in the comparison plot
