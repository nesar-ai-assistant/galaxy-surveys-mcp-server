---
name: bao-dark-energy-hints
description: Fit w0waCDM to DESI DR2 BAO with CMB and SN priors to probe dark energy evolution — does w deviate from -1?
---

# BAO dark energy hints (w0wa with priors)

Goal: fit w0waCDM to DESI DR2 BAO data using CMB and SN Gaussian priors.
DESI DR2 + CMB + SN hints at w0 > -1, wa < 0 (time-varying dark energy).

## Steps

1. `fetch_bao_measurements(output_dir, dataset="desi_dr2")`

2. BAO-only LCDM baseline:
   `fit_bao_cosmology(output_dir, dataset="desi_dr2")`
   Record Om, H0*rd.

3. w0waCDM with CMB prior (Planck 2018):
   ```
   fit_bao_cosmology(output_dir, dataset="desi_dr2", fit_w0wa=True,
       cmb_prior='{"Om":[0.3153,0.007],"H0rd":[99.17,0.76]}',
       n_grid=60)
   ```
   CMB prior: Om = 0.3153 +/- 0.007 (Planck 2018), H0*rd = 99.17 +/- 0.76 km/s.

4. w0waCDM with CMB + SN prior:
   ```
   fit_bao_cosmology(output_dir, dataset="desi_dr2", fit_w0wa=True,
       cmb_prior='{"Om":[0.3153,0.007],"H0rd":[99.17,0.76]}',
       sn_prior='{"w0":[-0.9,0.1],"wa":[-0.3,0.3]}',
       n_grid=60)
   ```
   SN prior values are approximate (Pantheon+/DES-SN5YR region).

5. `plot_bao_contours(output_dir, grid_file=..., dataset_label="DESI DR2 + CMB + SN")`

## Sanity checks
- LCDM should give Om ~ 0.297
- w0waCDM without priors: w0 near -1, wa near 0 (BAO alone cannot constrain well)
- With priors: w0 may shift to ~ -0.7 to -0.9, wa ~ -0.3 to -1.0

## Report
- Best-fit w0, wa values
- Whether LCDM (w0=-1, wa=0) is within 2-sigma
- Comparison of chi2 between LCDM and w0waCDM fits
