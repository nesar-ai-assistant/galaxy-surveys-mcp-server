---
name: galaxy-tour
description: Quick tour of all three tool families (BAO, SDSS, GSWLC) - one call per family, all local, < 5 min
---

# Galaxy Surveys Tour

Quick end-to-end demo exercising every tool family. Use one `output_dir`.

## 1. Discovery

```
list_galaxy_datasets()
list_galaxy_skills()
```

## 2. BAO family

```
fetch_bao_measurements(output_dir, dataset="desi_dr2")
```
Expect 13 data points (7 tracers, DM/rd + DH/rd + DV/rd).

```
plot_bao_distances(output_dir, datasets="desi_dr2,desi_dr1,boss_dr12_bao,sixdf,sdss_mgs")
```
Check: Planck LCDM curve threads through data points within ~2 sigma.

```
fit_bao_cosmology(output_dir, dataset="desi_dr2", n_grid=60)
```
Sanity check: Om ~ 0.297 +/- 0.009, H0*rd ~ 101.5 km/s (DESI DR2 BAO-only).

```
plot_bao_contours(output_dir, grid_file="<output_dir>/bao_fit_desi_dr2_grid.npz", dataset_label="DESI DR2")
```

## 3. SDSS family

```
fetch_sdss_galaxies(output_dir, limit=500, min_z=0.02, max_z=0.08)
```
Expect ~500 galaxies with ugriz photometry.

```
compute_sdss_absolute_magnitudes(output_dir, input_file="<output_dir>/sdss_galaxies.csv")
```

```
plot_color_magnitude(output_dir, input_file="<output_dir>/sdss_absolute_mags.csv")
```
Check: bimodal distribution with red sequence (g-r ~ 0.8) and blue cloud (g-r ~ 0.4).

```
plot_sdss_redshift_distribution(output_dir, input_file="<output_dir>/sdss_galaxies.csv")
```

## 4. GSWLC family

```
describe_gswlc_catalog()
```

```
query_gswlc_catalog(output_dir, filters='{"logMstar":{"min":9.5,"max":11.5},"flag_mgs":{"eq":1}}', limit=5000)
```

```
plot_star_forming_main_sequence(output_dir, input_file="<output_dir>/gswlc_query.csv")
```
Check: clear main sequence with slope ~0.76, green valley visible, quenched galaxies at logSFR < -1.

## Report

- BAO fit Om and H0*rd values
- Number of galaxies fetched from SDSS and GSWLC
- Whether CMD shows bimodality
- Whether SFMS shows the expected tight main sequence
