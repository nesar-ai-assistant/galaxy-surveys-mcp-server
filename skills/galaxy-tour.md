---
name: galaxy-tour
description: Quick tour of all four tool families (BAO, SDSS, GSWLC, LSS clustering) - one call per family, all local, < 10 min
---

# Galaxy Surveys Tour

Quick end-to-end demo exercising every tool family. Use one `output_dir`.

Data tiers: tier 1 = raw catalogs, tier 2 = summary statistics, tier 3 = published
inference (comparison only). Every tool has `@tier(n)`.

## 1. Discovery

```
list_galaxy_datasets()
list_galaxy_skills()
```

## 2. BAO family (tier 2, compare_bao_surveys is tier 3)

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

## 3. SDSS family (tier 1)

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

## 4. GSWLC family (tier 1)

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

## 5. LSS clustering family (tier 1 + tier 2)

### 5a. Fetch clustering catalogs (tier 1)

```
fetch_galaxy_lss_catalog(output_dir, catalog="eboss_dr16_lrg_ngc", kind="data", max_rows=20000)
fetch_galaxy_lss_catalog(output_dir, catalog="eboss_dr16_lrg_ngc", kind="random", max_rows=40000)
```
Expect ~20k data, ~40k randoms with RA, DEC, Z, weights.

### 5b. Measure correlation multipoles (tier 2, HEAVY)

```
measure_correlation_multipoles(output_dir,
    data_file="<output_dir>/eboss_dr16_lrg_ngc_data.csv",
    random_file="<output_dir>/eboss_dr16_lrg_ngc_random.csv",
    s_max=200, n_bins=40)
```
Sanity: xi_0 should show a BAO bump near s ~ 100-110 Mpc/h.
For BOSS CMASS: xi_0(s=20 Mpc/h) should be of order 0.1-0.3.

### 5c. n(z) histogram (tier 2)

```
compute_nz_histogram(output_dir, catalog_file="<output_dir>/eboss_dr16_lrg_ngc_data.csv", n_bins=30)
```

### 5d. Published multipoles (tier 2)

```
fetch_published_multipoles(output_dir, dataset="boss_dr12_pk_z1_ngc")
```
Expect Beutler+17 P_0, P_2, P_4 with ~120 k-bins and covariance.

```
fetch_published_multipoles(output_dir, dataset="eboss_dr16_lrg_xi")
```
Expect eBOSS DR16 LRG xi_0, xi_2 with covariance.

### 5e. Plot multipoles (tier 2)

```
plot_clustering_multipoles(output_dir,
    data_file="<output_dir>/xi_multipoles.csv",
    label="eBOSS LRG (measured)",
    published_file="<output_dir>/published_eboss_dr16_lrg_xi_data.csv",
    published_label="eBOSS DR16 (official)",
    plot_type="xi")
```

```
plot_clustering_multipoles(output_dir,
    data_file="<output_dir>/published_boss_dr12_pk_z1_ngc_data.csv",
    label="Beutler+17 z1",
    plot_type="pk")
```

## Report

- BAO fit Om and H0*rd values
- Number of galaxies fetched from SDSS and GSWLC
- Whether CMD shows bimodality
- Whether SFMS shows the expected tight main sequence
- LSS: number of data+random galaxies fetched
- Whether xi_0 shows a BAO bump near s ~ 100 Mpc/h
- n(z) histogram shape matches survey's published n(z)
- Published multipoles: citation and number of bins
