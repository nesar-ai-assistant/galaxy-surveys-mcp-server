---
name: clustering-measurement-recipe
description: Measure galaxy clustering multipoles from raw catalogs - end-to-end recipe from tier-1 download through tier-2 xi(s) measurement, comparison with published results
---

# Clustering Measurement Recipe

Step-by-step recipe for measuring galaxy 2-point correlation function multipoles
from public survey catalogs and comparing with published results.

## Prerequisites

Available catalogs: `eboss_dr16_lrg_ngc`, `eboss_dr16_qso_ngc`,
`eboss_dr16_elg_ngc`, `boss_dr12_cmass_north`, `boss_dr12_lowz_north`,
`desi_dr1_lrg_ngc`, `desi_dr1_elg_ngc`, `desi_dr1_qso_ngc`.

## 1. Fetch data and random catalogs (tier 1)

```
fetch_galaxy_lss_catalog(output_dir, catalog="boss_dr12_cmass_north", kind="data", max_rows=30000)
fetch_galaxy_lss_catalog(output_dir, catalog="boss_dr12_cmass_north", kind="random", max_rows=60000, random_fraction=0.1)
```
Use 2x more randoms than data. For large catalogs, subsample randoms with
`random_fraction` to keep pair-counting tractable (< 100k total objects).
Weight columns are combined into WEIGHT_TOTAL automatically.

## 2. Check n(z) (tier 2)

```
compute_nz_histogram(output_dir, catalog_file="<data csv>", n_bins=30)
```
Verify the redshift distribution matches the expected survey range.

## 3. Measure multipoles (tier 2, HEAVY)

```
measure_correlation_multipoles(output_dir,
    data_file="<data csv>", random_file="<random csv>",
    s_max=200, n_bins=40)
```
Uses Corrfunc if available, scipy cKDTree fallback otherwise.
Can be dispatched to HPC via `set_dispatch("polaris")`.

### Sanity checks:
- **BOSS CMASS**: BAO bump near s ~ 100-110 Mpc/h in xi_0;
  xi_0(s=20 Mpc/h) of order 0.1-0.3
- **eBOSS LRG**: similar BAO feature but noisier for small subsets
- **DESI LRG**: strong BAO signal with larger sample

## 4. Compare with published multipoles (tier 2)

```
fetch_published_multipoles(output_dir, dataset="boss_dr12_pk_z1_ngc")
# or for xi(s):
fetch_published_multipoles(output_dir, dataset="eboss_dr16_lrg_xi")
```

Available published datasets:
- `boss_dr12_pk_z1_ngc`: Beutler+17 P(k) z=0.38 NGC (monopole + quadrupole + hexadecapole + covariance)
- `boss_dr12_pk_z3_ngc`: Beutler+17 P(k) z=0.61 NGC
- `eboss_dr16_lrg_xi`: eBOSS DR16 LRG xi(s) NGC+SGC (monopole + quadrupole + covariance)
- `eboss_dr16_qso_xi`: eBOSS DR16 QSO xi(s) NGC+SGC

## 5. Plot comparison (tier 2)

```
plot_clustering_multipoles(output_dir,
    data_file="<xi_multipoles.csv>",
    label="My measurement",
    published_file="<published csv>",
    published_label="Beutler+17 (official)",
    plot_type="xi")
```

## Pitfalls

- **Subset size**: With < 10k galaxies, the BAO bump will be noisy. Use ≥ 20k
  data + 40k randoms for a recognizable feature.
- **Random fraction**: eBOSS/BOSS random catalogs have 50x the data; subsample
  with `random_fraction=0.05` to keep computation feasible.
- **Weight columns**: BOSS/eBOSS use WEIGHT_SYSTOT × WEIGHT_CP × WEIGHT_NOZ × WEIGHT_FKP;
  DESI uses WEIGHT × WEIGHT_FKP. All handled automatically.
- **Corrfunc vs cKDTree**: Corrfunc is ~10x faster but requires macOS x86_64
  or Linux. The cKDTree fallback is correct but slower for > 50k objects.
