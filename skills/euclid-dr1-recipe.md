---
name: euclid-dr1-recipe
description: What Euclid DR1-Foundation offers here (Wide-Survey photo-z / spec-z galaxy samples, ECLIPSE strong lenses) and what is not released yet
---

# Euclid DR1 recipe

Euclid DR1-Foundation went public on 2026-11-12. It covers ~1900 deg² of the Wide Survey,
with more than 900M objects, ~160M NISP spectra, ~15k strong-lens candidates and ~60M
quasar candidates. This server queries the **ESA Euclid Science Archive** TAP service
(`https://eas.esac.esa.int/tap-server/tap`, `astroquery.esa.euclid` environment `PDR`).
Queries are anonymous; no login is needed.

## Available (data kind: catalog)

| Table (schema `dr1`) | Content |
|---|---|
| `mer_catalogue_wide_survey` | MER photometry, positions, `spurious_flag`, `point_like_prob`, fluxes in µJy |
| `phz_photo_z_wide` | `phz_median`, `phz_mode_1`, `phz_70_int`, `phz_90_int`, `phz_flags`, `phz_classification` (2 = galaxy), `tom_bin_id` |
| `phz_classification_wide` | galaxy / star / QSO probabilities |
| `spectro_z_spe_galaxy_candidates_wide_survey` | NISP `spe_z`, `spe_z_err`, `spe_z_prob`, `spe_rank` (0 = best) |
| `spectro_z_spe_classification_wide_survey` | `spe_class` (galaxy/qso/star), class probabilities |
| `spectro_z_spe_quality_wide_survey` | spectroscopic quality / error flags |

The deep fields use the `*_deep_survey` / `*_deep` tables. All tables join on `object_id`.

Each tool returns:
- `fetch_euclid_catalog(output_dir, sample="spectroscopic", ra=, dec=, radius_deg=)` or `bounds=[ra0, ra1, dec0, dec1]`:
  a CSV/parquet table, a summary JSON (N, z percentiles, the exact ADQL) and an n(z) plot.
- `fetch_euclid_strong_lenses(output_dir, grades=["A","B"])`: the ECLIPSE catalogue
  (Zenodo 19654904, `ECLIPSE_object_catalogue.zip`). It holds 26 662 graded systems:
  6073 A, 8999 B, 9010 C and 2580 X.

## Steps (n(z) of a spectroscopic sample)

1. `fetch_euclid_catalog(out, sample="spectroscopic", ra=80.3, dec=-51.25, radius_deg=0.3)`
   Sanity: ~4.4k galaxies with `spe_z_prob >= 0.99` and VIS < 24.5. The median z is ≈ 1.3,
   and most objects are Hα emitters at 0.9 < z < 1.8.
2. Run the same call with `sample="photometric"` to compare with photo-z. Use a small cone
   (≤ 0.2 deg); the photo-z join is the slow part. Sanity: median `phz_median` is ≈ 0.85.
3. Run `compute_nz_histogram` (LSS family) on the CSV with `z_col="spe_z"` or `"phz_median"`.

## Pitfalls

- `spe_rank = 0` picks the most probable redshift solution. Without `min_spec_z_prob` cuts,
  many rank-0 candidates are unreliable. Keep the threshold at 0.99 for clean samples.
- Magnitudes are computed client-side from `flux_vis_sersic`, and the VIS limit is applied
  as a flux cut. A server-side `LOG10()` aborts the whole query, because the database
  evaluates it on negative fluxes.
- Large regions: queries are capped at 30 deg² and `max_rows`. Use `subsample_mod=N` for a
  deterministic 1-in-N sample, rather than relying on TOP, which returns objects in storage order.
- IRSA mirrors DR1 (doi:10.26131/IRSA674) if the ESA archive is down.

## Not in DR1-Foundation

- Weak-lensing shear catalogues ready for cosmology, and shear 2-pt functions.
- Galaxy-clustering 2-pt measurements, BAO / RSD results, covariances and likelihoods.

These products are expected with **"DR1 Complete" (mid-2027)**. Until then there is no
Euclid summary statistic to fit: use Euclid catalogs for n(z), samples and lens
selections only. For BAO fits use DESI and DES Y6 (`fit_bao_cosmology`).
