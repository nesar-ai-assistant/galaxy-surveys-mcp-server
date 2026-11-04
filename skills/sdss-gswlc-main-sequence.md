---
name: sdss-gswlc-main-sequence
description: Galaxy population study combining SDSS photometry with GSWLC SFRs — CMD bimodality and star-forming main sequence
---

# SDSS + GSWLC galaxy population study

Goal: fetch an SDSS galaxy sample, compute absolute magnitudes and colours,
cross-match with GSWLC for stellar masses and SFRs, and produce the classic
colour-magnitude diagram and star-forming main sequence.

## Steps

1. `fetch_sdss_galaxies(output_dir, limit=2000, min_z=0.02, max_z=0.08, clean_photometry=True)`

2. `compute_sdss_absolute_magnitudes(output_dir, input_file="<output_dir>/sdss_galaxies.csv")`

3. `plot_color_magnitude(output_dir, input_file="<output_dir>/sdss_absolute_mags.csv")`
   Sanity: bimodal g-r distribution. Red sequence at g-r ~ 0.7-0.9, M_r ~ -20 to -22.
   Blue cloud at g-r ~ 0.3-0.5, M_r ~ -18 to -21.

4. `plot_sdss_redshift_distribution(output_dir, input_file="<output_dir>/sdss_galaxies.csv")`

5. `query_gswlc_catalog(output_dir, filters='{"z":{"min":0.02,"max":0.08},"flag_mgs":{"eq":1}}', limit=50000)`

6. `plot_star_forming_main_sequence(output_dir, input_file="<output_dir>/gswlc_query.csv")`
   Sanity: tight main sequence with logSFR ~ 0.76 * logM* - 7.64.

7. `crossmatch_sdss_gswlc(output_dir, sdss_file="<output_dir>/sdss_galaxies.csv")`

## Report
- Number of galaxies in SDSS sample and GSWLC match
- Whether CMD shows bimodality
- Main sequence slope (should be ~ 0.7-0.8 dex/dex)
- Match fraction between SDSS and GSWLC
