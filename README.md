# galaxy-surveys-mcp-server

MCP server for spectroscopic galaxy surveys: BAO/RSD compilations (DESI DR2/DR1,
eBOSS DR16, BOSS DR12, 6dFGS, SDSS MGS) with cosmological fits and plots, plus
SDSS SkyServer catalog queries, GSWLC-A2 galaxy properties, and large-scale
structure clustering catalogs with correlation function measurements.

Part of the `cosmo-*-mcp-server` family (port 8011).

## Data products

| Catalogs | Summary statistics | Published results (for comparison only) |
|----------|--------------------|-----------------------------------------|
| SDSS SkyServer queries, GSWLC-A2 (640k galaxies), eBOSS/BOSS/DESI clustering catalogs (data + randoms) | BAO D/rd, f·σ₈, P(k)/ξ(s) multipoles with covariance (Beutler+17, eBOSS DR16 datavecs), n(z), user-measured ξ₀(s), ξ₂(s) | Planck ΛCDM reference curves in BAO plots |

Every tool declares the kind of data it works on with `@data_kind(...)` (`catalog`,
`summary_statistic`, `published_result`; see `docs/CONVENTIONS.md`). Heavy reductions
(from catalogs to correlation-function multipoles in `measure_correlation_multipoles`, and the
`fit_bao_cosmology` grid) are dispatchable kernels.

## Tool Families

### BAO (`tools/bao/`) — summary statistics
- `fetch_bao_measurements` — load BAO distance data from bundled files
- `plot_bao_distances` — Hubble diagram (D/rd vs z) with Planck LCDM
- `fit_bao_cosmology` — chi2 grid fit for Om, H0*rd (LCDM or w0waCDM) **[HEAVY]**
- `plot_bao_contours` — Om-H0*rd contour plot from fit grid
- `plot_fsigma8` — f*sigma8(z) growth-rate measurements
- `compare_bao_surveys` — residuals vs Planck LCDM across surveys *(published Planck reference, for comparison only)*
- `fetch_des_y6_bao_likelihood` — DES Y6 angular BAO Δχ²(α) profile (non-Gaussian) at z_eff = 0.851:
  consensus / ACF / APS / PCF and the DESI-overlap-removed (`nodesi_dr1tiles`, `nodesi_deccut`) variants
- `fetch_des_y6_clustering` — DES Y6 BAO w(θ), C_ℓ, ξ_p(s_⊥) data vectors (6 z-bins) + covariance
  (downloaded from the DES public release, cached)

Joint fits: `fit_bao_cosmology(datasets=["desi_dr2", "des_y6_nodesi"])` combines the
Gaussian DESI likelihood with the exact DES Y6 χ²(α) profile. DES Y6 overlaps the DESI
footprint, so use `des_y6_nodesi` with DESI. `des_y6` with DESI still runs but returns a warning.

### Euclid DR1 (`tools/euclid/`) — galaxy catalogs
- `fetch_euclid_catalog` — photometric (MER × photo-z) or spectroscopic (MER × NISP z) Wide-Survey
  galaxies in a cone or box. Queries the ESA Euclid Science Archive anonymously via `astroquery.esa.euclid`,
  with a plain TAP fallback. Quality cuts: VIS mag, galaxy class, `spe_z_prob`, `phz_flags`.
  Server-side TOP plus a deterministic 1-in-N subsample. Outputs CSV/parquet, a summary JSON with the ADQL, and n(z).
- `fetch_euclid_strong_lenses` — ECLIPSE DR1 strong-lens catalogue (Zenodo 19654904; 26 662 graded
  systems), filtered by expert grade/score, ensemble score, DESI deflector z, and cone.

DR1-Foundation (2026-11-12) holds **catalogues only**. Shear 2-pt functions, clustering
measurements and covariances are expected with "DR1 Complete" (mid-2027).

### SDSS (`tools/sdss/`) — galaxy catalogs
- `sdss_sql_query` — execute SQL on SDSS SkyServer DR18
- `fetch_sdss_galaxies` — spectro-photometric galaxy catalog
- `compute_sdss_absolute_magnitudes` — extinction correction + absolute mags + colours
- `plot_color_magnitude` — colour-magnitude diagram (red sequence / blue cloud)
- `plot_sdss_redshift_distribution` — N(z) histogram

### GSWLC (`tools/gswlc/`) — galaxy property catalog
- `query_gswlc_catalog` — query GSWLC-A2 with filters (640k galaxies)
- `describe_gswlc_catalog` — column descriptions and filter syntax
- `plot_star_forming_main_sequence` — SFR vs M* density plot
- `crossmatch_sdss_gswlc` — match SDSS with GSWLC on ObjID or sky coords

### LSS Clustering (`tools/lss/`) — clustering catalogs and multipoles
- `fetch_galaxy_lss_catalog` — download subsets of eBOSS DR16 / BOSS DR12 / DESI DR1
  clustering catalogs (data or randoms, row-limited, z-filtered) via HTTP Range
- `measure_correlation_multipoles` — Landy-Szalay ξ₀(s), ξ₂(s) from data+randoms
  with FKP/systematic weights (Corrfunc or scipy cKDTree fallback) **[HEAVY]**
- `compute_nz_histogram` — redshift distribution n(z) from a catalog
- `fetch_published_multipoles` — released P(k) (Beutler+17) or ξ(s) (eBOSS DR16) with covariance
- `plot_clustering_multipoles` — plot ξ(s) or P(k) multipoles with optional published comparison

**Catalogs available:**
| Key | Survey | z range | Weight columns |
|-----|--------|---------|----------------|
| `eboss_dr16_lrg_ngc` | eBOSS DR16 LRG NGC | 0.6–1.0 | WEIGHT_SYSTOT × CP × NOZ × FKP |
| `eboss_dr16_qso_ngc` | eBOSS DR16 QSO NGC | 0.8–2.2 | same |
| `eboss_dr16_elg_ngc` | eBOSS DR16 ELG NGC | 0.6–1.1 | same |
| `boss_dr12_cmass_north` | BOSS DR12 CMASS N | 0.43–0.7 | same |
| `boss_dr12_lowz_north` | BOSS DR12 LOWZ N | 0.15–0.43 | same |
| `desi_dr1_lrg_ngc` | DESI DR1 LRG NGC | 0.4–1.1 | WEIGHT × WEIGHT_FKP |
| `desi_dr1_elg_ngc` | DESI DR1 ELG NGC | 0.8–1.6 | same |
| `desi_dr1_qso_ngc` | DESI DR1 QSO NGC | 0.8–2.1 | same |

## Quick Start

```bash
cd galaxy-surveys-mcp-server
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
MPLBACKEND=Agg .venv/bin/python -m pytest -q
.venv/bin/python tests/smoke_server.py
```

## Sample Figures

See `docs/figures/` for example outputs from the tour.

## Data

BAO data from [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data)
(bundled, < 500 KB). SDSS queries go to SkyServer DR18 (network). GSWLC-A2
downloads on first use from https://salims.pages.iu.edu/gswlc/ (~35 MB).
LSS clustering catalogs download on demand from SDSS SAS (DR12/DR16) and
DESI data server (DR1); only subsets are fetched (row-limited, z-range filtered).
Published multipoles from BOSS DR12 (Beutler et al. 2017) and eBOSS DR16
(Alam et al. 2021) data vector archives.

DES Y6 BAO comes from the DES public release
(https://desdr-server.ncsa.illinois.edu/despublic/y6a2_files/y6_bao/; DES Collaboration,
arXiv:2402.10696). The χ²(α) profiles are bundled in `data/bao/des_y6/`, including the no-DESI
profiles of Mena-Fernandez et al., arXiv:2601.14864. Data vectors and covariances download on demand.
Euclid DR1 is queried live from https://eas.esac.esa.int/tap-server/tap.
ECLIPSE downloads once from Zenodo (3.7 MB).
