# galaxy-surveys-mcp-server

MCP server for spectroscopic galaxy surveys: BAO/RSD compilations (DESI DR2/DR1,
eBOSS DR16, BOSS DR12, 6dFGS, SDSS MGS) with cosmological fits and plots, plus
SDSS SkyServer catalog queries, GSWLC-A2 galaxy properties, and large-scale
structure clustering catalogs with correlation function measurements.

Part of the `cosmo-*-mcp-server` family (port 8011).

## Data Tiers

| Tier | Description | Examples |
|------|-------------|----------|
| **1** | Raw / catalog-level data | SDSS SkyServer queries, GSWLC-A2 (640k galaxies), eBOSS/BOSS/DESI clustering catalogs + randoms |
| **2** | Summary statistics with covariance | BAO D/rd, f·σ₈, P(k)/ξ(s) multipoles (Beutler+17, eBOSS DR16 datavecs), n(z), user-measured ξ₀(s), ξ₂(s) |
| **3** | Published inference (comparison only) | Planck LCDM reference curves in BAO plots |

Every tool is decorated `@tier(n)`. Heavy tier 1→2 reductions (`measure_correlation_multipoles`,
`fit_bao_cosmology`) are dispatchable kernels.

## Tool Families

### BAO (`tools/bao/`) — tier 2/3
- `fetch_bao_measurements` — load BAO distance data from bundled files
- `plot_bao_distances` — Hubble diagram (D/rd vs z) with Planck LCDM
- `fit_bao_cosmology` — chi2 grid fit for Om, H0*rd (LCDM or w0waCDM) **[HEAVY]**
- `plot_bao_contours` — Om-H0*rd contour plot from fit grid
- `plot_fsigma8` — f*sigma8(z) growth-rate measurements
- `compare_bao_surveys` — residuals vs Planck LCDM across surveys *(tier 3)*

### SDSS (`tools/sdss/`) — tier 1
- `sdss_sql_query` — execute SQL on SDSS SkyServer DR18
- `fetch_sdss_galaxies` — spectro-photometric galaxy catalog
- `compute_sdss_absolute_magnitudes` — extinction correction + absolute mags + colours
- `plot_color_magnitude` — colour-magnitude diagram (red sequence / blue cloud)
- `plot_sdss_redshift_distribution` — N(z) histogram

### GSWLC (`tools/gswlc/`) — tier 1
- `query_gswlc_catalog` — query GSWLC-A2 with filters (640k galaxies)
- `describe_gswlc_catalog` — column descriptions and filter syntax
- `plot_star_forming_main_sequence` — SFR vs M* density plot
- `crossmatch_sdss_gswlc` — match SDSS with GSWLC on ObjID or sky coords

### LSS Clustering (`tools/lss/`) — tier 1 + 2
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
