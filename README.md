# galaxy-surveys-mcp-server

MCP server for spectroscopic galaxy surveys: BAO/RSD compilations (DESI DR2/DR1,
eBOSS DR16, BOSS DR12, 6dFGS, SDSS MGS) with cosmological fits and plots, plus
SDSS SkyServer catalog queries and GSWLC-A2 galaxy properties.

Part of the `cosmo-*-mcp-server` family (port 8011).

## Tool Families

### BAO (`tools/bao/`)
- `fetch_bao_measurements` — load BAO distance data from bundled files
- `plot_bao_distances` — Hubble diagram (D/rd vs z) with Planck LCDM
- `fit_bao_cosmology` — chi2 grid fit for Om, H0*rd (LCDM or w0waCDM) **[HEAVY]**
- `plot_bao_contours` — Om-H0*rd contour plot from fit grid
- `plot_fsigma8` — f*sigma8(z) growth-rate measurements
- `compare_bao_surveys` — residuals vs Planck LCDM across surveys

### SDSS (`tools/sdss/`)
- `sdss_sql_query` — execute SQL on SDSS SkyServer DR18
- `fetch_sdss_galaxies` — spectro-photometric galaxy catalog
- `compute_sdss_absolute_magnitudes` — extinction correction + absolute mags + colours
- `plot_color_magnitude` — colour-magnitude diagram (red sequence / blue cloud)
- `plot_sdss_redshift_distribution` — N(z) histogram

### GSWLC (`tools/gswlc/`)
- `query_gswlc_catalog` — query GSWLC-A2 with filters (640k galaxies)
- `describe_gswlc_catalog` — column descriptions and filter syntax
- `plot_star_forming_main_sequence` — SFR vs M* density plot
- `crossmatch_sdss_gswlc` — match SDSS with GSWLC on ObjID or sky coords

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
