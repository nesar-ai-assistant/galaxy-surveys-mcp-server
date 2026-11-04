# Conventions — galaxy-surveys-mcp-server

Follows the `cosmo-*-mcp-server` template. See the template's
`docs/CONVENTIONS.md` for the full specification; here is the delta.

## Server identity

| key | value |
|---|---|
| `[project].name` | `galaxy-surveys-mcp-server` |
| `skill_prefix` | `galaxy` |
| port | 8011 |
| server display name | Galaxy Surveys (BAO/RSD + SDSS) MCP Server |

## Tool families

| family | module | weight | notes |
|---|---|---|---|
| bao | `tools.bao` | `fit_bao_cosmology` is heavy | BAO distance and growth-rate data (DESI, BOSS, eBOSS, 6dFGS, MGS) |
| sdss | `tools.sdss` | all light | SDSS SkyServer queries (network required) |
| gswlc | `tools.gswlc` | all light | GSWLC-A2 catalog (downloads on first use) |
| meta | `tools.meta` | all light | Discovery: `list_galaxy_*`, `describe_galaxy_*` |
| dispatch | `mcp_server.dispatch` | N/A | HPC dispatch management |

## Data provenance

- **BAO data** (< 500 KB): bundled in `data/bao/` from
  [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data)
  commit `bb0c1c9`.
- **SDSS**: live queries to SkyServer DR18 via `astroquery.sdss`.
  Results are cached in `~/.cache/cosmo-mcp/galaxy-surveys-mcp-server/sdss/`.
- **GSWLC-A2**: downloaded from https://salims.pages.iu.edu/gswlc/GSWLC-A2.dat.gz
  on first use (~35 MB compressed). Cached in
  `~/.cache/cosmo-mcp/galaxy-surveys-mcp-server/gswlc/`.

## Testing

```bash
# BAO tests (offline, fast):
MPLBACKEND=Agg .venv/bin/python -m pytest tests/test_bao.py -v

# Network tests (SDSS + GSWLC — run once, then cached):
MPLBACKEND=Agg .venv/bin/python -m pytest tests/test_sdss_gswlc.py -v -m network

# All tests:
MPLBACKEND=Agg .venv/bin/python -m pytest -v

# Smoke test (MCP server):
.venv/bin/python tests/smoke_server.py
```

## Sanity checks

- DESI DR2 BAO-only LCDM fit: Om ≈ 0.297 ± 0.009, h*rd ≈ 101.5 Mpc
- DESI DR1 BAO-only: Om ≈ 0.295 ± 0.015
- BOSS DR12 DM/rd values should be 10-20 range, DH/rd 20-25 range
