"""Single source of truth for the datasets this server exposes.

Keys are short dataset ids used by every tool's ``dataset`` argument.
"access" is "bundled" (file in data/), "download" (fetched on first use
into the data cache) or "remote-query" (live service).

Products are split by data tier:
  tier1 = raw / catalog-level (per-object catalogs)
  tier2 = processed summary statistics measured from tier-1 data
  tier3 = published inference (comparison only)
"""

DATASETS: dict[str, dict] = {
    # ---- BAO / RSD distance measurements (tier 2) ----
    "desi_dr2": {
        "title": "DESI Data Release 2 BAO (2025)",
        "products": {
            "tier2": ["DM/rd", "DH/rd", "DV/rd", "covariance"],
        },
        "access": "bundled",
        "urls": ["https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2"],
        "citation": "DESI Collaboration, arXiv:2503.14738, arXiv:2503.14739 (2025)",
        "tools": ["fetch_bao_measurements", "plot_bao_distances", "fit_bao_cosmology",
                  "plot_bao_contours", "compare_bao_surveys"],
    },
    "desi_dr1": {
        "title": "DESI Data Release 1 BAO (2024)",
        "products": {
            "tier2": ["DM/rd", "DH/rd", "DV/rd", "covariance"],
        },
        "access": "bundled",
        "urls": ["https://github.com/CobayaSampler/bao_data"],
        "citation": "DESI Collaboration, arXiv:2404.03002 (2024)",
        "tools": ["fetch_bao_measurements", "plot_bao_distances", "fit_bao_cosmology",
                  "plot_bao_contours", "compare_bao_surveys"],
    },
    "eboss_dr16": {
        "title": "eBOSS DR16 BAO + full-shape + clustering catalogs (LRG, QSO, ELG, Lya)",
        "products": {
            "tier1": ["LRG/QSO/ELG NGC clustering catalogs + randoms"],
            "tier2": ["DM/rd", "DH/rd", "DV/rd", "f*sigma8", "covariance",
                      "xi(s) multipoles (LRG, QSO)", "P(k) multipoles"],
        },
        "access": "bundled (BAO) / download (catalogs, multipoles)",
        "urls": ["https://github.com/CobayaSampler/bao_data",
                 "https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/",
                 "https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/dataveccov/"],
        "citation": "Alam et al., Phys. Rev. D 103, 083533 (2021), arXiv:2007.08991",
        "tools": ["fetch_bao_measurements", "plot_bao_distances", "plot_fsigma8",
                  "fit_bao_cosmology", "compare_bao_surveys",
                  "fetch_galaxy_lss_catalog", "measure_correlation_multipoles",
                  "fetch_published_multipoles", "plot_clustering_multipoles",
                  "compute_nz_histogram"],
    },
    "boss_dr12": {
        "title": "BOSS DR12 consensus BAO + full-shape + clustering catalogs",
        "products": {
            "tier1": ["CMASS/LOWZ North clustering catalogs + randoms"],
            "tier2": ["DM (Mpc)", "H*rs (km/s)", "f*sigma8", "covariance",
                      "P(k) multipoles (Beutler+17)"],
        },
        "access": "bundled (BAO) / download (catalogs, multipoles)",
        "urls": ["https://github.com/CobayaSampler/bao_data",
                 "https://data.sdss.org/sas/dr12/boss/lss/",
                 "https://data.sdss.org/sas/dr12/boss/papers/clustering/"],
        "citation": "Alam et al., MNRAS 470, 2617 (2017), arXiv:1607.03155",
        "notes": "rs_fid = 147.78 Mpc; DM stored in Mpc, H in km/s * rs_fid",
        "tools": ["fetch_bao_measurements", "plot_bao_distances", "plot_fsigma8",
                  "fit_bao_cosmology", "compare_bao_surveys",
                  "fetch_galaxy_lss_catalog", "measure_correlation_multipoles",
                  "fetch_published_multipoles", "plot_clustering_multipoles",
                  "compute_nz_histogram"],
    },
    "sdss_mgs": {
        "title": "SDSS DR7 Main Galaxy Sample BAO (z=0.15)",
        "products": {
            "tier2": ["DV/rd 1D probability"],
        },
        "access": "bundled",
        "urls": ["https://github.com/CobayaSampler/bao_data"],
        "citation": "Ross et al., MNRAS 449, 835 (2015), arXiv:1409.3242",
        "tools": ["fetch_bao_measurements", "plot_bao_distances"],
    },
    "sixdf": {
        "title": "6dF Galaxy Survey BAO (z=0.106)",
        "products": {
            "tier2": ["rs/DV point measurement"],
        },
        "access": "bundled",
        "urls": [],
        "citation": "Beutler et al., MNRAS 416, 3017 (2011), arXiv:1106.3366",
        "notes": "Single Gaussian point: rs/DV = 0.336 +/- 0.015 at z=0.106",
        "tools": ["fetch_bao_measurements", "plot_bao_distances"],
    },
    # ---- DESI DR1 LSS clustering catalogs (tier 1) ----
    "desi_dr1_lss": {
        "title": "DESI DR1 LSS clustering catalogs (LRG, ELG, QSO, NGC)",
        "products": {
            "tier1": ["LRG/ELG/QSO NGC clustering catalogs + randoms (v1.5)",
                      "official n(z)"],
            "tier2": ["n(z) histogram", "xi(s) multipoles (user-measured)"],
        },
        "access": "download (HTTP Range, row-limited subsets)",
        "urls": ["https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/"],
        "citation": "DESI Collaboration, arXiv:2306.06308 (2023)",
        "tools": ["fetch_galaxy_lss_catalog", "measure_correlation_multipoles",
                  "compute_nz_histogram", "plot_clustering_multipoles"],
    },
    # ---- SDSS SkyServer (tier 1) ----
    "sdss_skyserver": {
        "title": "SDSS SkyServer SQL (DR18)",
        "products": {
            "tier1": ["spectroscopy", "photometry", "galaxy catalogs"],
        },
        "access": "remote-query",
        "urls": ["https://skyserver.sdss.org/dr18/"],
        "citation": "Abdurro'uf et al., ApJS 267, 44 (2023)",
        "tools": ["sdss_sql_query", "fetch_sdss_galaxies",
                  "compute_sdss_absolute_magnitudes",
                  "plot_color_magnitude", "plot_sdss_redshift_distribution"],
    },
    # ---- GSWLC (tier 1) ----
    "gswlc_a2": {
        "title": "GALEX-SDSS-WISE Legacy Catalog A2 (GSWLC-A2)",
        "products": {
            "tier1": ["stellar mass", "SFR", "dust attenuation", "640k galaxies"],
        },
        "access": "download",
        "urls": ["https://salims.pages.iu.edu/gswlc/GSWLC-A2.dat.gz"],
        "citation": "Salim et al., ApJS 227, 2 (2016); Salim, Boquien & Lee, ApJ 859, 11 (2018)",
        "tools": ["query_gswlc_catalog", "describe_gswlc_catalog",
                  "plot_star_forming_main_sequence", "crossmatch_sdss_gswlc"],
    },
}
