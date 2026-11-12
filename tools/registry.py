"""Single source of truth for the datasets this server exposes.

Keys are short dataset ids used by every tool's ``dataset`` argument.
"access" is "bundled" (file in data/), "download" (fetched on first use
into the data cache) or "remote-query" (live service).

Products are grouped by kind of data:
  catalog           = object-level data (galaxy catalogs, clustering catalogs + randoms)
  summary_statistic = measurements reduced from catalogs, with covariance
  published_result  = published inference products (comparison only)
"""

DATASETS: dict[str, dict] = {
    # ---- BAO / RSD distance measurements (summary statistics) ----
    "desi_dr2": {
        "title": "DESI Data Release 2 BAO (2025)",
        "products": {
            "summary_statistic": ["DM/rd", "DH/rd", "DV/rd", "covariance"],
        },
        "access": "bundled",
        "urls": ["https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2"],
        "citation": "DESI Collaboration, arXiv:2503.14738, arXiv:2503.14739 (2025)",
        "tools": ["fetch_bao_measurements", "plot_bao_distances", "fit_bao_cosmology",
                  "plot_bao_contours", "compare_bao_surveys"],
        "notes": "Joint fits: fit_bao_cosmology(datasets=['desi_dr2','des_y6_nodesi']).",
    },
    "desi_dr1": {
        "title": "DESI Data Release 1 BAO (2024)",
        "products": {
            "summary_statistic": ["DM/rd", "DH/rd", "DV/rd", "covariance"],
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
            "catalog": ["LRG/QSO/ELG NGC clustering catalogs + randoms"],
            "summary_statistic": ["DM/rd", "DH/rd", "DV/rd", "f*sigma8", "covariance",
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
            "catalog": ["CMASS/LOWZ North clustering catalogs + randoms"],
            "summary_statistic": ["DM (Mpc)", "H*rs (km/s)", "f*sigma8", "covariance",
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
            "summary_statistic": ["DV/rd 1D probability"],
        },
        "access": "bundled",
        "urls": ["https://github.com/CobayaSampler/bao_data"],
        "citation": "Ross et al., MNRAS 449, 835 (2015), arXiv:1409.3242",
        "tools": ["fetch_bao_measurements", "plot_bao_distances"],
    },
    "sixdf": {
        "title": "6dF Galaxy Survey BAO (z=0.106)",
        "products": {
            "summary_statistic": ["rs/DV point measurement"],
        },
        "access": "bundled",
        "urls": [],
        "citation": "Beutler et al., MNRAS 416, 3017 (2011), arXiv:1106.3366",
        "notes": "Single Gaussian point: rs/DV = 0.336 +/- 0.015 at z=0.106",
        "tools": ["fetch_bao_measurements", "plot_bao_distances"],
    },
    "des_y6": {
        "title": "DES Y6 angular BAO (photometric red galaxies, z_eff = 0.851)",
        "products": {
            "summary_statistic": [
                "alpha = (D_M/r_d)/(D_M/r_d)_Planck chi2 profile (non-Gaussian; ACF+APS+PCF consensus "
                "and per-statistic)",
                "D_M/r_d(z=0.851) = 19.51 +/- 0.41",
                "w(theta), C_ell, xi_p(s_perp) data vectors in 6 z-bins + covariances (Planck data-like, MICE mock-like)",
                "no-DESI chi2 profiles (DESI DR1 tiles removed / Dec cut) for joint fits with DESI",
            ],
        },
        "access": "bundled (chi2 profiles) / download (data vectors, covariances)",
        "urls": ["https://desdr-server.ncsa.illinois.edu/despublic/y6a2_files/y6_bao/",
                 "https://des.ncsa.illinois.edu/releases/y6a2/Y6bao",
                 "https://github.com/cosmosis-developers/cosmosis-standard-library/tree/main/likelihood/bao/des"],
        "citation": "DES Collaboration, arXiv:2402.10696 (2024); no-DESI: Mena-Fernandez et al., "
                    "arXiv:2601.14864, Phys. Rev. D 113, 123521 (2026)",
        "notes": "BAO keys 'des_y6' (full footprint) and 'des_y6_nodesi' (DESI DR1 tiles removed; use with DESI "
                 "DR1/DR2). Catalogs, randoms, mask, n(z), templates and mocks are in the same release (not wrapped).",
        "tools": ["fetch_des_y6_bao_likelihood", "fetch_des_y6_clustering", "fetch_bao_measurements",
                  "fit_bao_cosmology", "plot_bao_distances"],
    },
    # ---- Euclid DR1 (galaxy catalogs, strong lenses) ----
    "euclid_dr1": {
        "title": "Euclid DR1-Foundation (2026-11-12): Wide Survey galaxy catalogs",
        "products": {
            "catalog": ["MER photometry (dr1.mer_catalogue_wide_survey)",
                        "photometric redshifts (dr1.phz_photo_z_wide)",
                        "NISP spectroscopic redshifts + classes (dr1.spectro_z_spe_*_wide_survey)"],
        },
        "access": "remote-query (ESA Euclid Science Archive TAP, anonymous; astroquery.esa.euclid PDR)",
        "urls": ["https://eas.esac.esa.int/tap-server/tap", "https://easidr.esac.esa.int/sas/",
                 "https://irsa.ipac.caltech.edu (mirror, doi:10.26131/IRSA674)"],
        "citation": "Euclid Collaboration, Euclid DR1 (2026)",
        "notes": "~1900 deg^2, >900M objects, ~160M spectra. NOT in DR1-Foundation: shear 2-pt functions, "
                 "clustering measurements, covariances (expected with 'DR1 Complete', mid-2027).",
        "tools": ["fetch_euclid_catalog"],
    },
    "euclid_eclipse": {
        "title": "ECLIPSE: Euclid DR1 strong-lens catalogue (26 662 expert-graded systems)",
        "products": {
            "catalog": ["lens candidates with expert grade A/B/C/X, scores, DESI deflector z"],
        },
        "access": "download (Zenodo, 3.7 MB object catalogue)",
        "urls": ["https://zenodo.org/records/19654904"],
        "citation": "Lines et al. (Euclid Collaboration), doi:10.5281/zenodo.19654904 (2026)",
        "tools": ["fetch_euclid_strong_lenses"],
    },
    # ---- DESI DR1 LSS clustering catalogs (data + randoms) ----
    "desi_dr1_lss": {
        "title": "DESI DR1 LSS clustering catalogs (LRG, ELG, QSO, NGC)",
        "products": {
            "catalog": ["LRG/ELG/QSO NGC clustering catalogs + randoms (v1.5)",
                        "official n(z)"],
            "summary_statistic": ["n(z) histogram", "xi(s) multipoles (user-measured)"],
        },
        "access": "download (HTTP Range, row-limited subsets)",
        "urls": ["https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/"],
        "citation": "DESI Collaboration, arXiv:2306.06308 (2023)",
        "tools": ["fetch_galaxy_lss_catalog", "measure_correlation_multipoles",
                  "compute_nz_histogram", "plot_clustering_multipoles"],
    },
    # ---- SDSS SkyServer (galaxy catalogs) ----
    "sdss_skyserver": {
        "title": "SDSS SkyServer SQL (DR18)",
        "products": {
            "catalog": ["spectroscopy", "photometry", "galaxy catalogs"],
        },
        "access": "remote-query",
        "urls": ["https://skyserver.sdss.org/dr18/"],
        "citation": "Abdurro'uf et al., ApJS 267, 44 (2023)",
        "tools": ["sdss_sql_query", "fetch_sdss_galaxies",
                  "compute_sdss_absolute_magnitudes",
                  "plot_color_magnitude", "plot_sdss_redshift_distribution"],
    },
    # ---- GSWLC (galaxy property catalog) ----
    "gswlc_a2": {
        "title": "GALEX-SDSS-WISE Legacy Catalog A2 (GSWLC-A2)",
        "products": {
            "catalog": ["stellar mass", "SFR", "dust attenuation", "640k galaxies"],
        },
        "access": "download",
        "urls": ["https://salims.pages.iu.edu/gswlc/GSWLC-A2.dat.gz"],
        "citation": "Salim et al., ApJS 227, 2 (2016); Salim, Boquien & Lee, ApJ 859, 11 (2018)",
        "tools": ["query_gswlc_catalog", "describe_gswlc_catalog",
                  "plot_star_forming_main_sequence", "crossmatch_sdss_gswlc"],
    },
}
