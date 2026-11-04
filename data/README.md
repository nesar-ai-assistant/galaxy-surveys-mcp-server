# Bundled data

Small canonical tables from the CobayaSampler/bao_data repository
(commit bb0c1c9, 2025-06-26). Total < 500 KB.

## Source

- **Repository**: https://github.com/CobayaSampler/bao_data
- **DESI DR2**: arXiv:2503.14738, arXiv:2503.14739
- **DESI DR1**: arXiv:2404.03000, arXiv:2404.03001, arXiv:2404.03002
- **eBOSS DR16**: arXiv:2007.08991
- **BOSS DR12**: arXiv:1607.03155
- **SDSS DR7 MGS**: arXiv:1409.3242
- **6dFGS**: arXiv:1106.3366 (hardcoded, not bundled)

## Files

### DESI DR2
- `desi_dr2_ALL_GCcomb_mean.txt` — 7 tracers, DM/rd DH/rd DV/rd
- `desi_dr2_ALL_GCcomb_cov.txt` — 13×13 covariance

### DESI DR1
- `desi_2024_gaussian_bao_ALL_GCcomb_mean.txt`
- `desi_2024_gaussian_bao_ALL_GCcomb_cov.txt`

### eBOSS DR16
- `sdss_DR16_LRG_BAO_DMDH.dat` + `_covtot.txt` — LRG BAO-only
- `sdss_DR16_QSO_BAO_DMDH.txt` + `_covtot.txt` — QSO BAO-only
- `sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8.dat` + `_covtot.txt` — LRG BAO+FS (includes f*sigma8)
- `sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8.dat` + `_covtot.txt` — QSO BAO+FS
- `sdss_DR16_ELG_BAO_DVtable.txt` — ELG DV/rd 1D grid likelihood
- `sdss_DR16_LYAUTO_BAO_DMDHgrid.txt` — Lya auto 2D grid
- `sdss_DR16_LYxQSO_BAO_DMDHgrid.txt` — Lya×QSO cross 2D grid

### BOSS DR12
- `sdss_DR12Consensus_bao.dat` + `BAO_consensus_covtot_dM_Hz.txt` — BAO-only (DM, H*rs)
- `sdss_DR12Consensus_final.dat` + `final_consensus_covtot_dM_Hz_fsig.txt` — BAO+FS (DM, H*rs, f*sigma8)
  Note: rs_fid = 147.78 Mpc for DR12 consensus

### SDSS DR7 MGS
- `sdss_MGS_prob.txt` — 1D probability distribution for DV/rs at z=0.15
