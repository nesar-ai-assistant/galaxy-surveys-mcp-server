"""Tests for BAO tool family."""
import numpy as np
from pathlib import Path

from tools.bao import (
    fetch_bao_measurements, plot_bao_distances, plot_fsigma8, compare_bao_surveys,
)
from tools.bao.io import load_gaussian, SURVEY_MAP
from tools.bao.kernels import fit_bao, _DM, _DH, _DV, _Ez
from tools.common import ArtifactResult


def test_load_all_bao_datasets():
    """Every BAO dataset key loads without error."""
    for key in SURVEY_MAP:
        d = load_gaussian(key)
        assert "z" in d and "value" in d and "cov" in d
        assert len(d["z"]) == len(d["value"])
        assert d["cov"].shape[0] == d["cov"].shape[1] == len(d["z"])


def test_desi_dr2_shape():
    d = load_gaussian("desi_dr2")
    assert len(d["z"]) == 13  # 7 tracers, mix of DM/DH/DV


def test_desi_dr1_shape():
    d = load_gaussian("desi_dr1")
    assert len(d["z"]) == 12  # DR1 has 12 data points


def test_boss_dr12_conversion():
    """BOSS DR12 data should be in DM/rd and DH/rd units after conversion."""
    d = load_gaussian("boss_dr12")
    obs = d["observable"]
    assert "DM_over_rs" in obs
    assert "DH_over_rs" in obs
    # Values should be order 10-30 (DM/rd), not 1500+ (raw Mpc)
    for v, o in zip(d["value"], obs, strict=False):
        if o == "DM_over_rs":
            assert 5 < v < 25, f"DM/rd = {v} looks wrong"
        elif o == "DH_over_rs":
            assert 10 < v < 50, f"DH/rd = {v} looks wrong"


def test_cosmology_helpers():
    """E(z) and distances are physically sensible."""
    assert abs(_Ez(0, 0.3) - 1.0) < 1e-10
    assert _Ez(1.0, 0.3) > 1.0
    dm = _DM(0.5, 0.3, 70)
    assert 1000 < dm < 2000  # ~1300 Mpc
    dh = _DH(0.5, 0.3, 70)
    assert 2000 < dh < 4000
    dv = _DV(0.5, 0.3, 70)
    assert dv > 0


def test_fit_bao_lcdm():
    """Fit DESI DR2 and check Om is near 0.3."""
    d = load_gaussian("desi_dr2")
    bao_mask = [o in ("DM_over_rs", "DH_over_rs", "DV_over_rs") for o in d["observable"]]
    result = fit_bao(
        z_arr=d["z"][bao_mask].tolist(),
        obs_arr=[o for o, m in zip(d["observable"], bao_mask, strict=False) if m],
        data_vec=d["value"][bao_mask].tolist(),
        cov_mat=d["cov"][np.ix_(bao_mask, bao_mask)].tolist(),
        n_grid_Om=50, n_grid_H0rd=50,
    )
    assert result["model"] == "LCDM"
    Om = result["best_fit"]["Om"]
    assert 0.25 < Om < 0.35, f"Om = {Om} is outside expected range"
    H0rd = result["best_fit"]["H0rd"]
    assert 95 < H0rd < 108, f"H0rd = {H0rd} is outside expected range"


def test_fetch_bao_measurements(tmp_path):
    r = fetch_bao_measurements(str(tmp_path), dataset="desi_dr2")
    assert isinstance(r, ArtifactResult)
    assert r.metadata["n_points"] == 13
    assert len(r.files) == 2


def test_plot_bao_distances(tmp_path):
    r = plot_bao_distances(str(tmp_path))
    assert isinstance(r, ArtifactResult)
    assert len(r.files) == 1
    assert Path(r.files[0]).exists()


def test_plot_fsigma8(tmp_path):
    r = plot_fsigma8(str(tmp_path), datasets="boss_dr12")
    assert isinstance(r, ArtifactResult)
    assert Path(r.files[0]).exists()


def test_compare_bao_surveys(tmp_path):
    r = compare_bao_surveys(str(tmp_path), datasets="desi_dr2,boss_dr12_bao")
    assert isinstance(r, ArtifactResult)
    assert Path(r.files[0]).exists()
