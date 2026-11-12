"""DES Y6 angular BAO: likelihood profiles (bundled), joint fits, released data vectors (network)."""
import json

import numpy as np
import pytest

from tools.bao import fetch_bao_measurements, fetch_des_y6_bao_likelihood, fetch_des_y6_clustering, fit_bao_cosmology
from tools.bao.des_y6 import DM_RD_FID, VARIANTS, load_profile
from tools.bao.kernels import _profile_chi2, fit_bao


def test_reproduces_published_des_y6_bao():
    """Consensus profile -> alpha = 0.9571 +/- 0.0201, D_M/r_d = 19.51 +/- 0.41 (arXiv:2402.10696)."""
    p = load_profile("fid")
    assert p["alpha_best"] == pytest.approx(0.9571, abs=2e-4)
    assert p["alpha_sigma"] == pytest.approx(0.0201, abs=5e-4)
    assert p["DM_over_rd"] == pytest.approx(19.51, abs=0.01)
    assert p["DM_over_rd_sigma"] == pytest.approx(0.41, abs=0.01)
    assert DM_RD_FID == pytest.approx(20.38, abs=0.01)
    assert p["z"] == 0.851


def test_all_variants_load_and_are_sensible():
    for v in VARIANTS:
        p = load_profile(v)
        assert p["chi2"].min() == 0.0
        assert 0.93 < p["alpha_best"] < 1.0, v
        assert p["alpha_lo"] < p["alpha_best"] < p["alpha_hi"]
    # removing the DESI area loses ~half the area -> wider constraint
    assert load_profile("nodesi_dr1tiles")["alpha_sigma"] > 1.3 * load_profile("fid")["alpha_sigma"]


def test_profile_chi2_matches_table():
    p = load_profile("fid")
    prof = {"z": p["z"], "observable": "DM_over_rs", "fid": p["fid"],
            "alpha": p["alpha"].tolist(), "chi2": p["chi2"].tolist()}
    from tools.bao.kernels import _theory_vec
    # Planck-like point: alpha ~ 1 -> Delta chi2 ~ (1-0.957)/0.020)^2 ~ 4.5
    th = _theory_vec([0.851], ["DM_over_rs"], 0.3153, 0.6736 * 147.8)[0]
    c2 = _profile_chi2([prof], 0.3153, 0.6736 * 147.8)
    assert c2 == pytest.approx(np.interp(th / p["fid"], p["alpha"], p["chi2"]))
    assert 3 < c2 < 6


def test_des_only_fit_is_degenerate_but_runs():
    p = load_profile("fid")
    prof = {"z": p["z"], "observable": "DM_over_rs", "fid": p["fid"],
            "alpha": p["alpha"].tolist(), "chi2": p["chi2"].tolist()}
    res = fit_bao([], [], [], [], profiles=[prof], n_grid_Om=20, n_grid_H0rd=20)
    assert res["chi2_min"] < 0.5 and res["ndof"] == -1


def test_joint_desi_dr2_des_y6_fit(tmp_path):
    r = fit_bao_cosmology(str(tmp_path), datasets=["desi_dr2", "des_y6_nodesi"], n_grid=40)
    m = r.metadata
    assert m["datasets"] == ["desi_dr2", "des_y6_nodesi"]
    assert m["n_gaussian_points"] == 13 and m["ndof"] == 12
    assert 0.28 < m["marginalised"]["Om"] < 0.31
    assert "warnings" not in m
    r2 = fit_bao_cosmology(str(tmp_path), dataset="desi_dr2,des_y6", n_grid=40)
    assert any("overlaps DESI" in w for w in r2.metadata["warnings"])
    # DES adds information: chi2_min grows relative to DESI alone
    r0 = fit_bao_cosmology(str(tmp_path), dataset="desi_dr2", n_grid=40)
    assert r2.metadata["chi2_min"] > r0.metadata["chi2_min"]
    with pytest.raises(ValueError):
        fit_bao_cosmology(str(tmp_path), datasets=["des_y6", "des_y6_nodesi"], n_grid=20)


def test_fetch_tools(tmp_path):
    r = fetch_des_y6_bao_likelihood(str(tmp_path), variant="fid", plot=True)
    assert len(r.files) == 3 and r.metadata["data_kind"] == "summary_statistic"
    s = json.loads(open(r.files[1]).read())
    assert s["DM_over_rd"] == pytest.approx(19.51, abs=0.01)
    r = fetch_bao_measurements(str(tmp_path), dataset="des_y6")
    assert r.metadata["n_points"] == 1
    with pytest.raises(ValueError):
        load_profile("nope")


@pytest.mark.network
@pytest.mark.parametrize("stat,n", [("acf", 144), ("aps", 210), ("pcf", 120)])
def test_des_y6_clustering_download(tmp_path, stat, n):
    r = fetch_des_y6_clustering(str(tmp_path), statistic=stat, plot=False)
    assert r.metadata["n_points"] == n
    cov = np.loadtxt(r.files[1], delimiter=",")
    assert cov.shape == (n, n)
    import pandas as pd
    df = pd.read_csv(r.files[0])
    assert np.allclose(np.sqrt(np.diag(cov)), df["error"], rtol=1e-4)
    assert sorted(df["zbin"].unique()) == [1, 2, 3, 4, 5, 6]
