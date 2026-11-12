"""Euclid DR1 tools: offline ADQL/region unit tests + live archive / Zenodo tests (network)."""
import json

import pytest

from tools.euclid.io import TABLES, build_query, region_clause, vis_flux_limit


def test_region_cone_and_box():
    c, a = region_clause(80.3, -51.25, 0.3, None)
    assert "CIRCLE('ICRS', 80.3, -51.25, 0.3)" in c and a == pytest.approx(0.2827, rel=1e-3)
    b, a = region_clause(None, None, None, [80.2, 80.4, -51.35, -51.15])
    assert "BETWEEN 80.2 AND 80.4" in b and 0.02 < a < 0.03
    w, _ = region_clause(None, None, None, [359.0, 1.0, -1, 1])   # wraps RA = 0
    assert "OR" in w
    for bad in [dict(ra=1, dec=1, radius_deg=5), dict(ra=None, dec=1, radius_deg=1)]:
        with pytest.raises(ValueError):
            region_clause(bad["ra"], bad["dec"], bad["radius_deg"], None)
    with pytest.raises(ValueError):
        region_clause(None, None, None, [1, 2, 3])


def test_build_queries():
    reg, _ = region_clause(10, 10, 0.1, None)
    q = build_query("spectroscopic", reg, max_rows=100, z_min=0.9, z_max=1.8, mag_vis_max=24.0)
    assert q.startswith("SELECT TOP 100 ")
    assert TABLES["spe_z"] in q and "s.spe_rank = 0" in q and "s.spe_z_prob >= 0.99" in q
    assert "c.spe_class = 'galaxy'" in q and "s.spe_z >= 0.9" in q
    assert "LOG10(" not in q.upper() and "LOG(" not in q.upper(), "no server-side log(): the archive aborts on negative fluxes"
    assert f"m.flux_vis_sersic > {vis_flux_limit(24.0):.6g}" in q
    q = build_query("photometric", reg, max_rows=5, sample_mod=10, require_phz_flags_zero=True)
    assert TABLES["phz"] in q and "p.phz_classification = 2" in q and "p.phz_flags = 0" in q
    assert "MOD(ABS(m.object_id), 10) = 0" in q
    assert vis_flux_limit(23.9) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        build_query("weak_lensing", reg, max_rows=1)


@pytest.mark.network
def test_live_euclid_spectroscopic_sample(tmp_path):
    from tools.euclid import fetch_euclid_catalog
    r = fetch_euclid_catalog(str(tmp_path), sample="spectroscopic", ra=80.3, dec=-51.25, radius_deg=0.05,
                             max_rows=2000)
    m = r.metadata
    assert m["data_kind"] == "catalog" and m["n_rows"] > 20
    assert 0.3 < m["z"]["median"] < 2.5
    assert m["mag_vis"]["max"] < 24.5


@pytest.mark.network
def test_live_euclid_photometric_box(tmp_path):
    from tools.euclid import fetch_euclid_catalog
    r = fetch_euclid_catalog(str(tmp_path), sample="photometric", bounds=[80.25, 80.3, -51.27, -51.23],
                             max_rows=2000, plot=False)
    assert r.metadata["n_rows"] > 50 and 0 < r.metadata["z"]["median"] < 3


@pytest.mark.network
def test_eclipse_strong_lenses(tmp_path):
    from tools.euclid import fetch_euclid_strong_lenses
    r = fetch_euclid_strong_lenses(str(tmp_path), grades=["A", "B"])
    m = r.metadata
    assert m["n_catalogue"] == 26662
    assert m["grade_counts_catalogue"] == {"A": 6073, "B": 8999, "C": 9010, "X": 2580}
    assert m["n_selected"] == 6073 + 8999
    s = json.loads(open(r.files[1]).read())
    assert s["source_file"] == "ECLIPSE_object_catalogue.zip"
    r = fetch_euclid_strong_lenses(str(tmp_path), grades=["A"], min_expert_score=2.9, plot=False)
    assert 0 < r.metadata["n_selected"] < 6073
