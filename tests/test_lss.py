"""Tests for the LSS clustering tool family.

Offline tests that don't require network access: unit-test the kernel,
verify tool registration, and mock-test the fetch/measure/plot pipeline.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


# ---------- kernel unit test (no network) -------------------------------------

def test_measure_xi_multipoles_synthetic():
    """Run measure_xi_multipoles on a small synthetic random catalog.

    Unclustered randoms should give xi_0 ~ 0 everywhere (within noise).
    """
    from tools.lss.kernels import measure_xi_multipoles

    rng = np.random.default_rng(42)
    n_data = 2000
    n_rand = 4000

    # Generate uniform points in a small patch
    ra = rng.uniform(180, 185, n_data)
    dec = rng.uniform(20, 25, n_data)
    z = rng.uniform(0.4, 0.7, n_data)
    w = np.ones(n_data)

    rra = rng.uniform(180, 185, n_rand)
    rdec = rng.uniform(20, 25, n_rand)
    rz = rng.uniform(0.4, 0.7, n_rand)
    rw = np.ones(n_rand)

    result = measure_xi_multipoles(
        ra.tolist(), dec.tolist(), z.tolist(), w.tolist(),
        rra.tolist(), rdec.tolist(), rz.tolist(), rw.tolist(),
        s_bins_min=10, s_bins_max=150, n_s_bins=14,
    )

    assert "s_mid" in result
    assert "xi0" in result
    assert "xi2" in result
    assert len(result["s_mid"]) == 14
    assert len(result["xi0"]) == 14
    assert result["n_data"] == n_data
    assert result["n_rand"] == n_rand

    # Unclustered data: xi0 should be near zero (within ~0.5 for small samples)
    xi0 = np.array(result["xi0"])
    # At least most bins should be within +/- 1.0
    assert np.median(np.abs(xi0)) < 1.0, f"xi0 too large for unclustered data: {xi0}"


# ---------- catalog registry --------------------------------------------------

def test_lss_catalog_registry():
    """All catalog keys have required fields."""
    from tools.lss.io import LSS_CATALOGS

    required = {"data_url", "random_url", "gzipped", "weight_cols", "z_col", "survey", "citation"}
    for key, cat in LSS_CATALOGS.items():
        for field in required:
            assert field in cat, f"Catalog {key} missing field {field}"
        assert cat["z_col"] == "Z"
        assert isinstance(cat["weight_cols"], list)
        assert len(cat["weight_cols"]) >= 2


def test_published_multipoles_registry():
    """All published multipole keys have required fields."""
    from tools.lss.io import PUBLISHED_MULTIPOLES

    for key, entry in PUBLISHED_MULTIPOLES.items():
        assert "type" in entry, f"{key} missing type"
        assert entry["type"] in ("pk", "xi"), f"{key} type={entry['type']}"
        assert "survey" in entry
        assert "citation" in entry
        assert "data_urls" in entry


# ---------- tool __all__ and data-kind decorators -----------------------------

def test_lss_tools_in_all():
    """All 5 LSS tools are listed in __all__."""
    from tools.lss import __all__
    expected = {
        "fetch_galaxy_lss_catalog",
        "measure_correlation_multipoles",
        "compute_nz_histogram",
        "fetch_published_multipoles",
        "plot_clustering_multipoles",
    }
    assert expected == set(__all__)


def test_lss_tools_have_data_kind():
    """Every LSS tool has a data_kind attribute from the @data_kind() decorator."""
    import tools.lss as lss
    for name in lss.__all__:
        fn = getattr(lss, name)
        assert hasattr(fn, "data_kind"), f"{name} missing @data_kind() decorator"
        assert fn.data_kind in ("catalog", "summary_statistic", "published_result"), \
            f"{name} has invalid data kind: {fn.data_kind}"


def test_lss_dispatch_kernels():
    """DISPATCH_KERNELS is declared for the heavy kernel."""
    from tools.lss import DISPATCH_KERNELS
    assert "measure_xi_multipoles" in DISPATCH_KERNELS
    k = DISPATCH_KERNELS["measure_xi_multipoles"]
    assert "function" in k
    assert "pip_deps" in k


# ---------- tool function contracts (mocked) ----------------------------------

def test_compute_nz_histogram_from_csv():
    """compute_nz_histogram on a small synthetic CSV."""
    from tools.lss import compute_nz_histogram

    rng = np.random.default_rng(123)
    n = 500
    df = pd.DataFrame({
        "RA": rng.uniform(0, 360, n),
        "DEC": rng.uniform(-30, 30, n),
        "Z": rng.uniform(0.4, 0.8, n),
        "WEIGHT_TOTAL": np.ones(n),
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = Path(tmpdir) / "test_data.csv"
        df.to_csv(csv_path, index=False)

        result = compute_nz_histogram(
            output_dir=tmpdir,
            catalog_file=str(csv_path),
            n_bins=20,
        )

    assert result.status == "success"
    assert len(result.files) >= 1
    assert result.metadata["n_bins"] == 20
    assert result.metadata["total_galaxies"] == n


def test_plot_clustering_multipoles_xi():
    """plot_clustering_multipoles on synthetic xi data."""
    from tools.lss import plot_clustering_multipoles

    s = np.linspace(10, 180, 30)
    xi0 = 0.1 * (100 / s) ** 1.5
    xi2 = -0.05 * (100 / s) ** 1.2

    with tempfile.TemporaryDirectory() as tmpdir:
        csv = Path(tmpdir) / "xi_test.csv"
        pd.DataFrame({"s_mid": s, "xi0": xi0, "xi2": xi2}).to_csv(csv, index=False)

        result = plot_clustering_multipoles(
            output_dir=tmpdir,
            data_file=str(csv),
            label="synthetic",
            plot_type="xi",
        )

    assert result.status == "success"
    assert any(f.endswith(".png") for f in result.files)


def test_plot_clustering_multipoles_pk():
    """plot_clustering_multipoles on synthetic P(k) data."""
    from tools.lss import plot_clustering_multipoles

    k = np.linspace(0.01, 0.3, 40)
    P0 = 5000 * np.exp(-k / 0.1)
    P2 = -1000 * np.exp(-k / 0.15)

    with tempfile.TemporaryDirectory() as tmpdir:
        csv = Path(tmpdir) / "pk_test.csv"
        pd.DataFrame({"k": k, "P0": P0, "P2": P2}).to_csv(csv, index=False)

        result = plot_clustering_multipoles(
            output_dir=tmpdir,
            data_file=str(csv),
            label="synthetic",
            plot_type="pk",
        )

    assert result.status == "success"
    assert any(f.endswith(".png") for f in result.files)
