"""Tests for SDSS and GSWLC tool families."""
import pytest
from pathlib import Path

from tools.common import ArtifactResult


# ---- SDSS tests (network) ----

@pytest.mark.network
class TestSDSSNetwork:
    def test_sdss_sql_query(self, tmp_path):
        from tools.sdss import sdss_sql_query
        r = sdss_sql_query(str(tmp_path), sql="SELECT TOP 5 objID, ra, dec, petroMag_r FROM PhotoObj WHERE petroMag_r < 17")
        assert isinstance(r, ArtifactResult)
        assert r.metadata.get("n_rows", 0) > 0

    def test_fetch_sdss_galaxies(self, tmp_path):
        from tools.sdss import fetch_sdss_galaxies
        r = fetch_sdss_galaxies(str(tmp_path), limit=10, min_z=0.02, max_z=0.05)
        assert isinstance(r, ArtifactResult)
        assert r.metadata.get("n_rows", 0) > 0

    def test_compute_absolute_magnitudes(self, tmp_path):
        from tools.sdss import fetch_sdss_galaxies, compute_sdss_absolute_magnitudes
        r1 = fetch_sdss_galaxies(str(tmp_path), limit=10, min_z=0.02, max_z=0.05)
        if r1.metadata.get("n_rows", 0) == 0:
            pytest.skip("No SDSS data available")
        r2 = compute_sdss_absolute_magnitudes(str(tmp_path), input_file=r1.files[0])
        assert isinstance(r2, ArtifactResult)
        assert len(r2.files) == 1

    def test_plot_color_magnitude(self, tmp_path):
        from tools.sdss import fetch_sdss_galaxies, compute_sdss_absolute_magnitudes, plot_color_magnitude
        r1 = fetch_sdss_galaxies(str(tmp_path), limit=50, min_z=0.02, max_z=0.05)
        if r1.metadata.get("n_rows", 0) == 0:
            pytest.skip("No SDSS data available")
        r2 = compute_sdss_absolute_magnitudes(str(tmp_path), input_file=r1.files[0])
        r3 = plot_color_magnitude(str(tmp_path), input_file=r2.files[0])
        assert isinstance(r3, ArtifactResult)
        assert Path(r3.files[0]).exists()

    def test_plot_redshift_distribution(self, tmp_path):
        from tools.sdss import fetch_sdss_galaxies, plot_sdss_redshift_distribution
        r1 = fetch_sdss_galaxies(str(tmp_path), limit=50, min_z=0.02, max_z=0.05)
        if r1.metadata.get("n_rows", 0) == 0:
            pytest.skip("No SDSS data available")
        r2 = plot_sdss_redshift_distribution(str(tmp_path), input_file=r1.files[0])
        assert isinstance(r2, ArtifactResult)
        assert Path(r2.files[0]).exists()


# ---- GSWLC tests (network for download, then local) ----

@pytest.mark.network
class TestGSWLCNetwork:
    def test_describe_gswlc(self):
        from tools.gswlc import describe_gswlc_catalog
        r = describe_gswlc_catalog()
        assert isinstance(r, ArtifactResult)
        assert "logMstar" in r.metadata["columns"]

    def test_query_gswlc(self, tmp_path):
        from tools.gswlc import query_gswlc_catalog
        r = query_gswlc_catalog(
            str(tmp_path),
            filters='{"logMstar":{"min":10,"max":11},"z":{"max":0.1}}',
            limit=100,
        )
        assert isinstance(r, ArtifactResult)
        assert r.metadata["n_rows"] > 0

    def test_plot_sfms(self, tmp_path):
        from tools.gswlc import query_gswlc_catalog, plot_star_forming_main_sequence
        r1 = query_gswlc_catalog(
            str(tmp_path),
            filters='{"logMstar":{"min":9},"flag_mgs":{"eq":1}}',
            limit=5000,
        )
        r2 = plot_star_forming_main_sequence(str(tmp_path), input_file=r1.files[0])
        assert isinstance(r2, ArtifactResult)
        assert Path(r2.files[0]).exists()
