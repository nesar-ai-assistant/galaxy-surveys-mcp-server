"""Euclid DR1 tool family: galaxy catalogs from the ESA Euclid Science Archive + ECLIPSE strong lenses.

Euclid DR1-Foundation (public 2026-11-12; ~1900 deg^2 Wide Survey) provides
object catalogues: MER photometry, photometric redshifts, NISP spectroscopic
redshifts and classifications. It does NOT contain cosmological summary
statistics (cosmic-shear / clustering 2-pt functions, covariances); those are
expected with 'DR1 Complete' (mid-2027). These tools therefore return
catalogs only.
"""
from __future__ import annotations

import json
from typing import Annotated, Literal

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, data_kind, mpl, ok, param_slug, resolve_outdir, save_fig

__all__ = ["fetch_euclid_catalog", "fetch_euclid_strong_lenses"]

NOT_IN_DR1 = ("Euclid DR1-Foundation has catalogues only; weak-lensing shear 2-pt functions, galaxy "
              "clustering measurements and covariances are not released (expected with 'DR1 Complete', mid-2027).")


def _zsummary(z) -> dict:
    z = np.asarray(z, dtype=float)
    z = z[np.isfinite(z)]
    if not len(z):
        return {"n": 0}
    q = np.percentile(z, [5, 16, 50, 84, 95])
    return {"n": int(len(z)), "min": float(z.min()), "max": float(z.max()), "mean": float(z.mean()),
            "p05": float(q[0]), "p16": float(q[1]), "median": float(q[2]), "p84": float(q[3]), "p95": float(q[4])}


@data_kind("catalog")
@validate_call
def fetch_euclid_catalog(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    sample: Annotated[Literal["photometric", "spectroscopic"],
                      Field(description="photometric = MER x photo-z (phz_median); "
                                        "spectroscopic = MER x NISP spectroscopic z (spe_rank=0)")] = "spectroscopic",
    ra: Annotated[float | None, Field(description="Cone centre RA [deg]")] = None,
    dec: Annotated[float | None, Field(description="Cone centre Dec [deg]")] = None,
    radius_deg: Annotated[float | None, Field(gt=0, le=3, description="Cone radius [deg] (<= 3)")] = None,
    bounds: Annotated[list[float] | None, Field(description="Box [ra_min, ra_max, dec_min, dec_max] in deg "
                                                           "(alternative to a cone; RA may wrap)")] = None,
    z_min: Annotated[float | None, Field(description="Minimum redshift (spe_z or phz_median)")] = None,
    z_max: Annotated[float | None, Field(description="Maximum redshift")] = None,
    mag_vis_max: Annotated[float | None, Field(description="VIS (Sersic) AB magnitude limit; null = none")] = 24.5,
    min_spec_z_prob: Annotated[float, Field(ge=0, le=1, description="spectroscopic only: minimum spe_z_prob "
                                                                    "(probability of the rank-0 redshift)")] = 0.99,
    require_galaxy_class: Annotated[bool, Field(description="Keep only objects classified as galaxies "
                                                            "(phz_classification=2 / spe_class='galaxy')")] = True,
    require_phz_flags_zero: Annotated[bool, Field(description="photometric only: require phz_flags = 0")] = False,
    max_rows: Annotated[int, Field(ge=1, le=200000, description="Server-side TOP limit")] = 20000,
    subsample_mod: Annotated[int, Field(ge=1, le=1000, description="Keep 1 in N objects (MOD(object_id, N)) "
                                                                   "server-side, a deterministic random-like "
                                                                   "subsample so TOP does not bias by position")] = 1,
    output_format: Annotated[Literal["csv", "parquet"], Field(description="Table format")] = "csv",
    plot: Annotated[bool, Field(description="Write an n(z) histogram")] = True,
) -> ArtifactResult:
    """Query a Euclid DR1 Wide-Survey galaxy sample in a sky region (ESA Euclid archive, anonymous TAP).

    Writes the table (CSV/parquet), a summary JSON (N, redshift distribution,
    the exact ADQL) and an optional n(z) plot. Tables used:
    dr1.mer_catalogue_wide_survey joined to dr1.phz_photo_z_wide (photometric)
    or dr1.spectro_z_spe_galaxy_candidates_wide_survey +
    dr1.spectro_z_spe_classification_wide_survey (spectroscopic).
    Quality cuts: spurious_flag = 0, VIS magnitude limit, galaxy class,
    spe_z_prob >= min_spec_z_prob (spectroscopic) or optional phz_flags = 0.
    A 0.3 deg cone holds ~25k spectroscopic candidates / ~5k with spe_z_prob > 0.99.
    """
    from .io import MAX_AREA_DEG2, TAP_URL, build_query, region_clause, run_tap
    region, area = region_clause(ra, dec, radius_deg, bounds)
    if area > MAX_AREA_DEG2:
        raise ValueError(f"Region area {area:.1f} deg^2 exceeds {MAX_AREA_DEG2} deg^2; split it into smaller queries")
    adql = build_query(sample, region, max_rows=max_rows, sample_mod=subsample_mod, z_min=z_min, z_max=z_max,
                       mag_vis_max=mag_vis_max, min_spec_z_prob=min_spec_z_prob,
                       require_galaxy_class=require_galaxy_class, require_phz_flags_zero=require_phz_flags_zero)
    df, backend = run_tap(adql)
    zcol = "spe_z" if sample == "spectroscopic" else "phz_median"
    out = resolve_outdir(output_dir)
    params = dict(sample=sample, ra=ra, dec=dec, r=radius_deg, b=bounds, z0=z_min, z1=z_max, m=mag_vis_max,
                  p=min_spec_z_prob, g=require_galaxy_class, f=require_phz_flags_zero, n=max_rows, s=subsample_mod)
    stem = f"euclid_dr1_{sample}_{param_slug(params)}"
    if output_format == "parquet":
        try:
            path = out / f"{stem}.parquet"
            df.to_parquet(path, index=False)
        except ImportError:
            path = out / f"{stem}.csv"
            df.to_csv(path, index=False)
            output_format = "csv (parquet engine not installed)"
    else:
        path = out / f"{stem}.csv"
        df.to_csv(path, index=False)
    truncated = len(df) >= max_rows
    summary = {
        "sample": sample, "n_rows": int(len(df)), "truncated_at_max_rows": truncated,
        "region_area_deg2": round(area, 4), "redshift_column": zcol,
        "z": _zsummary(df[zcol]) if zcol in df else {"n": 0},
        "mag_vis": _zsummary(df["mag_vis"]) if "mag_vis" in df else {},
        "adql": adql, "service": TAP_URL, "backend": backend, "release": "Euclid DR1 (DR1-Foundation)",
        "output_format": output_format, "columns": list(df.columns), "not_in_dr1": NOT_IN_DR1,
        "citation": "Euclid Collaboration, Euclid DR1 (2026); ESA Euclid Science Archive",
    }
    js = out / f"{stem}_summary.json"
    js.write_text(json.dumps(summary, indent=2, default=float))
    files = [path, js]
    if plot and len(df):
        plt = mpl()
        fig, ax = plt.subplots(figsize=(6.5, 4))
        z = df[zcol].to_numpy(dtype=float)
        z = z[np.isfinite(z)]
        ax.hist(z, bins=min(60, max(10, int(np.sqrt(len(z))))), histtype="stepfilled", alpha=0.6)
        ax.set_xlabel("spectroscopic z (spe_z)" if sample == "spectroscopic" else "photo-z (phz_median)")
        ax.set_ylabel("N")
        ax.set_title(f"Euclid DR1 {sample} sample: N = {len(df)}" + (" (TOP limit hit)" if truncated else ""))
        files.append(save_fig(fig, out, f"{stem}_nz"))
    zs = summary["z"]
    msg = f"Euclid DR1 {sample}: {len(df)} galaxies"
    if zs.get("n"):
        msg += f", z in [{zs['min']:.3f}, {zs['max']:.3f}], median {zs['median']:.3f}"
    if truncated:
        msg += f" (hit max_rows={max_rows}; increase it or use subsample_mod)"
    return ok(msg, files, **summary)


@data_kind("catalog")
@validate_call
def fetch_euclid_strong_lenses(
    output_dir: Annotated[str, Field(description="Directory for output files.")],
    grades: Annotated[list[Literal["A", "B", "C", "X"]] | None,
                      Field(description="Expert grades to keep: A sure, B probable, C possible, X not a lens "
                                        "(default ['A','B'])")] = None,
    min_expert_score: Annotated[float | None, Field(ge=0, le=3, description="Minimum Expert_Score (0-3)")] = None,
    min_ensemble_score: Annotated[float | None, Field(ge=0, le=1, description="Minimum neural-network Ensemble_score")] = None,
    require_desi_specz: Annotated[bool, Field(description="Keep only systems with a DESI DR1 deflector redshift")] = False,
    ra: Annotated[float | None, Field(description="Optional cone centre RA [deg]")] = None,
    dec: Annotated[float | None, Field(description="Optional cone centre Dec [deg]")] = None,
    radius_deg: Annotated[float | None, Field(gt=0, le=180, description="Optional cone radius [deg]")] = None,
    plot: Annotated[bool, Field(description="Sky map + score histogram")] = True,
) -> ArtifactResult:
    """ECLIPSE: the Euclid DR1 strong-lens candidate catalogue (Zenodo record 19654904).

    Downloads once (3.7 MB ECLIPSE_object_catalogue.zip; the record's multi-GB
    cutout archives are not fetched) and filters by expert grade/score,
    neural-network ensemble score, DESI deflector redshift and sky cone.
    Full catalogue: 26 662 expert-graded systems (6073 A, 8999 B, 9010 C, 2580 X).
    """
    from .io import ECLIPSE_CITATION, load_eclipse
    grades = list(grades) if grades else ["A", "B"]
    df, csv_path, info = load_eclipse()
    n_total = len(df)
    counts_all = {k: int(v) for k, v in df["Expert_Grade"].value_counts().sort_index().items()}
    m = df["Expert_Grade"].isin(grades)
    if min_expert_score is not None:
        m &= df["Expert_Score"] >= min_expert_score
    if min_ensemble_score is not None:
        m &= df["Ensemble_score"] >= min_ensemble_score
    if require_desi_specz:
        m &= df["desi_deflector_specz"].notna()
    if radius_deg is not None:
        if ra is None or dec is None:
            raise ValueError("cone filter needs ra, dec and radius_deg")
        r1, d1 = np.radians(df["RA_obj"].to_numpy()), np.radians(df["DEC_obj"].to_numpy())
        r0, d0 = np.radians(ra), np.radians(dec)
        cosd = np.sin(d1) * np.sin(d0) + np.cos(d1) * np.cos(d0) * np.cos(r1 - r0)
        m &= np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= radius_deg
    sel = df[m]
    out = resolve_outdir(output_dir)
    stem = "euclid_eclipse_" + param_slug(dict(g="".join(sorted(grades)), e=min_expert_score, n=min_ensemble_score,
                                               d=require_desi_specz, ra=ra, dec=dec, r=radius_deg))
    path = out / f"{stem}.csv"
    sel.to_csv(path, index=False)
    zd = sel["desi_deflector_specz"].dropna()
    summary = {
        "n_selected": int(len(sel)), "n_catalogue": int(n_total), "grade_counts_catalogue": counts_all,
        "grade_counts_selected": {k: int(v) for k, v in sel["Expert_Grade"].value_counts().sort_index().items()},
        "n_with_desi_deflector_z": int(len(zd)), "desi_deflector_z": _zsummary(zd),
        "n_previously_published": int(sel["earliest_publication"].notna().sum()),
        "n_large_system_warning": int(sel["large_system_warning"].sum()),
        "filters": {"grades": grades, "min_expert_score": min_expert_score,
                    "min_ensemble_score": min_ensemble_score, "require_desi_specz": require_desi_specz,
                    "cone": [ra, dec, radius_deg] if radius_deg else None},
        "source_file": info.get("key"), "zenodo_doi": info.get("doi"), "cached_csv": str(csv_path),
        "citation": ECLIPSE_CITATION,
    }
    js = out / f"{stem}_summary.json"
    js.write_text(json.dumps(summary, indent=2, default=float))
    files = [path, js]
    if plot and len(sel):
        plt = mpl()
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
        for g, c in zip("ABCX", ["C3", "C1", "C0", "gray"], strict=False):
            s = sel[sel["Expert_Grade"] == g]
            if len(s):
                a1.scatter(s["RA_obj"], s["DEC_obj"], s=2, c=c, label=f"{g} ({len(s)})")
                a2.hist(s["Expert_Score"], bins=30, range=(0, 3), color=c, alpha=0.6, label=g)
        a1.invert_xaxis()
        a1.set_xlabel("RA [deg]")
        a1.set_ylabel("Dec [deg]")
        a1.legend(markerscale=4, fontsize=7)
        a2.set_xlabel("Expert_Score")
        a2.set_ylabel("N")
        a2.legend(fontsize=7)
        fig.suptitle(f"ECLIPSE Euclid DR1 strong lenses: {len(sel)} selected")
        files.append(save_fig(fig, out, f"{stem}_overview"))
    return ok(f"ECLIPSE: {len(sel)} of {n_total} systems selected "
              f"(grades {','.join(grades)}; selected counts {summary['grade_counts_selected']}).", files, **summary)
