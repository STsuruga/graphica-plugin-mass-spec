import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mass_spec.chemistry import FormulaError, parse_adduct, parse_formula  # noqa: E402
from mass_spec.matching import TABLE_COLUMNS, match_adduct, match_all, results_table  # noqa: E402
from synthetic import ion_spectrum, tof_grid  # noqa: E402

GLUCOSE = parse_formula("C6H12O6")


def _spectrum(shift_ppm=3.0, resolution=12000, noise=0.0, formula="C6H12O6", adduct="[M+Na]+"):
    grid = tof_grid(150.0, 450.0)
    y, pattern = ion_spectrum(grid, formula, adduct, height=1000.0, resolution=resolution, shift_ppm=shift_ppm)
    if noise:
        y = y + np.random.default_rng(0).normal(0, noise, len(y)).clip(0)
    return grid, y, pattern


def test_found_adduct_reports_ppm_relative_intensity_and_resolution():
    mz, y, pattern = _spectrum(shift_ppm=3.0, noise=1.0)
    result = match_adduct(GLUCOSE, parse_adduct("[M+Na]+"), mz, y)
    assert result.found and result.status == "検出"
    mono = result.monoisotopic
    assert mono.error_ppm == pytest.approx(3.0, abs=1.0)
    m1 = next(p for p in result.peaks if p.label == "M+1")
    assert m1.measured_relative == pytest.approx(m1.calc_relative, abs=1.0)
    assert result.resolution_source == "実測"
    assert result.resolution == pytest.approx(12000, rel=0.15)
    assert result.scale == pytest.approx(1000.0, rel=0.03)


def test_fixed_resolution_is_used_as_given():
    mz, y, _ = _spectrum()
    result = match_adduct(GLUCOSE, parse_adduct("[M+Na]+"), mz, y, resolution=8000)
    assert (result.resolution, result.resolution_source) == (8000, "指定")


def test_absent_and_out_of_range_adducts():
    mz, y, _ = _spectrum()
    absent = match_adduct(GLUCOSE, parse_adduct("[M+K]+"), mz, y)
    assert (absent.found, absent.status) == (False, "未検出")
    assert np.isnan(absent.peaks[0].measured_mz)
    out = match_adduct(GLUCOSE, parse_adduct("[3M+Na]+"), mz, y)  # m/z 563 は 150〜450 の外
    assert out.status == "測定範囲外"


def test_weak_peak_below_detection_limit_is_not_found():
    mz, y, _ = _spectrum()
    big, _ = ion_spectrum(mz, "C8H10N4O2", "[M+H]+", height=1e6)
    result = match_adduct(GLUCOSE, parse_adduct("[M+Na]+"), mz, y + big, min_detect_percent=0.5)
    assert not result.found


def test_profile_and_sticks_are_scaled_to_the_measurement():
    mz, y, _ = _spectrum()
    result = match_adduct(GLUCOSE, parse_adduct("[M+Na]+"), mz, y)
    px, py = result.profile()
    sx, sy = result.sticks()
    assert py.max() == pytest.approx(sy.max(), rel=0.01)
    assert sy.max() == pytest.approx(y.max(), rel=0.03)
    assert px.min() < sx.min() < sx.max() < px.max()


def test_results_table_includes_errors_and_all_rows():
    mz, y, _ = _spectrum()
    results = match_all(GLUCOSE, ["[M+Na]+", "[M+K]+", "[M-H2O+X]+"], mz, y)
    assert isinstance(results[2], FormulaError)
    table = results_table(results, ["[M+Na]+", "[M+K]+", "[M-H2O+X]+"])
    assert list(table.columns) == TABLE_COLUMNS
    assert set(table["状態"].iloc[:-1]) == {"検出", "未検出"}
    assert table["状態"].iloc[-1].startswith("エラー")
    assert table.loc[0, "ピーク"] == "M"


def test_iron_pattern_labels_m_minus_2():
    grid = tof_grid(150.0, 450.0)
    y, _ = ion_spectrum(grid, "C10H10Fe", "[M]+•", resolution=12000)
    result = match_adduct(parse_formula("C10H10Fe"), parse_adduct("[M]+•"), grid, y)
    labels = [p.label for p in result.peaks]
    assert labels[:2] == ["M-2", "M"]
    assert result.found
