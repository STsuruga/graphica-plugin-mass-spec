import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mass_spec.chemistry import (  # noqa: E402
    ELECTRON_MASS, FormulaError, average_mass, format_formula, gaussian_profile, ion_pattern,
    isotope_pattern, monoisotopic_mass, parse_adduct, parse_formula, split_adducts,
)

H = 1.00782503223
NA = 22.9897692820
GLUCOSE_MONO = 180.06338810


@pytest.mark.parametrize("text, expected", [
    ("C6H12O6", {"C": 6, "H": 12, "O": 6}),
    ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}),
    ("(CH3)3COH", {"C": 4, "H": 10, "O": 1}),
    ("C6D6", {"C": 6, "D": 6}),
    ("[Cu(NH3)4]SO4", {"Cu": 1, "N": 4, "H": 12, "S": 1, "O": 4}),
    (" C 2 H 6 O ", {"C": 2, "H": 6, "O": 1}),
    ("CH3CH2OH", {"C": 2, "H": 6, "O": 1}),
])
def test_parse_formula(text, expected):
    assert parse_formula(text) == expected


@pytest.mark.parametrize("text", ["", "Xx2", "C6H12O6)", "(C6H12", "2H2O", "Tc", "c6h6"])
def test_parse_formula_rejects(text):
    with pytest.raises(FormulaError):
        parse_formula(text)


def test_format_formula_uses_hill_order():
    assert format_formula({"O": 6, "H": 12, "C": 6}) == "C6H12O6"
    assert format_formula({"Na": 1, "Cl": 1}) == "ClNa"
    assert format_formula({"N": 4, "C": 8, "O": 2, "H": 10}) == "C8H10N4O2"


def test_masses_of_glucose_and_caffeine():
    assert monoisotopic_mass(parse_formula("C6H12O6")) == pytest.approx(GLUCOSE_MONO, abs=1e-6)
    assert average_mass(parse_formula("C6H12O6")) == pytest.approx(180.156, abs=0.002)
    assert monoisotopic_mass(parse_formula("C8H10N4O2")) == pytest.approx(194.08037557, abs=1e-6)


def test_monoisotopic_mass_uses_the_most_abundant_isotope():
    # Fe は最も軽い 54Fe ではなく 56Fe
    assert monoisotopic_mass({"Fe": 1}) == pytest.approx(55.93493633, abs=1e-6)


@pytest.mark.parametrize("text, multimer, delta, charge", [
    ("[M+H]+", 1, {"H": 1}, 1),
    ("[M+Na]+", 1, {"Na": 1}, 1),
    ("[M+NH4]+", 1, {"N": 1, "H": 4}, 1),
    ("[M]+•", 1, {}, 1),
    ("M+.", 1, {}, 1),
    ("[M-H]-", 1, {"H": -1}, -1),
    ("[M−H]−", 1, {"H": -1}, -1),
    ("[M+HCOO]-", 1, {"H": 1, "C": 1, "O": 2}, -1),
    ("[M-H-H2O]-", 1, {"H": -3, "O": -1}, -1),
    ("[M+2H]2+", 1, {"H": 2}, 2),
    ("[M-2H]2-", 1, {"H": -2}, -2),
    ("[2M+Na]+", 2, {"Na": 1}, 1),
])
def test_parse_adduct(text, multimer, delta, charge):
    adduct = parse_adduct(text)
    assert (adduct.multimer, adduct.delta, adduct.charge) == (multimer, delta, charge)


@pytest.mark.parametrize("text", ["M+H", "[M+H]", "[M+Xx]+", "H+", "[M+H]0+"])
def test_parse_adduct_rejects(text):
    with pytest.raises(FormulaError):
        parse_adduct(text)


def test_split_adducts():
    assert split_adducts("[M+H]+, [M+Na]+ [M+NH4]+;[2M+Na]+") == ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[2M+Na]+"]


@pytest.mark.parametrize("adduct, expected", [
    ("[M+H]+", GLUCOSE_MONO + H - ELECTRON_MASS),
    ("[M+Na]+", GLUCOSE_MONO + NA - ELECTRON_MASS),
    ("[M-H]-", GLUCOSE_MONO - H + ELECTRON_MASS),
    ("[M]+•", GLUCOSE_MONO - ELECTRON_MASS),
    ("[M+2H]2+", (GLUCOSE_MONO + 2 * H - 2 * ELECTRON_MASS) / 2),
    ("[2M+Na]+", 2 * GLUCOSE_MONO + NA - ELECTRON_MASS),
])
def test_ion_monoisotopic_mz(adduct, expected):
    pattern = ion_pattern(parse_formula("C6H12O6"), parse_adduct(adduct))
    assert pattern.monoisotopic_mz == pytest.approx(expected, abs=1e-6)
    # 軽い元素だけなら最も軽いピークがモノアイソトピック
    assert pattern.mz[0] == pytest.approx(expected, abs=1e-6)


def test_ion_pattern_rejects_removing_missing_atoms():
    with pytest.raises(FormulaError):
        ion_pattern(parse_formula("CH4"), parse_adduct("[M-H2O]+"))


def test_chlorine_and_bromine_patterns():
    mz, rel, offsets = isotope_pattern({"Cl": 2})
    assert list(offsets) == [0, 2, 4]
    np.testing.assert_allclose(rel, [100, 2 * 0.2424 / 0.7576 * 100, (0.2424 / 0.7576) ** 2 * 100], rtol=1e-3)
    _mz, rel, _ = isotope_pattern({"Br": 1})
    # 79Br の方が多いので M が 100
    np.testing.assert_allclose(rel, [100, 100 * 0.4931 / 0.5069], rtol=1e-3)


def test_carbon_m_plus_1():
    _mz, rel, offsets = isotope_pattern({"C": 60})
    assert offsets[1] == 1
    ratio = 0.0107 / 0.9893
    assert rel[1] / rel[0] == pytest.approx(60 * ratio, rel=1e-6)
    assert rel[2] / rel[0] == pytest.approx(60 * 59 / 2 * ratio ** 2, rel=1e-6)


def test_bundled_peak_mass_is_the_abundance_weighted_mean():
    # C1H4 の M+1 は 13C と 2H の寄与の加重平均
    mz, rel, offsets = isotope_pattern({"C": 1, "H": 4})
    p13c = 0.0107 * 0.999885 ** 4
    p2h = 0.9893 * 4 * 0.000115 * 0.999885 ** 3
    m13c = 13.00335483507 + 4 * H
    m2h = 12.0 + 3 * H + 2.01410177812
    assert mz[1] == pytest.approx((p13c * m13c + p2h * m2h) / (p13c + p2h), abs=1e-9)


def test_pattern_for_a_large_molecule_is_fast_and_normalized():
    mz, rel, offsets = isotope_pattern(parse_formula("C200H300N50O60S5"))
    assert rel.max() == pytest.approx(100.0)
    assert np.all(np.diff(mz) > 0.9)
    assert len(mz) < 40


def test_gaussian_profile_width_matches_the_given_fwhm():
    x, y = gaussian_profile([500.0], [100.0], fwhm=0.05)
    assert x[np.argmax(y)] == pytest.approx(500.0, abs=0.05 / 20)
    above = x[y >= 50.0]
    assert above[-1] - above[0] == pytest.approx(0.05, rel=0.1)
