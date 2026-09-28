"""計算した同位体パターンと実測スペクトルの照合。GUI にも Graphica 本体にも依存しない。"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .chemistry import FormulaError, gaussian_profile, ion_pattern, parse_adduct
from .spectra import find_peaks

DEFAULT_FWHM = 0.1   # 計算パターンの半値全幅(m/z)
TABLE_MIN_RELATIVE = 1.0  # 表に出す計算ピークの相対強度の下限(%)


@dataclass
class PeakMatch:
    label: str                  # M, M+1, ...
    calc_mz: float
    calc_relative: float
    measured_mz: float = float("nan")
    measured_height: float = float("nan")
    measured_relative: float = float("nan")

    @property
    def error_ppm(self):
        return (self.measured_mz - self.calc_mz) / self.calc_mz * 1e6


@dataclass
class AdductMatch:
    pattern: object             # chemistry.IonPattern
    found: bool
    status: str                 # 検出 / 未検出 / 測定範囲外
    fwhm: float                 # 計算パターンの半値全幅(m/z)
    scale: float                # 計算の相対強度 100 に対応する実測の強度
    peaks: list = field(default_factory=list)

    @property
    def base(self):
        """計算で最も強いピークの照合結果。"""
        return max(self.peaks, key=lambda p: p.calc_relative) if self.peaks else None

    @property
    def monoisotopic(self):
        return next((p for p in self.peaks if p.label == "M"), None)

    def profile(self, points_per_fwhm=20):
        """実測の強度に合わせたガウスのプロファイル (m/z, 強度)。"""
        x, y = gaussian_profile(self.pattern.mz, self.pattern.relative, self.fwhm, points_per_fwhm)
        return x, y * self.scale / 100.0

    def sticks(self):
        return self.pattern.mz, self.pattern.relative * self.scale / 100.0


def _peak_labels(pattern):
    mono_index = int(np.argmin(np.abs(pattern.mz - pattern.monoisotopic_mz)))
    base_offset = pattern.offsets[mono_index]
    labels = []
    for offset in pattern.offsets:
        k = int(offset - base_offset)
        labels.append("M" if k == 0 else f"M{k:+d}")
    return labels


def _tallest_near(peaks_mz, peaks_h, target, tol_ppm):
    if len(peaks_mz) == 0:
        return None
    window = target * tol_ppm * 1e-6
    lo, hi = np.searchsorted(peaks_mz, [target - window, target + window])
    if hi <= lo:
        return None
    return lo + int(np.argmax(peaks_h[lo:hi]))


def match_adduct(counts, adduct, mz, y, fwhm=DEFAULT_FWHM, tolerance_ppm=50.0, min_detect_percent=0.5):
    """1つの付加イオンについて、計算パターンと実測を照合する。

    fwhm は計算パターンの半値全幅(m/z)。min_detect_percent は、計算で最も強いピークに対応する実測ピークが、
    スペクトル全体の最大に対して何 % 以上あれば「検出」とするか。
    """
    fwhm = float(fwhm)
    pattern = ion_pattern(counts, adduct)
    mz = np.asarray(mz, dtype=float)
    y = np.asarray(y, dtype=float)
    peaks = find_peaks(mz, y)
    peaks_mz = np.array([p.mz for p in peaks])
    peaks_h = np.array([p.height for p in peaks])
    spectrum_max = float(np.nanmax(y)) if len(y) else 0.0

    labels = _peak_labels(pattern)
    matches = [PeakMatch(label, float(m), float(r))
               for label, m, r in zip(labels, pattern.mz, pattern.relative) if r >= TABLE_MIN_RELATIVE]
    base = max(matches, key=lambda p: p.calc_relative)

    if len(mz) == 0 or not (mz.min() <= base.calc_mz <= mz.max()):
        return AdductMatch(pattern, False, "測定範囲外", fwhm, 0.0, matches)

    j = _tallest_near(peaks_mz, peaks_h, base.calc_mz, tolerance_ppm)
    found = j is not None and spectrum_max > 0 and peaks_h[j] >= spectrum_max * min_detect_percent / 100.0
    if not found:
        return AdductMatch(pattern, False, "未検出", fwhm, 0.0, matches)

    base_height = peaks_h[j]
    for match in matches:
        k = _tallest_near(peaks_mz, peaks_h, match.calc_mz, tolerance_ppm)
        if k is None:
            continue
        match.measured_mz = float(peaks_mz[k])
        match.measured_height = float(peaks_h[k])
        match.measured_relative = float(peaks_h[k] / base_height * base.calc_relative)
    return AdductMatch(pattern, True, "検出", fwhm, float(base_height / base.calc_relative * 100.0), matches)


def match_all(counts, adducts, mz, y, **kwargs):
    """表記の並びそれぞれを照合する。読めない表記や組成の誤りは、その位置に FormulaError を入れて返す。"""
    results = []
    for adduct in adducts:
        try:
            parsed = parse_adduct(adduct) if isinstance(adduct, str) else adduct
            results.append(match_adduct(counts, parsed, mz, y, **kwargs))
        except FormulaError as e:
            results.append(e)
    return results


TABLE_COLUMNS = ["付加イオン", "イオンの組成", "ピーク", "計算 m/z", "実測 m/z", "誤差 (ppm)",
                 "計算 相対強度 (%)", "実測 相対強度 (%)", "実測 強度", "半値全幅", "状態"]


def results_table(results, adducts_text=None):
    """照合結果の表。数値は丸めない(CSV で精度を落とさないため)。"""
    from .chemistry import format_formula
    rows = []
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            name = adducts_text[i] if adducts_text else ""
            rows.append({"付加イオン": name, "状態": f"エラー: {result}"})
            continue
        for match in result.peaks:
            rows.append({
                "付加イオン": result.pattern.adduct.notation,
                "イオンの組成": format_formula(result.pattern.composition),
                "ピーク": match.label,
                "計算 m/z": match.calc_mz,
                "実測 m/z": match.measured_mz,
                "誤差 (ppm)": match.error_ppm,
                "計算 相対強度 (%)": match.calc_relative,
                "実測 相対強度 (%)": match.measured_relative,
                "実測 強度": match.measured_height,
                "半値全幅": result.fwhm,
                "状態": result.status,
            })
    return pd.DataFrame(rows, columns=TABLE_COLUMNS)
