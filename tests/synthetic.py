"""テスト用の合成データ。ユーザーの測定データは使わない。"""
import base64
import zlib

import numpy as np

from mass_spec.chemistry import ion_pattern, parse_adduct, parse_formula


def tof_grid(lo=100.0, hi=1000.0, ppm_at_500=20.0):
    """√(m/z) が等間隔の格子(TOF と同じく、m/z が大きいほど間隔が広い)。"""
    step = ppm_at_500 * 1e-6 * np.sqrt(500.0) / 2.0
    return np.arange(np.sqrt(lo), np.sqrt(hi), step) ** 2


def gaussian_peaks(grid, centers, heights, resolution):
    y = np.zeros_like(grid)
    for c, h in zip(centers, heights):
        sigma = c / resolution / (2 * np.sqrt(2 * np.log(2)))
        near = np.abs(grid - c) < 8 * sigma
        y[near] += h * np.exp(-0.5 * ((grid[near] - c) / sigma) ** 2)
    return y


def ion_spectrum(grid, formula, adduct, height=1000.0, resolution=12000, shift_ppm=0.0):
    pattern = ion_pattern(parse_formula(formula), parse_adduct(adduct))
    centers = pattern.mz * (1 + shift_ppm * 1e-6)
    return gaussian_peaks(grid, centers, pattern.relative / 100 * height, resolution), pattern


def _encode(values, dtype, compress):
    raw = np.asarray(values, dtype=dtype).tobytes()
    if compress:
        raw = zlib.compress(raw)
    return base64.b64encode(raw).decode("ascii")


def mzml_text(scans, mz_dtype="<f8", int_dtype="<f4", compress=True, time_unit="second",
              indexed=False, polarity="positive", profile=True, with_tic=True, param_group=False):
    """scans: [(時刻, m/z の配列, 強度の配列), ...]。時刻は time_unit の単位。"""
    acc = {"<f8": ("MS:1000523", "64-bit float"), "<f4": ("MS:1000521", "32-bit float")}
    comp = ("MS:1000574", "zlib compression") if compress else ("MS:1000576", "no compression")
    pol = {"positive": ("MS:1000130", "positive scan"), "negative": ("MS:1000129", "negative scan")}[polarity]
    kind = ("MS:1000128", "profile spectrum") if profile else ("MS:1000127", "centroid spectrum")
    unit = {"second": ("UO:0000010", "second"), "minute": ("UO:0000031", "minute")}[time_unit]
    lines = ['<?xml version="1.0" encoding="utf-8"?>']
    if indexed:
        lines.append('<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">')
    lines.append('<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">')
    if param_group:
        lines.append('<referenceableParamGroupList count="1"><referenceableParamGroup id="common">'
                     f'<cvParam cvRef="MS" accession="{pol[0]}" name="{pol[1]}" value=""/>'
                     f'<cvParam cvRef="MS" accession="{kind[0]}" name="{kind[1]}" value=""/>'
                     '</referenceableParamGroup></referenceableParamGroupList>')
    lines.append(f'<run id="r"><spectrumList count="{len(scans)}">')
    for i, (t, mz, y) in enumerate(scans):
        lines.append(f'<spectrum index="{i}" id="scan={i + 1}" defaultArrayLength="{len(mz)}">')
        lines.append('<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>')
        if param_group:
            lines.append('<referenceableParamGroupRef ref="common"/>')
        else:
            lines.append(f'<cvParam cvRef="MS" accession="{pol[0]}" name="{pol[1]}" value=""/>')
            lines.append(f'<cvParam cvRef="MS" accession="{kind[0]}" name="{kind[1]}" value=""/>')
        if with_tic:
            lines.append(f'<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{float(np.sum(y)):.6g}"/>')
        lines.append('<scanList count="1"><scan>'
                     f'<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{t}" '
                     f'unitCvRef="UO" unitAccession="{unit[0]}" unitName="{unit[1]}"/></scan></scanList>')
        lines.append('<binaryDataArrayList count="2">')
        for values, dtype, array_acc in ((mz, mz_dtype, ("MS:1000514", "m/z array")),
                                         (y, int_dtype, ("MS:1000515", "intensity array"))):
            text = _encode(values, dtype, compress)
            lines.append(f'<binaryDataArray encodedLength="{len(text)}">'
                         f'<cvParam cvRef="MS" accession="{acc[dtype][0]}" name="{acc[dtype][1]}" value=""/>'
                         f'<cvParam cvRef="MS" accession="{comp[0]}" name="{comp[1]}" value=""/>'
                         f'<cvParam cvRef="MS" accession="{array_acc[0]}" name="{array_acc[1]}" value=""/>'
                         f'<binary>{text}</binary></binaryDataArray>')
        lines.append('</binaryDataArrayList></spectrum>')
    lines.append('</spectrumList></run></mzML>')
    if indexed:
        lines.append('<indexListOffset>0</indexListOffset></indexedmzML>')
    return "\n".join(lines)


def write_mzml(path, scans, **kwargs):
    with open(path, "w", encoding="utf-8") as f:
        f.write(mzml_text(scans, **kwargs))
    return str(path)


def run_scans(n=10, peak_scans=range(3, 7), shift_ppm_per_scan=0.02, formula="C6H12O6", adduct="[M+Na]+",
              lo=150.0, hi=450.0):
    """n スキャンの合成測定。peak_scans のスキャンだけに試料のイオンが出て、ほかは弱い背景だけ。"""
    grid = tof_grid(lo, hi)
    background = gaussian_peaks(grid, [300.0], [50.0], 12000)
    ion, _pattern = ion_spectrum(grid, formula, adduct)
    scans = []
    for k in range(n):
        mz = grid * (1 + k * shift_ppm_per_scan * 1e-6)
        y = background + (ion if k in peak_scans else 0.0)
        scans.append((float(k), mz, y))
    return scans
