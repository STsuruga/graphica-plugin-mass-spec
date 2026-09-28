"""analyzer「同位体パターンと照合」と、mzML・BAF の importer。"""
import os

import numpy as np
import pandas as pd

from graphica.plugin import AnalysisResult

from .chemistry import (
    NEGATIVE_STANDARD, POSITIVE_STANDARD, FormulaError, average_mass, format_formula,
    monoisotopic_mass, parse_formula, split_adducts,
)
from .datasets import INTENSITY_COL, MZ_COL, TIC_COL, TIME_COL, fallback_colors, pattern_datasets
from .matching import match_all, results_table
from .baf import BafUnsupported, read_baf
from .mzml import MzmlError, read_mzml
from .spectra import LabelFormat, compact_zeros

ANALYZER_NAME = "同位体パターンと照合"
PRESET_POSITIVE = "正の標準 (" + " ".join(POSITIVE_STANDARD) + ")"
PRESET_NEGATIVE = "負の標準 (" + " ".join(NEGATIVE_STANDARD) + ")"
PRESET_NONE = "下の欄だけ"

PARAM_SCHEMA = [
    {"name": "formula", "label": "組成式(中性分子)", "type": "str", "default": ""},
    {"name": "preset", "label": "付加イオン", "type": "choice", "default": PRESET_POSITIVE,
     "choices": [PRESET_POSITIVE, PRESET_NEGATIVE, PRESET_NONE]},
    {"name": "extra_adducts", "label": "追加の付加イオン(例: [2M+Na]+, [M+2H]2+)", "type": "str", "default": ""},
    {"name": "fwhm", "label": "計算パターンの半値全幅 (m/z)", "type": "float", "default": 0.1,
     "min": 0.0001, "max": 100.0, "decimals": 4},
    {"name": "tolerance_ppm", "label": "探す幅 (±ppm)", "type": "float", "default": 50.0,
     "min": 1.0, "max": 5000.0, "decimals": 1},
    {"name": "detect_percent", "label": "検出の下限(最大ピークに対する %)", "type": "float", "default": 0.5,
     "min": 0.0, "max": 100.0, "decimals": 2},
    {"name": "mz_decimals", "label": "m/z の小数点以下の桁数", "type": "int", "default": 4, "min": 0, "max": 6},
    {"name": "ppm_decimals", "label": "ppm の小数点以下の桁数", "type": "int", "default": 1, "min": 0, "max": 3},
    {"name": "overlay", "label": "計算パターンを重ねる", "type": "bool", "default": True},
]


def adducts_from_params(params):
    preset = params.get("preset", PRESET_POSITIVE)
    adducts = list(POSITIVE_STANDARD if preset == PRESET_POSITIVE else
                   NEGATIVE_STANDARD if preset == PRESET_NEGATIVE else ())
    extra = split_adducts(params.get("extra_adducts", "") or "")
    return adducts + [a for a in extra if a not in adducts]


def analyze(dataset, params):
    formula_text = (params.get("formula") or "").strip()
    try:
        counts = parse_formula(formula_text)
    except FormulaError as e:
        raise ValueError(str(e)) from e
    adducts = adducts_from_params(params)
    if not adducts:
        raise ValueError("付加イオンがありません。標準のセットを選ぶか、表記を入力してください")
    x = np.asarray(pd.to_numeric(pd.Series(dataset.x_data), errors="coerce"), dtype=float)
    y = np.asarray(pd.to_numeric(pd.Series(dataset.y_data), errors="coerce"), dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        raise ValueError("数値の点が少なすぎます。X が m/z、Y が強度のスペクトルを選んでください")
    order = np.argsort(x[ok], kind="stable")
    x, y = x[ok][order], y[ok][order]

    results = match_all(counts, adducts, x, y, fwhm=float(params.get("fwhm", 0.1) or 0.1),
                        tolerance_ppm=float(params.get("tolerance_ppm", 50.0)),
                        min_detect_percent=float(params.get("detect_percent", 0.5)))
    fmt = LabelFormat(mz_decimals=int(params.get("mz_decimals", 4)), ppm_decimals=int(params.get("ppm_decimals", 1)))
    table = results_table(results, adducts)
    summary = pd.DataFrame([{
        "付加イオン": "M(中性分子)", "イオンの組成": format_formula(counts), "ピーク": "M",
        "計算 m/z": monoisotopic_mass(counts), "状態": f"平均分子量 {average_mass(counts):.4f}",
    }], columns=table.columns)
    table = pd.concat([summary, table], ignore_index=True)

    provenance = {"plugin": "mass_spec", "analyzer": ANALYZER_NAME, "formula": format_formula(counts),
                  "adducts": adducts, "source_dataset": dataset.name}
    colors = [c for c in fallback_colors() if c.lower() != (dataset.color or "").lower()]
    new_datasets = []
    annotations = []
    found = [r for r in results if not isinstance(r, Exception) and r.found]
    for i, result in enumerate(found):
        color = colors[i % len(colors)]
        if params.get("overlay", True):
            new_datasets.extend(pattern_datasets(result, format_formula(counts), color, fmt,
                                                 provenance=dict(provenance, adduct=result.pattern.adduct.notation,
                                                                 fwhm=result.fwhm)))
        anchor = result.monoisotopic if result.monoisotopic and np.isfinite(result.monoisotopic.measured_mz) \
            else result.base
        text = f"{result.pattern.adduct.notation} {fmt.mz(anchor.measured_mz)} ({fmt.ppm(anchor.error_ppm)})"
        annotations.append({"type": "text", "xy": [anchor.measured_mz, anchor.measured_height], "text": text,
                            "color": color})
    return AnalysisResult(table=table, annotations=annotations or None, new_datasets=new_datasets or None)


def load_mzml_file(path):
    """importer: スキャンが1つならそのスペクトル、複数なら TIC。"""
    try:
        run = read_mzml(path)
    except MzmlError as e:
        raise ValueError(str(e)) from e
    return _run_to_frame(run)


def load_baf_file(path):
    """importer: .d の中の analysis.baf を選んだとき。フォルダごと読む。"""
    try:
        run = read_baf(os.path.dirname(os.path.abspath(path)))
    except BafUnsupported as e:
        raise ValueError(f"{e}。MS ビューアの「.d を開く」なら msconvert で変換して開けます") from e
    return _run_to_frame(run)


def _run_to_frame(run):
    scans = run.ms1_scans() or run.scans
    if len(scans) == 1:
        mz, y = compact_zeros(*scans[0].arrays())
        return pd.DataFrame({MZ_COL: mz, INTENSITY_COL: y})
    return pd.DataFrame({TIME_COL: [s.time_min for s in scans], TIC_COL: [s.tic for s in scans]})
