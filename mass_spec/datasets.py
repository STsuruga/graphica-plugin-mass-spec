"""本体へ渡す Dataset を作る。パネル・analyzer・importer で共通に使う。"""
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import to_hex

from graphica.plugin import Dataset

from .plot_types import STICK
from .spectra import LabelFormat, cluster_heads, compact_zeros, select_label_peaks

PLUGIN_NAME = "mass_spec"
MZ_COL = "m/z"
INTENSITY_COL = "強度"
TIME_COL = "時間 (min)"
TIC_COL = "TIC"
LABEL_COL = "ラベル"
# 本体は点数が上限(既定 1000)を超えるデータセットに点のラベルを描かないので、ラベル付きはこれ以下にする
MAX_LABELED_POINTS = 1000


def fallback_colors():
    """ctx が無いところ(analyzer)で使う色。Qt に渡すこともあるので #rrggbb にしておく。"""
    return [to_hex(c) for c in colormaps["tab10"].colors]


def _dataset(name, df, x_col, y_col, color, source_file=None, provenance=None, **extra):
    ds = Dataset(name=name, df=df, x_col_name=x_col, y_col_name=y_col)
    ds.color = color
    ds.source_plugin = PLUGIN_NAME
    ds.source_file = source_file
    ds.provenance = provenance
    for key, value in extra.items():
        setattr(ds, key, value)
    return ds


def tic_dataset(run, color, provenance=None):
    df = pd.DataFrame({TIME_COL: run.times(), TIC_COL: run.tics()})
    return _dataset(f"{run.name} TIC", df, TIME_COL, TIC_COL, color, run.path, provenance)


def spectrum_dataset(name, mz, y, color, source_file=None, provenance=None, x_range=None, compact=True):
    """profile の線。x_range があればその m/z 範囲だけにする。0 が続く区間は形を変えずに詰める。"""
    mz = np.asarray(mz, dtype=float)
    y = np.asarray(y, dtype=float)
    if x_range is not None:
        keep = (mz >= min(x_range)) & (mz <= max(x_range))
        mz, y = mz[keep], y[keep]
    if compact:
        mz, y = compact_zeros(mz, y)
    df = pd.DataFrame({MZ_COL: mz, INTENSITY_COL: y})
    return _dataset(name, df, MZ_COL, INTENSITY_COL, color, source_file, provenance)


def stick_dataset(name, mz, heights, labels, color, source_file=None, provenance=None, labels_only=False):
    """棒と点のラベル。labels は各行の文字列(空文字の行には本体がラベルを描かない)。"""
    df = pd.DataFrame({MZ_COL: np.asarray(mz, dtype=float), INTENSITY_COL: np.asarray(heights, dtype=float),
                       LABEL_COL: list(labels)})
    return _dataset(name, df, MZ_COL, INTENSITY_COL, color, source_file, provenance,
                    plot_type=STICK, show_point_labels=any(labels), point_label_col_name=LABEL_COL,
                    linewidth=0.0 if labels_only else 1.5)


def centroid_dataset(name, peaks, color, fmt=None, label_top_n=None, label_min_relative=20.0,
                     min_relative=1.0, source_file=None, provenance=None, labels_only=False, pinned=()):
    """ピークの一覧を Stick にする。min_relative(最大に対する %)未満は省き、多すぎれば強い順に上限まで。"""
    fmt = fmt or LabelFormat()
    if not peaks:
        raise ValueError("ピークがありません")
    base = max(p.height for p in peaks)
    kept = [p for p in peaks if p.height >= base * min_relative / 100.0]
    if len(kept) > MAX_LABELED_POINTS:
        kept = sorted(sorted(kept, key=lambda p: p.height, reverse=True)[:MAX_LABELED_POINTS], key=lambda p: p.mz)
    # 本体の点のラベルは重なりを避けないので、同位体でまとめたシグナルごとにいちばん強いピーク(と固定したもの)だけ
    candidates = select_label_peaks(kept, top_n=label_top_n, min_relative=label_min_relative, pinned=pinned)
    pinned_peaks = [min(kept, key=lambda p: abs(p.mz - m)) for m in pinned]
    labeled = {id(p) for p in cluster_heads(candidates)} | {id(p) for p in pinned_peaks}
    labels = [fmt.label(p.mz, p.height, base) if id(p) in labeled else "" for p in kept]
    return stick_dataset(name, [p.mz for p in kept], [p.height for p in kept], labels, color,
                         source_file, provenance, labels_only)


def pattern_datasets(match, formula_text, color, fmt=None, label_min_relative=5.0, provenance=None):
    """照合した付加イオンの計算パターン: 実測に合わせた profile の線と、ラベル付きの棒。"""
    fmt = fmt or LabelFormat()
    notation = match.pattern.adduct.notation
    base_name = f"{formula_text} {notation} 計算"
    px, py = match.profile()
    profile = spectrum_dataset(f"{base_name}(FWHM {match.fwhm:g})", px, py, color,
                               provenance=provenance, compact=False)
    profile.linestyle = "--"
    sx, sy = match.sticks()
    labels = [fmt.mz(m) if rel >= label_min_relative else "" for m, rel in zip(sx, match.pattern.relative)]
    sticks = stick_dataset(base_name, sx, sy, labels, color, provenance=provenance)
    return [profile, sticks]


def range_text(t0, t1):
    return f"{min(t0, t1):.2f}–{max(t0, t1):.2f} min"
