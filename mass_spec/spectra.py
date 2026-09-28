"""スペクトルの平均・背景の差し引き・ピーク検出・ラベル。GUI にも Graphica 本体にも依存しない。"""
from dataclasses import dataclass

import numpy as np

# 点の番号ごとに平均してよい、スキャン間の m/z の格子のずれの上限(ppm)。
# micrOTOF の実測では 0.4 ppm 以下で、点の間隔(約 20 ppm)より十分小さい。
GRID_TOLERANCE_PPM = 2.0
# 半値を探す片側の点数の上限。micrOTOF の profile では半値全幅は 2〜6 点。上限が無いと、
# 持ち上がったベースライン上の極大ごとに数万点を歩いて、読み込みが 10 秒を超えた。
FWHM_SEARCH_POINTS = 100


def average_spectrum(scans):
    """スキャンの平均スペクトル (m/z, 強度)。

    点数が同じで格子のずれが小さければ点の番号ごとに平均し、そうでなければ最初のスキャンの格子へ線形補間してから平均する。
    """
    if not scans:
        raise ValueError("平均するスキャンがありません")
    ref_mz, ref_y = scans[0].arrays()
    mz_sum = ref_mz.copy()
    y_sum = ref_y.copy()
    for scan in scans[1:]:
        mz, y = scan.arrays()
        if len(mz) == len(ref_mz) and len(mz) and np.max(np.abs(mz - ref_mz) / ref_mz) * 1e6 <= GRID_TOLERANCE_PPM:
            mz_sum += mz
            y_sum += y
        else:
            mz_sum += ref_mz
            y_sum += np.interp(ref_mz, mz, y, left=0.0, right=0.0) if len(mz) else 0.0
    n = len(scans)
    return mz_sum / n, y_sum / n


def subtract_background(mz, y, bg_mz, bg_y, clip_negative=True):
    """背景のスペクトルを試料の格子に合わせて引く。"""
    mz, y, bg_mz, bg_y = (np.asarray(a, dtype=float) for a in (mz, y, bg_mz, bg_y))
    if len(bg_mz) == len(mz) and np.allclose(bg_mz, mz, rtol=GRID_TOLERANCE_PPM * 1e-6, atol=0):
        bg = bg_y
    else:
        bg = np.interp(mz, bg_mz, bg_y, left=0.0, right=0.0)
    result = y - bg
    if clip_negative:
        result = np.clip(result, 0.0, None)
    return result


def compact_zeros(mz, y):
    """0 が続く区間の内側の点を省く(ピークの裾の 0 は残す)。線で描いた形は変わらない。"""
    y = np.asarray(y)
    if len(y) < 3:
        return np.asarray(mz), y
    zero = y == 0
    inner = zero.copy()
    inner[1:] &= zero[:-1]
    inner[:-1] &= zero[1:]
    inner[0] = inner[-1] = False
    keep = ~inner
    return np.asarray(mz)[keep], y[keep]


@dataclass
class Peak:
    mz: float           # 頂点の m/z(ガウスの3点補間)
    height: float       # 頂点の強度
    fwhm: float         # 半値全幅(求められなければ nan)
    index: int          # 元の配列での最大点の番号

    @property
    def resolution(self):
        return self.mz / self.fwhm if self.fwhm and np.isfinite(self.fwhm) and self.fwhm > 0 else float("nan")


def is_centroid_like(mz, y):
    """棒のスペクトル(点の間隔がピークの幅より広い)らしいか。"""
    mz = np.asarray(mz, dtype=float)
    if len(mz) < 3:
        return True
    rel = np.diff(mz) / mz[1:]
    return float(np.median(rel)) > 1e-4


def _apex(mz, y, i):
    """最大点とその両隣から、ガウス(対数の放物線)で頂点を補間する。"""
    if 0 < i < len(y) - 1 and y[i - 1] > 0 and y[i + 1] > 0 and y[i] > 0:
        a, b, c = np.log(y[i - 1]), np.log(y[i]), np.log(y[i + 1])
        denom = a - 2 * b + c
        if denom < 0:
            offset = 0.5 * (a - c) / denom
            if abs(offset) <= 1:
                step = (mz[i + 1] - mz[i - 1]) / 2
                x = mz[i] + offset * step
                h = np.exp(b - 0.25 * (a - c) * offset)
                return float(x), float(h)
    return float(mz[i]), float(y[i])


def _fwhm(mz, y, i, height, max_points=FWHM_SEARCH_POINTS):
    """半値全幅。半値まで下がらないピーク(持ち上がったベースライン上など)は nan。"""
    half = height / 2.0
    lo = i
    while lo > 0 and y[lo] > half and i - lo < max_points:
        lo -= 1
    hi = i
    while hi < len(y) - 1 and y[hi] > half and hi - i < max_points:
        hi += 1
    if y[lo] > half or y[hi] > half:
        return float("nan")
    left = np.interp(half, [y[lo], y[lo + 1]], [mz[lo], mz[lo + 1]])
    right = np.interp(half, [y[hi], y[hi - 1]], [mz[hi], mz[hi - 1]])
    return float(right - left)


def find_peaks(mz, y, min_height=0.0, centroid=None):
    """ピークの一覧(m/z の昇順)。centroid が None なら点の間隔から判断する。"""
    mz = np.asarray(mz, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(mz) == 0:
        return []
    if centroid is None:
        centroid = is_centroid_like(mz, y)
    if centroid:
        idx = np.nonzero((y > 0) & (y >= min_height))[0]
        return [Peak(float(mz[i]), float(y[i]), float("nan"), int(i)) for i in idx]
    if len(y) < 3:
        return []
    # 頂点が平らな場合も1つのピークにするため、左は >、右は >= で比べる
    core = (y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > 0) & (y[1:-1] >= min_height)
    peaks = []
    for i in np.nonzero(core)[0] + 1:
        x, h = _apex(mz, y, i)
        peaks.append(Peak(x, h, _fwhm(mz, y, i, y[i]), int(i)))
    return peaks


def centroid_arrays(peaks):
    return (np.array([p.mz for p in peaks], dtype=float),
            np.array([p.height for p in peaks], dtype=float))


def select_label_peaks(peaks, x_range=None, top_n=10, min_relative=None, min_separation_ppm=500.0, pinned=()):
    """ラベルを付けるピーク。表示範囲の中で強い順に選び、近すぎるものは強い方を残す。

    top_n と min_relative(表示範囲の最大に対する %)はどちらか一方でも両方でもよい。pinned の m/z に最も近いピークは必ず入れる。
    """
    candidates = [p for p in peaks if x_range is None or x_range[0] <= p.mz <= x_range[1]]
    if not candidates:
        return []
    base = max(p.height for p in candidates)
    ordered = sorted(candidates, key=lambda p: p.height, reverse=True)
    chosen = []
    for target in pinned:
        near = min(candidates, key=lambda p: abs(p.mz - target))
        if near not in chosen:
            chosen.append(near)
    for p in ordered:
        if top_n is not None and len(chosen) >= top_n + len(pinned):
            break
        if min_relative is not None and p.height < base * min_relative / 100.0:
            break
        if any(abs(p.mz - q.mz) / q.mz * 1e6 < min_separation_ppm for q in chosen):
            continue
        chosen.append(p)
    return sorted(chosen, key=lambda p: p.mz)


@dataclass
class LabelFormat:
    mz_decimals: int = 4
    intensity: str = "none"     # "none" / "relative" / "absolute"
    intensity_decimals: int = 1
    ppm_decimals: int = 1

    def mz(self, value):
        return f"{value:.{self.mz_decimals}f}"

    def ppm(self, value):
        return f"{value:+.{self.ppm_decimals}f} ppm"

    def label(self, mz, height=None, base=None):
        text = self.mz(mz)
        if self.intensity == "relative" and height is not None and base:
            text += f" ({height / base * 100:.{self.intensity_decimals}f}%)"
        elif self.intensity == "absolute" and height is not None:
            text += f" ({height:.{max(self.intensity_decimals, 1)}e})"
        return text
