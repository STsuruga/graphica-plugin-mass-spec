"""plot type「Stick」(棒のスペクトル)。"""
import numpy as np

STICK = "Stick"

_NO_LINE = {"None", "none", "", " "}


def draw_stick(dataset, ax, x, y):
    """0 から各点まで縦線を引く。線の太さが 0 か線種が「なし」なら棒を描かず、点のラベルだけが残る。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if dataset.linewidth <= 0 or dataset.linestyle in _NO_LINE:
        # 凡例と自動の表示範囲のために、見えない点だけ置く
        (artist,) = ax.plot(x[ok], y[ok], linestyle="none", marker="none", color=dataset.color, label=dataset.name)
        return artist
    return ax.vlines(x[ok], 0.0, y[ok], colors=dataset.color, linewidths=dataset.linewidth,
                     linestyles=dataset.linestyle, alpha=dataset.alpha, label=dataset.name)
