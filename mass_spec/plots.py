"""MS ビューアのグラフ(TIC・スペクトル・3段目の枠)。マウス操作は README の表のとおり。

matplotlib の標準のツールバーは使わず、イベントを自前で処理する。本体のプロットにはプラグインから操作を足せないので、
範囲の選択や拡大はここで行い、結果を本体へ転送する。
"""
import numpy as np
from matplotlib.backend_bases import MouseButton
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QCursor, QPalette

from .spectra import LabelFormat, Peak, select_label_peaks

WHEEL_FACTOR = 1.25
CLICK_PIXELS = 4       # これより動かなければドラッグでなくクリック
EDGE_PIXELS = 6        # 範囲の端をつかめる距離
PICK_PIXELS = 10       # ピークを選べる距離
ISOTOPE_SPACING = 1.00335
# 本体はグラフの文字ごとにこの候補を渡している(matplotlib 全体の設定は変えていない)。日本語を含む文字に使う
JP_FONTS = ["DejaVu Sans", "Yu Gothic", "Hiragino Sans", "Hiragino Kaku Gothic ProN", "Meiryo", "MS Gothic",
            "Noto Sans CJK JP"]
HIGHLIGHT = "#2563EB"
MEASURE_COLOR = "#1F6F78"


def _qcolor_hex(color):
    return to_hex((color.redF(), color.greenF(), color.blueF()))


class InteractivePlot(FigureCanvasQTAgg):
    """軸の上: ホイールで拡大縮小、左ドラッグでずらす。グラフの中: ホイールで拡大縮小(Shift で縦)、中ボタンでパン。
    ダブルクリックで全体、Backspace で1つ前の表示。"""

    context_requested = Signal(object)   # 右クリックの画面座標(QPoint)
    view_changed = Signal()
    x_range_changed = Signal(float, float)

    def __init__(self, height_inches=2.0, parent=None):
        self.figure = Figure(figsize=(5, height_inches), layout="constrained")
        super().__init__(self.figure)
        self.setParent(parent)
        self.ax = self.figure.add_subplot(111)
        self._home = None
        self._history = []
        self._pan = None
        self._axis_drag = None
        self._last_wheel = False
        self._cursor_region = None
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumHeight(int(height_inches * 45))
        self.mpl_connect("button_press_event", self._on_press)
        self.mpl_connect("button_release_event", self._on_release)
        self.mpl_connect("motion_notify_event", self._on_motion)
        self.mpl_connect("scroll_event", self._on_scroll)
        self.mpl_connect("key_press_event", self._on_key)
        self.apply_palette()

    # ---- 見た目
    def apply_palette(self):
        """本体のテーマ(Qt のパレット)に合わせる。"""
        pal = self.palette()
        bg = _qcolor_hex(pal.color(QPalette.ColorRole.Base))
        fg = _qcolor_hex(pal.color(QPalette.ColorRole.Text))
        self.figure.set_facecolor(bg)
        self.ax.set_facecolor(bg)
        for spine in self.ax.spines.values():
            spine.set_color(fg)
        self.ax.tick_params(colors=fg, labelsize=8)
        self.ax.xaxis.label.set_color(fg)
        self.ax.yaxis.label.set_color(fg)
        self.foreground = fg

    # ---- 表示範囲
    def set_home(self, xlim, ylim):
        self._home = (tuple(xlim), tuple(ylim))
        self._history.clear()
        self.ax.set_xlim(*xlim)
        self.ax.set_ylim(*ylim)

    def _push_history(self):
        self._history.append((self.ax.get_xlim(), self.ax.get_ylim()))
        del self._history[:-50]

    def set_view(self, xlim=None, ylim=None, remember=True):
        if remember:
            self._push_history()
        if xlim is not None:
            self.ax.set_xlim(*xlim)
        if ylim is not None:
            self.ax.set_ylim(*ylim)
        elif xlim is not None:
            self.fit_y()
        self._view_updated()

    def reset_view(self):
        if self._home is None:
            return
        self._push_history()
        self.ax.set_xlim(*self._home[0])
        self.ax.set_ylim(*self._home[1])
        self._view_updated()

    def back(self):
        if not self._history:
            return
        xlim, ylim = self._history.pop()
        self.ax.set_xlim(*xlim)
        self.ax.set_ylim(*ylim)
        self._view_updated()

    def fit_y(self):
        """横の範囲を変えたあとの縦の合わせ方。既定は何もしない。"""

    def _view_updated(self):
        self.view_changed.emit()
        self.x_range_changed.emit(*self.ax.get_xlim())
        self.draw_idle()

    # ---- 位置
    def axis_region(self, event):
        """"x"(横軸の目盛りの帯)/ "y"(縦軸の帯)/ None。"""
        bbox = self.ax.bbox
        if bbox.x0 <= event.x <= bbox.x1 and event.y < bbox.y0:
            return "x"
        if bbox.y0 <= event.y <= bbox.y1 and event.x < bbox.x0:
            return "y"
        return None

    def _data_x(self, x_pixel):
        return self.ax.transData.inverted().transform((x_pixel, self.ax.bbox.y0))[0]

    def data_x_to_pixel(self, x):
        return self.ax.transData.transform((x, 0))[0]

    # ---- イベント
    def _zoom_x(self, center, factor):
        x0, x1 = self.ax.get_xlim()
        self.ax.set_xlim(center - (center - x0) * factor, center + (x1 - center) * factor)
        self.fit_y()

    def _zoom_y(self, factor):
        y0, y1 = self.ax.get_ylim()
        self.ax.set_ylim(y0, y0 + (y1 - y0) * factor)

    def _on_scroll(self, event):
        region = self.axis_region(event)
        if event.inaxes is not self.ax and region is None:
            return
        if not self._last_wheel:
            self._push_history()
        self._last_wheel = True
        factor = 1 / WHEEL_FACTOR if event.button == "up" else WHEEL_FACTOR
        if region == "y" or (region is None and event.key == "shift"):
            self._zoom_y(factor)
        else:
            self._zoom_x(event.xdata if event.inaxes is self.ax else self._data_x(event.x), factor)
        self._view_updated()

    def _on_press(self, event):
        self._last_wheel = False
        region = self.axis_region(event)
        if region is not None and event.button == MouseButton.LEFT:
            if event.dblclick:
                self.reset_view()
                return
            self._push_history()
            self._axis_drag = (region, event.x, event.y, self.ax.get_xlim(), self.ax.get_ylim())
            return
        if event.inaxes is not self.ax:
            return
        if event.button == MouseButton.MIDDLE:
            self._push_history()
            self._pan = (event.x, event.y, self.ax.get_xlim(), self.ax.get_ylim())
        elif event.button == MouseButton.RIGHT:
            self.context_requested.emit(QCursor.pos())
        elif event.button == MouseButton.LEFT and event.dblclick:
            self.reset_view()
        elif event.button == MouseButton.LEFT:
            self.on_left_press(event)

    def _shift(self, x0, y0, xlim, ylim, event, horizontal=True, vertical=True):
        bbox = self.ax.bbox
        if horizontal:
            dx = (event.x - x0) * (xlim[1] - xlim[0]) / bbox.width
            self.ax.set_xlim(xlim[0] - dx, xlim[1] - dx)
        if vertical:
            dy = (event.y - y0) * (ylim[1] - ylim[0]) / bbox.height
            self.ax.set_ylim(ylim[0] - dy, ylim[1] - dy)

    def _on_motion(self, event):
        if self._axis_drag is not None:
            region, x0, y0, xlim, ylim = self._axis_drag
            self._shift(x0, y0, xlim, ylim, event, horizontal=region == "x", vertical=region == "y")
            if region == "x":
                self.fit_y()
            self._view_updated()
            return
        if self._pan is not None:
            x0, y0, xlim, ylim = self._pan
            self._shift(x0, y0, xlim, ylim, event)
            self._view_updated()
            return
        self._update_cursor(self.axis_region(event))
        self.on_motion(event)

    def _update_cursor(self, region):
        if region == self._cursor_region:
            return
        self._cursor_region = region
        if region == "x":
            self.setCursor(Qt.CursorShape.SizeHorCursor)
        elif region == "y":
            self.setCursor(Qt.CursorShape.SizeVerCursor)
        else:
            self.unsetCursor()

    def _on_release(self, event):
        if self._axis_drag is not None and event.button == MouseButton.LEFT:
            self._axis_drag = None
            return
        if self._pan is not None and event.button == MouseButton.MIDDLE:
            self._pan = None
            return
        if event.button == MouseButton.LEFT:
            self.on_left_release(event)

    def _on_key(self, event):
        if event.key == "backspace":
            self.back()
        elif event.key == "escape":
            self.on_escape()

    # ---- 派生クラスで上書きする
    def on_left_press(self, event):
        pass

    def on_motion(self, event):
        pass

    def on_left_release(self, event):
        pass

    def on_escape(self):
        pass


class TicPlot(InteractivePlot):
    """TIC。左ドラッグで試料、Shift+左ドラッグで背景の時間範囲。帯の端で伸縮、帯の中で移動、クリックで1スキャン。"""

    range_changing = Signal(str, float, float)   # "sample" / "background"、ドラッグ中
    range_changed = Signal(str, float, float)    # ドラッグを終えたとき
    scan_clicked = Signal(float)

    SAMPLE_COLOR = HIGHLIGHT
    BACKGROUND_COLOR = "#7A837F"

    def __init__(self, parent=None):
        super().__init__(1.5, parent)
        self.ranges = {"sample": None, "background": None}
        self._patches = {}
        self._line = None
        self._scan_marker = None
        self._drag = None
        self.ax.set_xlabel("Time (min)", fontsize=8)

    def set_data(self, times, tics):
        self.ax.cla()
        self.apply_palette()
        self._patches = {}
        self._scan_marker = None
        self.ax.set_xlabel("Time (min)", fontsize=8)
        times = np.asarray(times, dtype=float)
        tics = np.asarray(tics, dtype=float)
        (self._line,) = self.ax.plot(times, tics, color=self.foreground, linewidth=1.0)
        self.ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 3))
        if len(times):
            span = max(times[-1] - times[0], 1e-6)
            top = float(np.nanmax(tics)) if np.any(np.isfinite(tics)) else 1.0
            self.set_home((times[0] - span * 0.02, times[-1] + span * 0.02), (0, top * 1.08 or 1.0))
        self.ranges = {"sample": None, "background": None}
        self.draw_idle()

    def set_range(self, kind, t0, t1):
        """範囲を描き直す(数値欄から変えたとき)。None で消す。"""
        if t0 is None or t1 is None:
            self.ranges[kind] = None
        else:
            self.ranges[kind] = (min(t0, t1), max(t0, t1))
        self._draw_range(kind)
        self.draw_idle()

    def _draw_range(self, kind):
        old = self._patches.pop(kind, None)
        if old is not None:
            old.remove()
        r = self.ranges[kind]
        if r is None:
            return
        color = self.SAMPLE_COLOR if kind == "sample" else self.BACKGROUND_COLOR
        self._patches[kind] = self.ax.axvspan(r[0], r[1], color=color, alpha=0.18, linewidth=0)

    def mark_scan(self, t):
        if self._scan_marker is not None:
            self._scan_marker.remove()
            self._scan_marker = None
        if t is not None:
            self._scan_marker = self.ax.axvline(t, color=self.SAMPLE_COLOR, linewidth=0.8, linestyle=":")
        self.draw_idle()

    def _hit(self, event):
        """押した位置にある範囲と、つかんだ部分("left" / "right" / "inside")。"""
        for kind in ("sample", "background"):
            r = self.ranges[kind]
            if r is None:
                continue
            left, right = self.data_x_to_pixel(r[0]), self.data_x_to_pixel(r[1])
            if abs(event.x - left) <= EDGE_PIXELS:
                return kind, "left"
            if abs(event.x - right) <= EDGE_PIXELS:
                return kind, "right"
            if left < event.x < right:
                return kind, "inside"
        return None, None

    def on_left_press(self, event):
        kind, part = (None, None) if event.key == "shift" else self._hit(event)
        if kind is None:
            kind, part = ("background" if event.key == "shift" else "sample"), "new"
        self._drag = {"kind": kind, "part": part, "x0": event.x, "t0": event.xdata,
                      "orig": self.ranges[kind], "moved": False}

    def _dragged_range(self, event):
        d = self._drag
        t = event.xdata if event.xdata is not None else d["t0"]
        if d["part"] == "new":
            return min(d["t0"], t), max(d["t0"], t)
        lo, hi = d["orig"]
        if d["part"] == "left":
            return min(t, hi), max(t, hi)
        if d["part"] == "right":
            return min(lo, t), max(lo, t)
        shift = t - d["t0"]
        return lo + shift, hi + shift

    def on_motion(self, event):
        if self._drag is None or event.inaxes is not self.ax:
            return
        if abs(event.x - self._drag["x0"]) < CLICK_PIXELS and not self._drag["moved"]:
            return
        self._drag["moved"] = True
        kind = self._drag["kind"]
        self.ranges[kind] = self._dragged_range(event)
        self._draw_range(kind)
        self.draw_idle()
        self.range_changing.emit(kind, *self.ranges[kind])

    def on_left_release(self, event):
        d, self._drag = self._drag, None
        if d is None:
            return
        if not d["moved"]:
            if d["t0"] is not None:
                self.mark_scan(d["t0"])
                self.scan_clicked.emit(float(d["t0"]))
            return
        kind = d["kind"]
        self.range_changed.emit(kind, *self.ranges[kind])

    def on_escape(self):
        if self._drag is not None:
            kind = self._drag["kind"]
            self.ranges[kind] = self._drag["orig"]
            self._draw_range(kind)
            self._drag = None
            self.draw_idle()


class SpectrumPlot(InteractivePlot):
    """スペクトル。左ドラッグで m/z 範囲に拡大、Ctrl+左ドラッグで矩形の拡大、Shift+左ドラッグで Δm/z の測定、
    ピークをクリックでラベルの固定、カーソルを合わせるとピークの値を hover_text で知らせる。"""

    hover_text = Signal(str)
    pins_changed = Signal()
    measured = Signal(str)

    def __init__(self, parent=None, height_inches=2.4):
        super().__init__(height_inches, parent)
        self.mz = np.zeros(0)
        self.y = np.zeros(0)
        self.peaks = []
        self.pinned = []
        self.label_mode = "top"        # "top" / "percent" / "none"
        self.label_top_n = 10
        self.label_percent = 5.0
        self.label_format = LabelFormat()
        self.auto_y = True
        self._overlays = []
        self._overlay_artists = []
        self._label_artists = []
        self._drag = None
        self._drag_artist = None
        self._measure_artists = []
        self._line = None
        self.view_changed.connect(self.update_labels)
        self.ax.set_xlabel("m/z", fontsize=8)

    def set_spectrum(self, mz, y, peaks, keep_view=False, color=None, linestyle="-"):
        xlim = self.ax.get_xlim() if keep_view and len(self.mz) else None
        self.mz = np.asarray(mz, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.peaks = list(peaks)
        self.ax.cla()
        self.apply_palette()
        self.ax.set_xlabel("m/z", fontsize=8)
        self._label_artists = []
        self._overlay_artists = []
        self._measure_artists = []
        (self._line,) = self.ax.plot(self.mz, self.y, color=color or self.foreground, linewidth=0.8,
                                     linestyle=linestyle)
        self.ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 3))
        if len(self.mz):
            self.set_home((self.mz[0], self.mz[-1]), (0, self._top(self.mz[0], self.mz[-1])))
            if xlim is not None:
                self.ax.set_xlim(*xlim)
                self.ax.set_ylim(0, self._top(*xlim))
        self._draw_overlays()
        self.after_set_spectrum()
        self.update_labels()

    def after_set_spectrum(self):
        pass

    def clear(self):
        self.set_spectrum([], [], [])

    def _top(self, x0, x1):
        lo, hi = min(x0, x1), max(x0, x1)
        visible = self.y[(self.mz >= lo) & (self.mz <= hi)]
        top = float(np.nanmax(visible)) if len(visible) else 0.0
        for x, y, *_ in self._overlays:
            inside = y[(x >= lo) & (x <= hi)]
            if len(inside):
                top = max(top, float(np.nanmax(inside)))
        return top * 1.15 if top > 0 else 1.0

    def fit_y(self):
        if self.auto_y and len(self.mz):
            self.ax.set_ylim(0, self._top(*self.ax.get_xlim()))

    # ---- 重ね描き(計算パターン)
    def set_overlays(self, overlays):
        """[(x, y, color, label, "stick" / "line"), ...]"""
        self._overlays = [(np.asarray(x, float), np.asarray(y, float), c, lab, kind)
                          for x, y, c, lab, kind in overlays]
        self._draw_overlays()
        self.draw_idle()

    def _draw_overlays(self):
        for artist in self._overlay_artists:
            artist.remove()
        self._overlay_artists = []
        for x, y, color, label, kind in self._overlays:
            if kind == "stick":
                self._overlay_artists.append(self.ax.vlines(x, 0, y, colors=color, linewidths=1.6, label=label))
            else:
                (line,) = self.ax.plot(x, y, color=color, linewidth=1.0, linestyle="--", label=label)
                self._overlay_artists.append(line)

    # ---- ラベル
    def update_labels(self):
        for artist in self._label_artists:
            artist.remove()
        self._label_artists = []
        if self.label_mode != "none" and self.peaks:
            x0, x1 = self.ax.get_xlim()
            chosen = select_label_peaks(
                self.peaks, (x0, x1),
                top_n=self.label_top_n if self.label_mode == "top" else None,
                min_relative=self.label_percent if self.label_mode == "percent" else None,
                pinned=[m for m in self.pinned if x0 <= m <= x1])
            visible = [p.height for p in self.peaks if x0 <= p.mz <= x1]
            base = max(visible) if visible else None
            for p in chosen:
                text = self.label_format.label(p.mz, p.height, base)
                pinned = any(abs(p.mz - m) / p.mz < 5e-6 for m in self.pinned)
                self._label_artists.append(self.ax.annotate(
                    text, (p.mz, p.height), textcoords="offset points", xytext=(0, 3), ha="center",
                    fontsize=7.5, color=self.foreground, fontweight="bold" if pinned else "normal", clip_on=True))
        self.draw_idle()

    def labeled_texts(self):
        return [a.get_text() for a in self._label_artists]

    def nearest_peak(self, event, pixels=PICK_PIXELS):
        return self.peak_near_pixel(event.x, pixels)

    def peak_near_pixel(self, x_pixel, pixels=PICK_PIXELS):
        if not self.peaks or x_pixel is None:
            return None
        xs = np.array([p.mz for p in self.peaks])
        pix = self.ax.transData.transform(np.column_stack([xs, np.zeros_like(xs)]))[:, 0]
        i = int(np.argmin(np.abs(pix - x_pixel)))
        return self.peaks[i] if abs(pix[i] - x_pixel) <= pixels else None

    # ---- マウス
    def on_left_press(self, event):
        mode = "box" if event.key in ("control", "ctrl") else "measure" if event.key == "shift" else "zoom"
        self._drag = {"mode": mode, "x0": event.x, "y0": event.y, "xd": event.xdata, "yd": event.ydata,
                      "moved": False}

    def _clear_drag_artist(self):
        if self._drag_artist is not None:
            self._drag_artist.remove()
            self._drag_artist = None

    def on_motion(self, event):
        if self._drag is None:
            if event.inaxes is self.ax:
                peak = self.nearest_peak(event)
                if peak is not None:
                    self.hover_text.emit(f"m/z {self.label_format.mz(peak.mz)}  強度 {peak.height:.4g}")
                elif event.xdata is not None:
                    self.hover_text.emit(f"m/z {self.label_format.mz(event.xdata)}")
            return
        d = self._drag
        if event.inaxes is not self.ax or event.xdata is None:
            return
        if not d["moved"] and abs(event.x - d["x0"]) < CLICK_PIXELS and abs(event.y - d["y0"]) < CLICK_PIXELS:
            return
        d["moved"] = True
        self._clear_drag_artist()
        if d["mode"] == "zoom":
            self._drag_artist = self.ax.axvspan(d["xd"], event.xdata, color=HIGHLIGHT, alpha=0.15, linewidth=0)
        elif d["mode"] == "box":
            self._drag_artist = self.ax.add_patch(Rectangle(
                (min(d["xd"], event.xdata), min(d["yd"], event.ydata)), abs(event.xdata - d["xd"]),
                abs(event.ydata - d["yd"]), fill=False, edgecolor=HIGHLIGHT, linewidth=1.0))
        else:
            (self._drag_artist,) = self.ax.plot([d["xd"], event.xdata], [d["yd"], d["yd"]], color=MEASURE_COLOR,
                                                linewidth=1.0)
            self.hover_text.emit(self.measure_text(d["xd"], event.xdata))
        self.draw_idle()

    def measure_text(self, x0, x1):
        a, b = sorted((x0, x1))
        delta = b - a
        text = f"Δm/z {delta:.4f}  ({delta / a * 1e6:.1f} ppm)" if a > 0 else f"Δm/z {delta:.4f}"
        if delta > 0:
            z = ISOTOPE_SPACING / delta
            n = round(z)
            if n >= 1 and abs(z - n) / n < 0.1:
                text += f"  → z = {n}"
        return text

    def _snap(self, x_data, x_pixel):
        """測定の端は近くのピークの頂点に吸着させる(同位体の間隔を正確に測るため)。"""
        peak = self.peak_near_pixel(x_pixel)
        return peak.mz if peak is not None else x_data

    def on_left_release(self, event):
        d, self._drag = self._drag, None
        self._clear_drag_artist()
        if d is None:
            return
        if not d["moved"]:
            peak = self.nearest_peak(event)
            if peak is not None:
                self.toggle_pin(peak.mz)
            self.draw_idle()
            return
        x1 = event.xdata if event.xdata is not None else self.ax.get_xlim()[1]
        if d["mode"] == "zoom":
            lo, hi = sorted((d["xd"], x1))
            self.set_view((lo, hi), (0, self._top(lo, hi)))
        elif d["mode"] == "box":
            y1 = event.ydata if event.ydata is not None else self.ax.get_ylim()[1]
            self.set_view(tuple(sorted((d["xd"], x1))), tuple(sorted((d["yd"], y1))))
        else:
            self._draw_measure(self._snap(d["xd"], d["x0"]), self._snap(x1, event.x), d["yd"])
        self.draw_idle()

    def _draw_measure(self, xa, xb, y):
        for artist in self._measure_artists:
            artist.remove()
        text = self.measure_text(xa, xb)
        (line,) = self.ax.plot([xa, xb], [y, y], color=MEASURE_COLOR, linewidth=1.0, marker="|", markersize=8)
        label = self.ax.annotate(text, ((xa + xb) / 2, y), textcoords="offset points", xytext=(0, 4),
                                 ha="center", fontsize=7.5, color=MEASURE_COLOR)
        self._measure_artists = [line, label]
        self.measured.emit(text)

    def on_escape(self):
        for artist in self._measure_artists:
            artist.remove()
        self._measure_artists = []
        self._clear_drag_artist()
        self._drag = None
        self.draw_idle()

    def toggle_pin(self, mz):
        for m in self.pinned:
            if abs(m - mz) / mz < 5e-6:
                self.pinned.remove(m)
                break
        else:
            self.pinned.append(mz)
        self.update_labels()
        self.pins_changed.emit()


class PanePlot(SpectrumPlot):
    """3段目の枠。実測のコピーか計算パターンを1つだけ持つ。"""

    def __init__(self, number, parent=None):
        super().__init__(parent, height_inches=1.5)
        self.number = number
        self.item = None
        self._title = None
        self.show_item(None)

    def show_item(self, item):
        """item: {"kind": "measured" / "calc", "name", "mz", "y", "peaks", "color", "sticks": (x, y) | None, ...}"""
        self.item = item
        self.pinned = []
        if item is None:
            self.set_spectrum([], [], [])
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.draw_idle()
            return
        if item["kind"] == "calc":
            sx, sy = item["sticks"]
            peaks = [Peak(float(m), float(h), float("nan"), i) for i, (m, h) in enumerate(zip(sx, sy))]
            self._overlays = [(np.asarray(sx, float), np.asarray(sy, float), item["color"], item["name"], "stick")]
            self.set_spectrum(item["mz"], item["y"], peaks, color=item["color"], linestyle="--")
        else:
            self._overlays = []
            self.set_spectrum(item["mz"], item["y"], item["peaks"])

    def after_set_spectrum(self):
        name = self.item["name"] if self.item else None
        if name:
            self._set_title(f"{self.number}: {name}")

    def _set_title(self, text):
        self._title = self.ax.text(0.005, 0.97, text, transform=self.ax.transAxes, va="top", ha="left",
                                   fontsize=7.5, color=self.foreground, fontfamily=JP_FONTS, clip_on=True)
