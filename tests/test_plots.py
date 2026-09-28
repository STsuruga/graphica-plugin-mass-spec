"""パネルのグラフのマウス操作。matplotlib のイベントを合成して、offscreen の Qt で確かめる。"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from matplotlib.backend_bases import KeyEvent, MouseButton, MouseEvent  # noqa: E402

from mass_spec.plots import SpectrumPlot, TicPlot  # noqa: E402
from mass_spec.spectra import find_peaks  # noqa: E402
from synthetic import gaussian_peaks, tof_grid  # noqa: E402


def _px(plot, x, y=None):
    if y is None:
        y = sum(plot.ax.get_ylim()) / 2
    return plot.ax.transData.transform((x, y))


def _mouse(plot, name, x, y=None, button=MouseButton.LEFT, key=None, dblclick=False, step=0):
    px, py = _px(plot, x, y)
    if name == "scroll_event":
        from matplotlib.backend_bases import MouseEvent as ME
        ME("scroll_event", plot, px, py, button="up" if step > 0 else "down", key=key, step=step)._process()
        return
    MouseEvent(name, plot, px, py, button=button, key=key, dblclick=dblclick)._process()


def _drag(plot, x0, x1, key=None, button=MouseButton.LEFT, y0=None, y1=None):
    _mouse(plot, "button_press_event", x0, y0, button=button, key=key)
    _mouse(plot, "motion_notify_event", (x0 + x1) / 2, y0 if y1 is None else (y0 + y1) / 2, button=button, key=key)
    _mouse(plot, "motion_notify_event", x1, y1 if y1 is not None else y0, button=button, key=key)
    _mouse(plot, "button_release_event", x1, y1 if y1 is not None else y0, button=button, key=key)


@pytest.fixture
def tic():
    plot = TicPlot()
    plot.resize(600, 200)
    plot.set_data(np.arange(0, 60) / 60.0, np.r_[np.ones(20), np.full(20, 10.0), np.ones(20)])
    plot.draw()
    return plot


def test_tic_drag_selects_sample_and_shift_drag_selects_background(tic):
    changes = []
    tic.range_changed.connect(lambda kind, a, b: changes.append((kind, a, b)))
    _drag(tic, 0.4, 0.6)
    assert changes[-1][0] == "sample"
    assert changes[-1][1:] == pytest.approx((0.4, 0.6), abs=0.01)
    _drag(tic, 0.9, 0.8, key="shift")
    assert changes[-1][0] == "background"
    assert changes[-1][1:] == pytest.approx((0.8, 0.9), abs=0.01)


def test_tic_drag_edge_resizes_and_inside_moves(tic):
    tic.set_range("sample", 0.4, 0.6)
    tic.draw()
    _drag(tic, 0.6, 0.7)
    assert tic.ranges["sample"] == pytest.approx((0.4, 0.7), abs=0.01)
    _drag(tic, 0.5, 0.55)
    assert tic.ranges["sample"] == pytest.approx((0.45, 0.75), abs=0.01)


def test_tic_click_picks_a_scan_and_escape_cancels_a_drag(tic):
    clicked = []
    tic.scan_clicked.connect(clicked.append)
    _mouse(tic, "button_press_event", 0.5)
    _mouse(tic, "button_release_event", 0.5)
    assert clicked == [pytest.approx(0.5, abs=0.01)]
    tic.set_range("sample", 0.1, 0.2)
    _mouse(tic, "button_press_event", 0.5)
    _mouse(tic, "motion_notify_event", 0.7)
    KeyEvent("key_press_event", tic, "escape")._process()
    assert tic.ranges["sample"] == (0.1, 0.2)


def test_wheel_zoom_pan_back_and_reset(tic):
    home = tic.ax.get_xlim()
    _mouse(tic, "scroll_event", 0.5, step=1)
    zoomed = tic.ax.get_xlim()
    assert zoomed[1] - zoomed[0] < home[1] - home[0]
    _drag(tic, 0.5, 0.4, button=MouseButton.MIDDLE)
    panned = tic.ax.get_xlim()
    assert panned[0] > zoomed[0]
    KeyEvent("key_press_event", tic, "backspace")._process()
    assert tic.ax.get_xlim() == pytest.approx(zoomed)
    _mouse(tic, "button_press_event", 0.5, dblclick=True)
    assert tic.ax.get_xlim() == pytest.approx(home)


def _axis_px(plot, x, below=12):
    """横軸の目盛りの帯(グラフの枠のすぐ下)の画面座標。"""
    px = plot.ax.transData.transform((x, plot.ax.get_ylim()[0]))[0]
    return px, plot.ax.bbox.y0 - below


def test_wheel_on_the_x_axis_zooms_and_drag_on_it_pans(tic):
    home = tic.ax.get_xlim()
    px, py = _axis_px(tic, 0.5)
    MouseEvent("scroll_event", tic, px, py, button="up", step=1)._process()
    zoomed = tic.ax.get_xlim()
    assert zoomed[1] - zoomed[0] < home[1] - home[0]
    assert (zoomed[0] + zoomed[1]) / 2 == pytest.approx(0.5, abs=0.05)
    ylim = tic.ax.get_ylim()
    start, _ = _axis_px(tic, 0.5)
    end, _ = _axis_px(tic, 0.4)
    MouseEvent("button_press_event", tic, start, py, button=MouseButton.LEFT)._process()
    MouseEvent("motion_notify_event", tic, end, py, button=MouseButton.LEFT)._process()
    MouseEvent("button_release_event", tic, end, py, button=MouseButton.LEFT)._process()
    panned = tic.ax.get_xlim()
    assert panned[0] == pytest.approx(zoomed[0] + 0.1, abs=0.01)
    assert tic.ax.get_ylim() == pytest.approx(ylim)
    assert tic.ranges["sample"] is None  # 軸の上のドラッグは範囲を選ばない


def test_wheel_on_the_y_axis_zooms_vertically(tic):
    ylim = tic.ax.get_ylim()
    bbox = tic.ax.bbox
    MouseEvent("scroll_event", tic, bbox.x0 - 12, (bbox.y0 + bbox.y1) / 2, button="up", step=1)._process()
    assert tic.ax.get_ylim()[1] < ylim[1] and tic.ax.get_ylim()[0] == ylim[0]


def test_cursor_changes_over_the_axes(tic):
    from PySide6.QtCore import Qt
    px, py = _axis_px(tic, 0.5)
    MouseEvent("motion_notify_event", tic, px, py)._process()
    assert tic.cursor().shape() == Qt.CursorShape.SizeHorCursor
    MouseEvent("motion_notify_event", tic, *_px(tic, 0.5))._process()
    assert tic.cursor().shape() == Qt.CursorShape.ArrowCursor


@pytest.fixture
def spectrum():
    plot = SpectrumPlot()
    plot.resize(700, 300)
    grid = tof_grid(400.0, 600.0)
    centers = [450.0, 451.00335, 500.0, 550.0]
    y = gaussian_peaks(grid, centers, [100.0, 30.0, 60.0, 5.0], 12000)
    plot.set_spectrum(grid, y, find_peaks(grid, y, min_height=0.5))
    plot.draw()
    return plot


def test_spectrum_labels_follow_the_view(spectrum):
    spectrum.label_top_n = 2
    spectrum.update_labels()
    assert spectrum.labeled_texts() == ["450.0000", "500.0000"]
    _drag(spectrum, 540.0, 560.0, y0=10.0)
    assert spectrum.ax.get_xlim() == pytest.approx((540.0, 560.0), abs=0.2)
    assert spectrum.ax.get_ylim()[1] == pytest.approx(5.0 * 1.15, rel=0.05)
    assert spectrum.labeled_texts() == ["550.0000"]
    spectrum.label_format.mz_decimals = 2
    spectrum.update_labels()
    assert spectrum.labeled_texts() == ["550.00"]


def test_spectrum_box_zoom_and_click_to_pin(spectrum):
    _drag(spectrum, 440.0, 520.0, key="control", y0=5.0, y1=80.0)
    assert spectrum.ax.get_ylim() == pytest.approx((5.0, 80.0), abs=1.0)
    spectrum.reset_view()
    spectrum.label_top_n = 1
    spectrum.update_labels()
    _mouse(spectrum, "button_press_event", 550.0, 3.0)
    _mouse(spectrum, "button_release_event", 550.0, 3.0)
    assert spectrum.pinned == [pytest.approx(550.0, abs=0.01)]
    assert "550.0000" in spectrum.labeled_texts()
    _mouse(spectrum, "button_press_event", 550.0, 3.0)
    _mouse(spectrum, "button_release_event", 550.0, 3.0)
    assert spectrum.pinned == []


def test_spectrum_shift_drag_measures_isotope_spacing(spectrum):
    spectrum.set_view((448.0, 453.0), (0, 110))
    spectrum.draw()
    measured = []
    spectrum.measured.connect(measured.append)
    _drag(spectrum, 450.02, 450.98, key="shift", y0=50.0)
    assert "Δm/z 1.003" in measured[-1] and "z = 1" in measured[-1]


def test_spectrum_hover_reports_the_nearest_peak(spectrum):
    texts = []
    spectrum.hover_text.connect(texts.append)
    _mouse(spectrum, "motion_notify_event", 500.0, 10.0, button=None)
    assert texts[-1].startswith("m/z 500.0000  強度 60")


def test_overlays_are_drawn_and_included_in_the_height(spectrum):
    spectrum.set_overlays([(np.array([575.0]), np.array([500.0]), "#ff0000", "calc", "stick")])
    spectrum.set_view((570.0, 580.0), (0, spectrum._top(570.0, 580.0)))
    assert spectrum.ax.get_ylim()[1] == pytest.approx(575.0, rel=0.01)


def test_y_axis_uses_superscript_powers_from_1e4():
    plot = SpectrumPlot()
    plot.resize(600, 250)
    grid = tof_grid(400.0, 600.0)
    plot.set_spectrum(grid, gaussian_peaks(grid, [500.0], [2.0e5], 12000), [])
    plot.draw()
    text = plot.ax.yaxis.get_offset_text().get_text()
    assert "times" in text and "10^{5}" in text
    plot.set_spectrum(grid, gaussian_peaks(grid, [500.0], [5000.0], 12000), [])
    plot.draw()
    assert plot.ax.yaxis.get_offset_text().get_text() == ""
