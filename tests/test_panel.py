import os
import sys
import time

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from conftest import GRAPHICA_AVAILABLE, requires_graphica  # noqa: E402
from synthetic import mzml_text, run_scans, write_mzml  # noqa: E402

if GRAPHICA_AVAILABLE:
    from PySide6.QtWidgets import QApplication

    from graphica.plugin.testing import FakePluginContext

pytestmark = requires_graphica

FAKE_MSCONVERT = '''\
import os, sys
args = sys.argv[1:]
out = os.path.join(args[args.index("-o") + 1], args[args.index("--outfile") + 1])
here = os.path.dirname(os.path.abspath(__file__))
if os.path.exists(os.path.join(here, "fail")):
    print("fake failure"); sys.exit(3)
with open(os.path.join(here, "calls.txt"), "a") as f:
    f.write(" ".join(args[1:]) + "\\n")
with open(os.path.join(here, "payload.mzML"), encoding="utf-8") as f:
    text = f.read()
with open(out, "w", encoding="utf-8") as f:
    f.write(text)
'''


def _wait(condition, timeout=20.0):
    end = time.monotonic() + timeout
    while not condition():
        QApplication.processEvents()
        if time.monotonic() > end:
            raise AssertionError("時間内に終わりませんでした")
        time.sleep(0.02)


@pytest.fixture
def ctx(tmp_path):
    return FakePluginContext(plugin_name="mass_spec", data_dir=str(tmp_path / "data"),
                             color_cycle=["#1f77b4", "#ff7f0e", "#2ca02c"])


@pytest.fixture
def panel(ctx):
    from mass_spec.panel import MassSpecPanel
    p = MassSpecPanel(ctx)
    p.resize(800, 1100)
    yield p
    p.close()
    p.deleteLater()


@pytest.fixture
def mzml(tmp_path):
    return write_mzml(tmp_path / "sample.mzML", run_scans(n=10))


def _calculate(panel, formula="C6H12O6", preset=0, extra=""):
    w = panel.calc_window
    w.formula_edit.setText(formula)
    w.preset_combo.setCurrentIndex(preset)
    w.extra_edit.setText(extra)
    panel.run_calculation()


def test_menu_bar_and_toolbar(panel):
    titles = [a.text() for a in panel.menu_bar.actions()]
    assert titles == ["ファイル", "表示", "解析", "転送", "ヘルプ"]
    assert panel.pane_count() == 3 and panel.pane_spin.value() == 3
    assert panel.splitter.count() == 3


def test_open_shows_tic_and_the_average_of_all_scans(panel, mzml):
    run = panel.open_mzml(mzml)
    assert run is not None and panel.run_combo.count() == 1
    assert "10 スキャン" in panel.run_combo.currentText()
    assert panel.range_window.sample_from.value() == pytest.approx(0.0)
    assert panel.range_window.sample_to.value() == pytest.approx(9 / 60, abs=1e-3)
    assert panel.spectrum is not None and panel.spectrum["provenance"]["scans"] == 10


def test_sample_and_background_ranges_subtract(panel, mzml):
    panel.open_mzml(mzml)
    panel._on_tic_range_changed("sample", 3 / 60, 6 / 60)
    with_bg_peak = panel.spectrum["y"].max()
    panel._on_tic_range_changed("background", 0.0, 2 / 60)
    s = panel.spectrum
    assert "背景" in s["name"] and s["provenance"]["background_scans"] == 3
    near_300 = np.abs(s["mz"] - 300.0) < 0.05
    assert s["y"][near_300].max() == pytest.approx(0.0, abs=1e-6)
    assert s["y"].max() == pytest.approx(with_bg_peak, rel=0.01)
    assert panel.range_window.bg_from.value() == pytest.approx(0.0) and panel.range_window.bg_to.value() == pytest.approx(2 / 60, abs=1e-3)
    panel.settings_window.subtract_check.setChecked(False)
    assert "背景" not in panel.spectrum["name"]


def test_clicking_a_scan_selects_that_scan(panel, mzml):
    panel.open_mzml(mzml)
    panel._on_scan_clicked(4.2 / 60)
    assert panel.spectrum["provenance"]["scans"] == 1
    assert panel.range_window.sample_from.value() == pytest.approx(4 / 60, abs=1e-3)


def test_calculation_window_matches_and_overlays(panel, mzml):
    panel.open_mzml(mzml)
    panel._on_tic_range_changed("sample", 3 / 60, 6 / 60)
    _calculate(panel)
    header, rows = panel.calc_rows()
    na = [r for r in rows if r[0] == "[M+Na]+" and r[1] == "M"][0]
    assert na[-1] == "検出" and na[3]
    assert {r[-1] for r in rows if r[0] == "[M+H]+"} == {"未検出"}
    assert panel.calc_window.table.rowCount() == len(rows)
    assert len(panel.spectrum_plot._overlays) == 1
    assert "C6H12O6" in panel.calc_window.summary.text()
    panel.calc_window.overlay_check.setChecked(False)
    assert panel.spectrum_plot._overlays == []


def test_calculation_without_a_measurement(panel, ctx):
    _calculate(panel)
    header, rows = panel.calc_rows()
    assert rows and all(r[-1] == "計算のみ" for r in rows)
    panel.transfer_calculation()
    assert all(d.plot_type == "Stick" for d in ctx.datasets())
    assert len(ctx.datasets()) == 5


def test_formula_errors_are_reported(panel, ctx):
    _calculate(panel, formula="C6H12Qq")
    assert ctx.messages[-1][0] == "error" and "元素記号" in ctx.messages[-1][2]


def test_spectrum_copies_into_panes(panel, mzml):
    panel.open_mzml(mzml)
    panel.copy_spectrum_to_pane()
    panel._on_scan_clicked(4 / 60)
    panel.copy_spectrum_to_pane()
    assert [p.item["kind"] for p in panel.panes[:2]] == ["measured", "measured"]
    assert panel.panes[2].item is None
    assert "10" not in panel.panes[1].item["name"] or panel.panes[1].item["provenance"]["scans"] == 1
    panel.copy_spectrum_to_pane(2)
    panel.copy_spectrum_to_pane()
    assert panel.first_empty_pane() is None


def test_found_patterns_fill_panes_and_one_item_per_pane(panel, mzml):
    panel.open_mzml(mzml)
    panel._on_tic_range_changed("sample", 3 / 60, 6 / 60)
    _calculate(panel, extra="[M+H]+")
    panel.found_patterns_to_panes()
    item = panel.panes[0].item
    assert item["kind"] == "calc" and "[M+Na]+" in item["name"]
    assert panel.panes[1].item is None
    panel.calc_window.table.selectRow(0)  # [M+H]+ の行(未検出でも計算パターンは出せる)
    panel.calc_window.pane_combo.setCurrentIndex(1)
    panel.calc_window._send_selected_to_pane()
    assert "[M+H]+" in panel.panes[1].item["name"]


def test_pane_count_spin(panel, mzml):
    panel.open_mzml(mzml)
    panel.copy_spectrum_to_pane(0)
    panel.pane_spin.setValue(5)
    assert panel.pane_count() == 5 and panel.calc_window.pane_combo.count() == 5
    panel.pane_spin.setValue(1)
    assert panel.pane_count() == 1 and panel.panes[0].item is not None
    assert panel.settings["pane_count"] == 1


def test_sync_mode_links_only_the_panes(panel, mzml):
    panel.open_mzml(mzml)
    panel.copy_spectrum_to_pane(0)
    panel.copy_spectrum_to_pane(1)
    panel.panes[0].set_view((200.0, 250.0))
    assert panel.panes[1].ax.get_xlim() == pytest.approx((200.0, 250.0))
    assert panel.spectrum_plot.ax.get_xlim() != pytest.approx((200.0, 250.0))
    panel.sync_check.setChecked(False)
    panel.panes[0].set_view((210.0, 220.0))
    assert panel.panes[1].ax.get_xlim() == pytest.approx((200.0, 250.0))
    panel.sync_check.setChecked(True)
    assert panel.panes[1].ax.get_xlim() == pytest.approx((210.0, 220.0))


def test_a_new_pane_item_follows_the_synced_range(panel, mzml):
    panel.open_mzml(mzml)
    panel.copy_spectrum_to_pane(0)
    panel.panes[0].set_view((200.0, 250.0))
    panel.copy_spectrum_to_pane(1)
    assert panel.panes[1].ax.get_xlim() == pytest.approx((200.0, 250.0))


def test_transfers_add_datasets_with_labels_colors_and_subplot(panel, mzml, ctx):
    panel.open_mzml(mzml)
    panel._on_tic_range_changed("sample", 3 / 60, 6 / 60)
    panel.settings_window.subplot_spin.setValue(2)
    panel.settings_window.mz_decimals_spin.setValue(2)
    panel.transfer_tic()
    panel.transfer_spectrum()
    panel.transfer_centroid()
    panel.transfer_centroid(labels_only=True)
    _calculate(panel)
    panel.transfer_calculation()
    datasets = ctx.datasets()
    assert [d.plot_type for d in datasets] == ["Line", "Line", "Stick", "Stick", "Line", "Stick"]
    assert all(d.subplot_target == 1 for d in datasets)
    assert all(d.color.startswith("#") and len(d.color) == 7 for d in datasets)
    tic, spectrum, centroid, labels, calc_profile, calc_sticks = datasets
    assert list(tic.df.columns) == ["時間 (min)", "TIC"]
    assert len(spectrum.df) < len(panel.spectrum["mz"])  # 0 の区間を詰めている
    labeled = [t for t in centroid.df["ラベル"] if t]
    assert labeled and all(len(t.split(".")[1]) == 2 for t in labeled)
    assert labels.linewidth == 0 and labels.show_point_labels
    assert "[M+Na]+" in calc_sticks.name
    assert ctx.undo_descriptions[0].startswith("[MS]")


def test_pane_contents_can_be_transferred(panel, mzml, ctx):
    panel.open_mzml(mzml)
    panel._on_tic_range_changed("sample", 3 / 60, 6 / 60)
    panel.copy_spectrum_to_pane(0)
    _calculate(panel)
    panel.found_patterns_to_panes()  # 枠 1 を計算パターンで置き換える
    panel.copy_spectrum_to_pane(1)
    panel.transfer_pane(panel.panes[0])
    panel.transfer_pane(panel.panes[1])
    names = [d.name for d in ctx.datasets()]
    assert any("[M+Na]+" in n for n in names)
    measured = ctx.datasets()[-1].df["m/z"]
    assert measured.min() < 151 and measured.max() > 449  # 実測は全範囲を送る


def test_view_only_transfer_crops_the_spectrum(panel, mzml, ctx):
    panel.open_mzml(mzml)
    panel.spectrum_plot.set_view((200.0, 250.0), (0, 100))
    panel.settings_window.view_only_check.setChecked(True)
    panel.transfer_spectrum()
    mz = ctx.datasets()[-1].df["m/z"]
    assert mz.min() >= 200.0 and mz.max() <= 250.0


def test_settings_persist_between_panels(ctx, mzml):
    from mass_spec.panel import MassSpecPanel
    first = MassSpecPanel(ctx)
    first.settings_window.mz_decimals_spin.setValue(3)
    first.settings_window.label_mode_combo.setCurrentIndex(1)
    first.settings_window.view_only_check.setChecked(True)
    first.pane_spin.setValue(4)
    first.sync_check.setChecked(False)
    first.close()
    second = MassSpecPanel(ctx)
    assert second.settings_window.mz_decimals_spin.value() == 3
    assert second.settings["label_mode"] == "percent"
    assert second.spectrum_plot.label_format.mz_decimals == 3
    assert second.panes[0].label_format.mz_decimals == 3
    assert second.settings_window.view_only_check.isChecked()
    assert second.pane_count() == 4 and not second.sync_check.isChecked()
    second.close()


def _fake_d(tmp_path, fail=False):
    d = tmp_path / "measure.d"
    d.mkdir()
    (d / "analysis.baf").write_bytes(b"raw")
    (d / "__main__.py").write_text(FAKE_MSCONVERT, encoding="utf-8")
    (d / "payload.mzML").write_text(mzml_text(run_scans(n=4)), encoding="utf-8")
    if fail:
        (d / "fail").write_text("")
    return d


def test_open_d_converts_once_then_uses_the_cache(panel, ctx, tmp_path):
    # Python にフォルダを渡すと __main__.py が動くので、偽の msconvert として使う
    panel.settings["msconvert_path"] = sys.executable
    d = _fake_d(tmp_path)
    panel.open_d(str(d))
    assert panel.cancel_button.isVisibleTo(panel) and panel.cancel_action.isEnabled()
    _wait(lambda: panel._process is None)
    assert panel.run_combo.count() == 1
    assert panel.run_combo.currentText().startswith("measure")
    calls = (d / "calls.txt").read_text().splitlines()
    assert len(calls) == 1 and "--zlib" in calls[0]
    cache = os.path.join(ctx.data_dir, "mzml_cache")
    assert [f for f in os.listdir(cache) if f.endswith(".mzML")]
    assert not [f for f in os.listdir(cache) if f.endswith(".part")]
    panel.open_d(str(d))
    assert panel._process is None and panel.run_combo.count() == 2
    assert len((d / "calls.txt").read_text().splitlines()) == 1


def test_failed_conversion_reports_the_output(panel, ctx, tmp_path):
    panel.settings["msconvert_path"] = sys.executable
    panel.open_d(str(_fake_d(tmp_path, fail=True)))
    _wait(lambda: panel._process is None)
    kind, _title, text = ctx.messages[-1]
    assert kind == "error" and "終了コード 3" in text and "fake failure" in text
    cache = os.path.join(ctx.data_dir, "mzml_cache")
    assert os.listdir(cache) == []


def test_open_d_rejects_other_folders_and_missing_msconvert(panel, ctx, tmp_path, monkeypatch):
    panel.open_d(str(tmp_path))
    assert "フォルダではありません" in ctx.messages[-1][2]
    from mass_spec import convert
    monkeypatch.setattr(convert, "find_msconvert", lambda configured=None: None)
    asked = []
    monkeypatch.setattr(panel, "ask_msconvert", lambda explain=True: asked.append(True))
    panel.open_d(str(_fake_d(tmp_path)))
    assert asked and panel._process is None


def test_broken_mzml_is_reported(panel, ctx, tmp_path):
    bad = tmp_path / "bad.mzML"
    bad.write_text("<mzML>", encoding="utf-8")
    assert panel.open_mzml(str(bad)) is None
    assert ctx.messages[-1][0] == "error"


def test_register_adds_panel_and_help(tmp_path):
    from graphica.plugin.testing import FakeGraphicaPluginAPI

    from mass_spec import HELP_MENU, PANEL_NAME, register
    api = FakeGraphicaPluginAPI(plugin_name="mass_spec")
    register(api)
    assert PANEL_NAME in api.panels
    ctx = FakePluginContext(plugin_name="mass_spec", data_dir=str(tmp_path))
    action = next(a for a in api.menu_actions if a.text == HELP_MENU)
    action.callback(ctx)
    assert "MS スペクトル" in ctx.messages[-1][2]
    widget = api.panels[PANEL_NAME]["widget_factory"](ctx)
    assert widget.run_combo.count() == 0
    widget.close()


def test_range_fields_live_in_the_detail_window(panel, mzml):
    panel.open_mzml(mzml)
    w = panel.range_window
    assert not w.isVisible() and w.parent() is panel
    w.sample_from.setValue(3 / 60)
    w.sample_to.setValue(6 / 60)
    panel.compute_spectrum()
    assert panel.spectrum["provenance"]["scans"] == 4
    assert panel.tic_plot.ranges["sample"] == pytest.approx((3 / 60, 6 / 60))
    panel.reset_sample_range()
    assert panel.spectrum["provenance"]["scans"] == 10
    assert w.sample_to.value() == pytest.approx(9 / 60, abs=1e-3)
