import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from conftest import GRAPHICA_AVAILABLE, requires_graphica  # noqa: E402
from synthetic import ion_spectrum, run_scans, tof_grid, write_mzml  # noqa: E402

if GRAPHICA_AVAILABLE:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    from graphica.plugin import AnalysisResult, Dataset
    from graphica.plugin.testing import FakeGraphicaPluginAPI, load_plugin_like_graphica

PLUGIN_DIR = os.path.join(ROOT, "mass_spec")


def _spectrum_dataset(formula="C6H12O6", adduct="[M+Na]+"):
    grid = tof_grid(150.0, 450.0)
    y, _ = ion_spectrum(grid, formula, adduct, shift_ppm=2.0)
    return Dataset(name="実測", df=pd.DataFrame({"m/z": grid, "強度": y}), x_col_name="m/z", y_col_name="強度")


@requires_graphica
def test_register_adds_the_hooks():
    from mass_spec import register
    api = FakeGraphicaPluginAPI(plugin_name="mass_spec")
    register(api)
    assert "Stick" in api.plot_types
    assert "同位体パターンと照合" in api.analyzers
    assert ".mzml" in api.importers  # 本体は拡張子を小文字にそろえて照合する


@requires_graphica
def test_loads_like_graphica(tmp_path):
    api, record = load_plugin_like_graphica(PLUGIN_DIR, work_dir=str(tmp_path))
    assert record["error"] is None
    assert api.get_plot_type("Stick") is not None
    assert api.get_importer_for_extension(".MZML") is not None
    # パネルの import は遅延なので、本番の経路で作って相対 import の誤りを捕まえる
    from graphica.plugin.testing import FakePluginContext
    panel = next(p for p in api.get_panels() if p.name == "MS スペクトル")
    widget = panel.widget_factory(FakePluginContext(data_dir=str(tmp_path / "data")))
    assert widget.spectrum_plot is not None
    widget.close()


@requires_graphica
def test_analyzer_matches_and_overlays():
    from mass_spec.analyzer import PRESET_POSITIVE, analyze
    ds = _spectrum_dataset()
    before = ds.df.copy()
    result = analyze(ds, {"formula": "C6H12O6", "preset": PRESET_POSITIVE, "extra_adducts": "[2M+Na]+",
                          "resolution": 0.0, "tolerance_ppm": 50.0, "detect_percent": 0.5,
                          "mz_decimals": 3, "ppm_decimals": 1, "overlay": True})
    assert isinstance(result, AnalysisResult)
    pd.testing.assert_frame_equal(ds.df, before)
    table = result.table
    assert table.loc[0, "付加イオン"] == "M(中性分子)"
    na = table[(table["付加イオン"] == "[M+Na]+") & (table["ピーク"] == "M")].iloc[0]
    assert na["状態"] == "検出"
    assert na["誤差 (ppm)"] == pytest.approx(2.0, abs=1.0)
    assert set(table[table["付加イオン"] == "[M+H]+"]["状態"]) == {"未検出"}
    assert len(result.new_datasets) == 2
    profile, sticks = result.new_datasets
    assert sticks.plot_type == "Stick" and sticks.show_point_labels
    assert sticks.df["ラベル"].iloc[0] == f"{sticks.df['m/z'].iloc[0]:.3f}"
    assert all(d.color.startswith("#") and len(d.color) == 7 for d in result.new_datasets)
    assert result.annotations[0]["text"].startswith("[M+Na]+ ")
    assert "ppm" in result.annotations[0]["text"]


@requires_graphica
@pytest.mark.parametrize("params, message", [
    ({"formula": "C6H12Q"}, "元素記号"),
    ({"formula": "C6H12O6", "preset": "下の欄だけ", "extra_adducts": ""}, "付加イオン"),
])
def test_analyzer_reports_input_errors(params, message):
    from mass_spec.analyzer import analyze
    with pytest.raises(ValueError, match=message):
        analyze(_spectrum_dataset(), params)


@requires_graphica
def test_importer_returns_tic_or_single_spectrum(tmp_path):
    from mass_spec.analyzer import load_mzml_file
    many = load_mzml_file(write_mzml(tmp_path / "many.mzML", run_scans(n=5)))
    assert list(many.columns) == ["時間 (min)", "TIC"] and len(many) == 5
    one = load_mzml_file(write_mzml(tmp_path / "one.mzML", run_scans(n=1)))
    assert list(one.columns) == ["m/z", "強度"]
    with pytest.raises(ValueError):
        bad = tmp_path / "bad.mzML"
        bad.write_text("<mzML>", encoding="utf-8")
        load_mzml_file(str(bad))


@requires_graphica
def test_stick_drawer_draws_bars_or_only_a_legend_entry():
    from mass_spec.plot_types import draw_stick
    ds = Dataset(name="s", df=pd.DataFrame({"x": [1.0, 2.0], "y": [3.0, np.nan]}), x_col_name="x", y_col_name="y")
    fig, ax = plt.subplots()
    artist = draw_stick(ds, ax, ds.x_data, ds.y_data)
    assert len(artist.get_segments()) == 1
    assert artist.get_label() == "s"
    ds.linewidth = 0.0
    hidden = draw_stick(ds, ax, ds.x_data, ds.y_data)
    assert hidden.get_linestyle() == "None"
    plt.close(fig)
