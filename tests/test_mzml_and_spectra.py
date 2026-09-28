import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mass_spec.mzml import MzmlError, read_mzml  # noqa: E402
from mass_spec.spectra import (  # noqa: E402
    LabelFormat, Peak, average_spectrum, compact_zeros, find_peaks, is_centroid_like,
    select_label_peaks, subtract_background,
)
from synthetic import gaussian_peaks, mzml_text, run_scans, tof_grid, write_mzml  # noqa: E402


@pytest.mark.parametrize("options", [
    {},
    {"mz_dtype": "<f4", "int_dtype": "<f8"},
    {"compress": False},
    {"indexed": True},
    {"param_group": True},
    {"time_unit": "minute"},
])
def test_read_mzml_variants(tmp_path, options):
    scans = run_scans(n=4)
    if options.get("time_unit") == "minute":
        scans = [(t / 60.0, mz, y) for t, mz, y in scans]
    path = write_mzml(tmp_path / "a.mzML", scans, **options)
    run = read_mzml(path)
    assert run.name == "a"
    assert len(run.scans) == 4
    np.testing.assert_allclose(run.times(), [0, 1 / 60, 2 / 60, 3 / 60], atol=1e-9)
    assert run.polarity() == 1
    mz, y = run.scans[2].arrays()
    rtol = 1e-6 if options.get("mz_dtype") == "<f4" else 1e-12
    np.testing.assert_allclose(mz, scans[2][1], rtol=rtol)
    np.testing.assert_allclose(y, scans[2][2], rtol=1e-6)
    assert run.scans[0].centroid is False


def test_read_mzml_polarity_and_computed_tic(tmp_path):
    scans = run_scans(n=3)
    run = read_mzml(write_mzml(tmp_path / "n.mzML", scans, polarity="negative", with_tic=False, profile=False))
    assert run.polarity() == -1
    assert run.scans[0].centroid is True
    assert run.tics()[0] == pytest.approx(np.sum(scans[0][2]), rel=1e-5)


def test_read_mzml_rejects_numpress_and_garbage(tmp_path):
    text = mzml_text(run_scans(n=1)).replace("MS:1000574", "MS:1002312", 1)
    path = tmp_path / "np.mzML"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(MzmlError, match="numpress"):
        read_mzml(str(path))
    bad = tmp_path / "bad.mzML"
    bad.write_text("<mzML><run>", encoding="utf-8")
    with pytest.raises(MzmlError):
        read_mzml(str(bad))


def test_scans_in_range_and_progress(tmp_path):
    calls = []
    run = read_mzml(write_mzml(tmp_path / "r.mzML", run_scans(n=120)), progress=lambda a, b: calls.append((a, b)))
    assert calls[-1] == (120, 120)
    assert [s.index for s in run.scans_in_range(3 / 60, 6 / 60)] == [3, 4, 5, 6]
    assert run.nearest_scan(4.4 / 60).index == 4


def test_average_spectrum_by_index_despite_small_grid_shifts(tmp_path):
    run = read_mzml(write_mzml(tmp_path / "r.mzML", run_scans(n=10)))
    mz, y = average_spectrum(run.scans_in_range(3 / 60, 6 / 60))
    grid = tof_grid(150.0, 450.0)
    assert len(mz) == len(grid)
    np.testing.assert_allclose(mz, grid, rtol=1e-6)
    _bg_mz, bg_y = average_spectrum(run.scans_in_range(0, 2 / 60))
    assert y.max() > bg_y.max() * 5


def test_average_spectrum_interpolates_other_grids():
    class S:
        def __init__(self, mz, y):
            self._a = (np.asarray(mz, float), np.asarray(y, float))

        def arrays(self):
            return self._a

    mz, y = average_spectrum([S([1, 2, 3], [0, 2, 0]), S([1.5, 2.5], [4, 4])])
    np.testing.assert_allclose(mz, [1, 2, 3])
    np.testing.assert_allclose(y, [0, 3, 0])
    with pytest.raises(ValueError):
        average_spectrum([])


def test_subtract_background_and_clip():
    mz = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(subtract_background(mz, [5, 1, 3], mz, [1, 2, 1]), [4, 0, 2])
    np.testing.assert_allclose(subtract_background(mz, [5, 1, 3], mz, [1, 2, 1], clip_negative=False), [4, -1, 2])
    np.testing.assert_allclose(subtract_background(mz, [5, 1, 3], [0.5, 1.5], [2, 2]), [3, 1, 3])  # 背景の範囲外は 0


def test_compact_zeros_keeps_the_shape():
    mz = np.arange(10.0)
    y = np.array([0, 0, 0, 1, 2, 0, 0, 0, 0, 0.0])
    cx, cy = compact_zeros(mz, y)
    np.testing.assert_array_equal(cx, [0, 2, 3, 4, 5, 9])
    np.testing.assert_allclose(np.interp(mz, cx, cy), y)


def test_find_peaks_on_a_coarse_tof_grid():
    grid = tof_grid(400.0, 600.0)
    centers = [450.123456, 500.2, 550.3]
    y = gaussian_peaks(grid, centers, [100.0, 50.0, 10.0], 12000)
    peaks = find_peaks(grid, y, min_height=1.0)
    assert len(peaks) == 3
    for p, c in zip(peaks, centers):
        assert (p.mz - c) / c * 1e6 == pytest.approx(0.0, abs=1.0)
        assert p.resolution == pytest.approx(12000, rel=0.15)
    assert peaks[0].height == pytest.approx(100.0, rel=0.02)


def test_find_peaks_treats_sparse_data_as_centroid():
    mz = np.array([100.0, 101.0, 102.0, 150.0])
    y = np.array([10.0, 5.0, 0.0, 3.0])
    assert is_centroid_like(mz, y)
    peaks = find_peaks(mz, y)
    assert [p.mz for p in peaks] == [100.0, 101.0, 150.0]


def test_select_label_peaks():
    peaks = [Peak(100.0, 10, 0.01, 0), Peak(100.02, 9, 0.01, 1), Peak(200.0, 5, 0.01, 2),
             Peak(300.0, 1, 0.01, 3), Peak(400.0, 0.1, 0.01, 4)]
    assert [p.mz for p in select_label_peaks(peaks, top_n=3)] == [100.0, 200.0, 300.0]
    assert [p.mz for p in select_label_peaks(peaks, top_n=None, min_relative=20)] == [100.0, 200.0]
    assert [p.mz for p in select_label_peaks(peaks, x_range=(150, 450), top_n=1)] == [200.0]
    assert [p.mz for p in select_label_peaks(peaks, top_n=1, pinned=[399.9])] == [100.0, 400.0]


def test_label_format():
    assert LabelFormat(mz_decimals=2).label(523.24512) == "523.25"
    assert LabelFormat(mz_decimals=4, intensity="relative").label(523.24512, 50, 200) == "523.2451 (25.0%)"
    assert LabelFormat(mz_decimals=1, intensity="absolute", intensity_decimals=2).label(1.0, 12345, 1) == "1.0 (1.23e+04)"
    assert LabelFormat(ppm_decimals=2).ppm(-1.234) == "-1.23 ppm"


def test_isotope_clusters_group_signals():
    from mass_spec.spectra import cluster_heads, isotope_clusters
    peaks = [Peak(450.0, 100, 0.05, 0), Peak(451.0034, 30, 0.05, 1), Peak(452.0068, 5, 0.05, 2),
             Peak(500.0, 60, 0.05, 3), Peak(500.30, 50, 0.05, 4),             # 0.3 離れは別のシグナル
             Peak(600.0, 40, 0.05, 5), Peak(600.5017, 45, 0.05, 6),           # 2 価の同位体
             Peak(700.0, 50, 0.05, 7), Peak(701.998, 49, 0.05, 8)]            # Br の M+2
    groups = [[p.mz for p in c] for c in isotope_clusters(peaks)]
    assert groups == [[450.0, 451.0034, 452.0068], [500.0], [500.30], [600.0, 600.5017], [700.0, 701.998]]
    assert [p.mz for p in cluster_heads(peaks)] == [450.0, 500.0, 500.30, 600.5017, 700.0]
