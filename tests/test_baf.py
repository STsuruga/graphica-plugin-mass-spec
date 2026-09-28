import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mass_spec.baf import BafUnsupported, decode_profile, is_baf_folder, mz_axis, read_baf  # noqa: E402
from synthetic_baf import CALIBRATION, encode_profile, peaks_profile, write_baf_d  # noqa: E402


@pytest.mark.parametrize("values", [
    [0, 0, 0, 4, 0, 0],
    [7] + [0] * 500 + [3, 9, 15, 13, 6, 4] + [0] * 20,
    [0, 48, 0, 110, 158, 165, 141, 102, 65, 40, 0],       # 中くらいの差(8・10 ビットの符号)
    [0, 5282, 4132, 1801, 928, 250000, 3, 0],            # 大きな差(値そのもの)
    list(range(0, 3000, 7)) + list(range(3000, 0, -13)),
])
def test_profile_round_trip(values):
    assert np.array_equal(decode_profile(encode_profile(values)), values)


def test_profile_round_trip_random():
    rng = np.random.default_rng(1)
    y = peaks_profile(20000)
    y[rng.integers(0, len(y), 50)] = rng.integers(0, 3_000_000, 50)
    assert np.array_equal(decode_profile(encode_profile(y)), y)


def test_decode_rejects_broken_payloads():
    good = encode_profile([0, 4, 0, 0, 9, 0])
    with pytest.raises(BafUnsupported):
        decode_profile(good[:12])
    with pytest.raises(BafUnsupported):
        decode_profile(good[:4] + (10 ** 6).to_bytes(4, "little") + good[8:])   # 点数が足りない
    with pytest.raises(BafUnsupported):
        decode_profile(good[:8] + (2).to_bytes(4, "little") + good[12:])


def test_mz_axis_solves_the_tof_equation():
    delay, interval, c1, c2, c3, offset = CALIBRATION
    mz = mz_axis(CALIBRATION, 5000)
    assert np.all(np.diff(mz) > 0)
    x = np.sqrt(mz + offset)
    tof = delay + interval * np.arange(5000)
    residual = c3 * x * x + np.sqrt(1e12 / c1) * x + (c2 - tof)
    assert np.max(np.abs(residual)) < 1e-6 * tof.max()
    linear = mz_axis((delay, interval, c1, c2, 0.0, offset), 10)
    np.testing.assert_allclose(np.sqrt(linear + offset), (tof[:10] - c2) / np.sqrt(1e12 / c1))


def test_read_baf(tmp_path):
    profiles = [peaks_profile(3000, noise_seed=s) for s in (0, 1, 2)]
    scans = [(1100, -1, profiles[0], 400.0), (2100, -1, profiles[1], 500.0), (3100, -1, profiles[2], 600.0)]
    d = write_baf_d(tmp_path / "sample neg.d", scans)
    assert is_baf_folder(d) and not is_baf_folder(str(tmp_path))
    run = read_baf(d)
    assert run.name == "sample neg"
    assert run.polarity() == -1
    np.testing.assert_allclose(run.times(), [1.1 / 60, 2.1 / 60, 3.1 / 60])
    np.testing.assert_allclose(run.tics(), [400.0, 500.0, 600.0])
    mz, y = run.scans[1].arrays()
    assert np.array_equal(y, profiles[1])
    np.testing.assert_allclose(mz, mz_axis(CALIBRATION, 3000))
    assert run.scans[0].centroid is False and run.scans[0].ms_level == 1


def test_read_baf_positive_and_scans_in_range(tmp_path):
    d = write_baf_d(tmp_path / "p.d", [(t, 1, peaks_profile(3000), 10.0) for t in (1000, 2000, 3000, 4000)])
    run = read_baf(d)
    assert run.polarity() == 1
    assert len(run.scans_in_range(1.5 / 60, 3.5 / 60)) == 2


def test_read_baf_rejects_other_formats(tmp_path):
    with pytest.raises(BafUnsupported):
        read_baf(str(tmp_path))
    d = write_baf_d(tmp_path / "other.d", [(1000, 1, [0, 1, 0], 1.0)], declare=False)
    with pytest.raises(BafUnsupported, match="圧縮方式"):
        read_baf(d)
