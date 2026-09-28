"""テスト用の合成 BAF(mass_spec.baf が読む形)。ユーザーの測定データは使わない。"""
import os
import struct

import numpy as np

from mass_spec.baf import COMPRESSOR_ID, COMPRESSOR_NAME

CALIBRATION = (10000.0, 0.5, 400000.0, 250.0, 1e-4, -0.01)   # 遅れ, 間隔, c1, c2, c3, ずれ


def _block(kind_hex, header_extra=b"", body=b""):
    hsize = 12 + len(header_extra)
    return struct.pack("<I4sI", hsize + len(body), bytes.fromhex(kind_hex), hsize) + header_extra + body


def encode_profile(values):
    """decode_profile の逆。values は 0 以上の整数(最初の値は 255 以下)。"""
    values = [int(v) for v in values]
    codes = []
    prev = values[0]
    for v in values[1:]:
        d = v - prev
        z = 2 * abs(d) + (1 if d < 0 else 0)
        code = z + 32
        k = code.bit_length() - 6
        codes.append("0" * k + format(code, "b") if k <= 6 else "00000001" + format(v, "032b"))
        prev = v
    bits = "".join(codes)
    bits += "0" * (-len(bits) % 8)
    stream = int(bits, 2).to_bytes(len(bits) // 8, "big") if bits else b""
    return struct.pack("<III", 0, len(values), 1) + bytes([values[0]]) + stream


def scan_header(time_ms, polarity, calibration=CALIBRATION):
    extra = bytearray(28)                      # 見出しの 12〜39 バイト目
    struct.pack_into("<I", extra, 24 - 12, time_ms)
    extra[32 - 12] = 0 if polarity > 0 else 1
    delay, interval, c1, c2, c3, offset = calibration
    doubles = struct.pack("<10d", delay, interval, c2, c1, c3, 0.0, offset, 0.0, 0.0, 0.0)
    sub = struct.pack("<I4sII", 16 + len(doubles), bytes.fromhex("0600a0bf"), 16, 10) + doubles
    return _block("0220a0bf", bytes(extra), bytes(60) + sub)


def write_baf_d(folder, scans, calibration=CALIBRATION, declare=True):
    """scans: [(時刻 ms, 極性 +1/-1, 強度の整数の並び, TIC), ...]。folder(.d)を作って返す。"""
    os.makedirs(folder, exist_ok=True)
    parts = [_block("0340a0bf", b"", (COMPRESSOR_ID + b"\0" + COMPRESSOR_NAME + b"\0") if declare else b"other\0")]
    for time_ms, polarity, values, tic in scans:
        parts.append(scan_header(time_ms, polarity, calibration))
        parts.append(_block("0110a0bf", bytes(32), encode_profile(values)))
        line_extra = bytearray(40)
        struct.pack_into("<f", line_extra, 48 - 12, float(tic))
        parts.append(_block("0210a0bf", bytes(line_extra), bytes(16)))
    with open(os.path.join(folder, "analysis.baf"), "wb") as f:
        f.write(b"".join(parts))
    return str(folder)


def peaks_profile(n=4000, peaks=((1000, 40), (1004, 900), (2500, 70000)), noise_seed=0):
    rng = np.random.default_rng(noise_seed)
    y = np.zeros(n, dtype=np.int64)
    for center, height in peaks:
        for k in range(-6, 7):
            y[center + k] += int(height * np.exp(-0.5 * (k / 2.0) ** 2))
    y[rng.integers(0, n, 30)] += rng.integers(4, 20, 30)
    return y
