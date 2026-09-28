"""Bruker の .d(analysis.baf)を直接読む。GUI にも Graphica 本体にも依存しない。

Bruker のライブラリも ProteoWizard も使わず、データファイルの形式だけを読む。確かめたのは micrOTOF(TOF、
「DPCM/VLC based compressor」で圧縮した profile)の測定で、msconvert の mzML と強度は全点一致、m/z は 0.03 ppm 以内。
違う圧縮方式や、読めない形の測定は BafUnsupported を上げる(呼び出し側は msconvert に切り替える)。

ファイルは「大きさ(4 バイト)・種類(4 バイト)・見出しの長さ(4 バイト)」で始まるブロックの並び。
スキャンごとに、見出し(0220)・profile(0110)・centroid(0210)の3つが続く。
"""
import os
import re
import struct

import numpy as np

from .mzml import MzmlRun

KIND_SCAN_HEADER = bytes.fromhex("0220a0bf")
KIND_PROFILE = bytes.fromhex("0110a0bf")
KIND_LINE = bytes.fromhex("0210a0bf")
KIND_CALIBRATION = bytes.fromhex("0600a0bf")
COMPRESSOR_NAME = b"DPCM/VLC based compressor"
COMPRESSOR_ID = b"BCC4_4_4_8"
DECLARATION_BYTES = 1 << 16   # 圧縮方式の宣言はファイルの先頭にある

# 見出し(0220)の中の位置
HEADER_TIME_MS = 24       # uint32、取り込み開始からのミリ秒
HEADER_POLARITY = 32      # 0 = 正、1 = 負
LINE_TIC = 48             # centroid(0210)の中の TIC(float32)

# 差分の符号: 先頭の 0 が k 個(0〜6)、1、(5 + k) ビット。合わせた値 - 32 が z で、差は z/2(z が奇数なら負)。
# 0 が 7 個の後は 32 ビットの値そのもの。0 の差(100000)の連続は1つにまとめて読む。
_ZERO_RUN = "(?:100000)+"
_ESCAPE = "00000001[01]{32}"
_CODES = "|".join(f"{'0' * k}1[01]{{{5 + k}}}" for k in range(7))
_TOKEN = re.compile(f"{_ZERO_RUN}|{_ESCAPE}|{_CODES}")
_BYTE_BITS = [format(i, "08b") for i in range(256)]


class BafUnsupported(ValueError):
    """この形式は読めない(msconvert など別の方法で開く)。"""


def is_baf_folder(path):
    return os.path.isfile(os.path.join(path, "analysis.baf"))


def decode_profile(payload):
    """profile のブロックの中身(見出しの後)を強度の配列にする。"""
    if len(payload) < 13:
        raise BafUnsupported("profile のデータが短すぎます")
    n, first_count = struct.unpack_from("<II", payload, 4)
    if n == 0 or n > 50_000_000 or first_count != 1:
        raise BafUnsupported("profile の形式が想定と違います")
    out = np.empty(n, dtype=np.float64)
    value = payload[12]
    out[0] = value
    i = 1
    bits = "".join([_BYTE_BITS[b] for b in payload[13:]])
    for token in _TOKEN.findall(bits):
        if token[0] == "1":
            if len(token) > 6:           # 0 の差の連続
                k = len(token) // 6
                out[i:i + k] = value
                i += k
            else:
                z = int(token, 2) - 32
                value += -(z >> 1) if z & 1 else z >> 1
                out[i] = value
                i += 1
        elif len(token) == 40:           # 値そのもの
            value = int(token[8:], 2)
            out[i] = value
            i += 1
        else:
            z = int(token, 2) - 32
            value += -(z >> 1) if z & 1 else z >> 1
            out[i] = value
            i += 1
        if i >= n:
            break
    if i < n:
        raise BafUnsupported(f"profile を最後まで読めませんでした({i}/{n} 点)")
    return out[:n]


def _calibration(header):
    """見出しの中の較正のブロックから、(遅れ, 点の間隔, c1, c2, c3, ずれ)。"""
    at = header.find(KIND_CALIBRATION, 40)
    if at < 4:
        raise BafUnsupported("m/z の較正値が見つかりません")
    start = at - 4
    size, _kind, hsize, count = struct.unpack_from("<I4sII", header, start)
    if count < 7 or start + hsize + 7 * 8 > len(header):
        raise BafUnsupported("m/z の較正値の形が想定と違います")
    delay, interval, c2, c1, c3, _unused, offset = struct.unpack_from("<7d", header, start + hsize)
    if not (interval > 0 and c1 > 0):
        raise BafUnsupported("m/z の較正値が正しくありません")
    return delay, interval, c1, c2, c3, offset


def mz_axis(calibration, n):
    """TOF の較正式: 飛行時間 t = 遅れ + 間隔 × 点の番号、c3·x² + √(10¹²/c1)·x + (c2 − t) = 0 の解 x について m/z = x² − ずれ。"""
    delay, interval, c1, c2, c3, offset = calibration
    tof = delay + interval * np.arange(n, dtype=np.float64)
    b = np.sqrt(1e12 / c1)
    c = c2 - tof
    if c3 == 0:
        x = -c / b
    else:
        x = (-b + np.sqrt(b * b - 4.0 * c3 * c)) / (2.0 * c3)
    return x * x - offset


class BafScan:
    """mzml.Scan と同じ使い方ができるスキャン。配列は使うときに読んで解く。"""

    ms_level = 1
    centroid = False

    def __init__(self, path, index, time_min, polarity, tic, profile_at, calibration):
        self.path = path
        self.index = index
        self.native_id = f"scan={index + 1}"
        self.time_min = time_min
        self.polarity = polarity
        self.tic = tic
        self._profile_at = profile_at      # (位置, 大きさ, 見出しの長さ)
        self._calibration = calibration

    def arrays(self):
        pos, size, hsize = self._profile_at
        with open(self.path, "rb") as f:
            f.seek(pos)
            block = f.read(size)
        intensity = decode_profile(block[hsize:])
        return mz_axis(self._calibration, len(intensity)), intensity


def _blocks(f, file_size):
    """(位置, 大きさ, 種類, 見出しの長さ, 見出しのバイト列) を順に返す。profile の中身は読まない。"""
    pos = 0
    while pos + 12 <= file_size:
        f.seek(pos)
        head = f.read(12)
        size, kind, hsize = struct.unpack("<I4sI", head)
        if size < 12 or pos + size > file_size:
            break
        body = f.read(min(size, 4096) - 12) if kind != KIND_PROFILE else b""
        yield pos, size, kind, hsize, head + body
        pos += size


def read_baf(d_path, progress=None):
    """.d フォルダを読み、MzmlRun と同じ形で返す(配列は使うときに解く)。"""
    path = os.path.join(d_path, "analysis.baf")
    if not os.path.isfile(path):
        raise BafUnsupported("analysis.baf がありません(timsTOF の analysis.tdf などは未対応)")
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        declaration = f.read(DECLARATION_BYTES)
        if COMPRESSOR_NAME not in declaration or COMPRESSOR_ID not in declaration:
            raise BafUnsupported("この測定の圧縮方式には対応していません")
        scans = []
        current = None
        for pos, size, kind, hsize, raw in _blocks(f, file_size):
            if kind == KIND_SCAN_HEADER:
                current = {"header": raw}
                scans.append(current)
            elif current is not None and kind == KIND_PROFILE:
                current["profile"] = (pos, size, hsize)
            elif current is not None and kind == KIND_LINE and len(raw) >= LINE_TIC + 4:
                current["tic"] = struct.unpack_from("<f", raw, LINE_TIC)[0]
            if progress is not None and len(scans) % 50 == 0 and kind == KIND_SCAN_HEADER:
                progress(pos, file_size)
    result = []
    for i, s in enumerate(scans):
        if "profile" not in s:
            continue
        header = s["header"]
        time_ms = struct.unpack_from("<I", header, HEADER_TIME_MS)[0]
        polarity = {0: 1, 1: -1}.get(header[HEADER_POLARITY], 0)
        result.append(BafScan(path, i, time_ms / 60000.0, polarity, float(s.get("tic", float("nan"))),
                              (s["profile"]), _calibration(header)))
    if not result:
        raise BafUnsupported("profile のスペクトルがありません")
    if progress is not None:
        progress(file_size, file_size)
    return MzmlRun(d_path, result)
