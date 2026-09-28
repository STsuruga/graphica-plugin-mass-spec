"""mzML の読み込み。GUI にも Graphica 本体にも依存しない。

スキャンの見出し(時刻・極性・TIC)は最初に全部読み、m/z と強度の配列は base64 の文字列のまま持って、
使うときに復号する。数百〜数千スキャンの profile を全部復号すると、GUI スレッドで待たせる時間とメモリが増えるため。
"""
import base64
import os
import zlib
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import numpy as np

_MZ_ARRAY = "MS:1000514"
_INTENSITY_ARRAY = "MS:1000515"
_FLOAT32 = "MS:1000521"
_FLOAT64 = "MS:1000523"
_ZLIB = "MS:1000574"
_NUMPRESS = {"MS:1002312", "MS:1002313", "MS:1002314", "MS:1002746", "MS:1002747", "MS:1002748"}
_MS_LEVEL = "MS:1000511"
_POSITIVE = "MS:1000130"
_NEGATIVE = "MS:1000129"
_CENTROID = "MS:1000127"
_PROFILE = "MS:1000128"
_TIC = "MS:1000285"
_SCAN_START_TIME = "MS:1000016"
_SECOND_UNITS = {"UO:0000010"}
_HOUR_UNITS = {"UO:0000032"}


class MzmlError(ValueError):
    pass


def _local(tag):
    return tag.rsplit("}", 1)[-1]


@dataclass
class EncodedArray:
    text: bytes
    dtype: str          # "<f4" / "<f8"
    compressed: bool
    length: int

    def decode(self):
        if not self.text:
            return np.zeros(0)
        raw = base64.b64decode(self.text)
        if self.compressed:
            raw = zlib.decompress(raw)
        values = np.frombuffer(raw, dtype=self.dtype)
        if self.length and len(values) != self.length:
            raise MzmlError(f"配列の長さが defaultArrayLength({self.length})と合いません({len(values)})")
        return values.astype(float)


@dataclass
class Scan:
    index: int
    native_id: str
    time_min: float
    ms_level: int
    polarity: int                 # +1 / -1 / 0(記載なし)
    centroid: bool | None         # None は記載なし
    tic: float
    mz_array: EncodedArray | None = field(default=None, repr=False)
    intensity_array: EncodedArray | None = field(default=None, repr=False)

    def arrays(self):
        """(m/z, 強度)。m/z の昇順。"""
        if self.mz_array is None or self.intensity_array is None:
            return np.zeros(0), np.zeros(0)
        mz = self.mz_array.decode()
        intensity = self.intensity_array.decode()
        if len(mz) != len(intensity):
            raise MzmlError(f"スキャン {self.native_id}: m/z と強度の点数が違います")
        if len(mz) > 1 and np.any(np.diff(mz) < 0):
            order = np.argsort(mz, kind="stable")
            mz, intensity = mz[order], intensity[order]
        return mz, intensity


@dataclass
class MzmlRun:
    path: str
    scans: list

    @property
    def name(self):
        base = os.path.basename(self.path.rstrip("\\/"))
        for ext in (".mzml", ".d"):
            if base.lower().endswith(ext):
                return base[:-len(ext)]
        return base

    def ms1_scans(self):
        return [s for s in self.scans if s.ms_level == 1]

    def times(self):
        return np.array([s.time_min for s in self.ms1_scans()])

    def tics(self):
        return np.array([s.tic for s in self.ms1_scans()])

    def polarity(self):
        """ファイル全体の極性。混在していれば 0。"""
        values = {s.polarity for s in self.ms1_scans()}
        return values.pop() if len(values) == 1 else 0

    def scans_in_range(self, t0, t1):
        lo, hi = min(t0, t1), max(t0, t1)
        return [s for s in self.ms1_scans() if lo <= s.time_min <= hi]

    def nearest_scan(self, t):
        scans = self.ms1_scans()
        if not scans:
            return None
        return min(scans, key=lambda s: abs(s.time_min - t))


def _params(element, groups):
    """要素直下の cvParam を {accession: (value, unitAccession)} にする。参照された paramGroup も展開する。"""
    params = {}
    for child in element:
        tag = _local(child.tag)
        if tag == "cvParam":
            params[child.get("accession")] = (child.get("value", ""), child.get("unitAccession", ""))
        elif tag == "referenceableParamGroupRef":
            params.update(groups.get(child.get("ref"), {}))
    return params


def _time_in_minutes(value, unit):
    t = float(value)
    if unit in _SECOND_UNITS:
        return t / 60.0
    if unit in _HOUR_UNITS:
        return t * 60.0
    return t


def _encoded_array(bda, groups, default_length):
    params = _params(bda, groups)
    if _FLOAT64 in params:
        dtype = "<f8"
    elif _FLOAT32 in params:
        dtype = "<f4"
    else:
        raise MzmlError("32/64 bit 浮動小数点以外の配列(整数・numpress など)には対応していません。"
                        "msconvert では numpress を使わずに変換してください")
    if _NUMPRESS & params.keys():
        raise MzmlError("numpress で圧縮された配列には対応していません。"
                        "msconvert では numpress を使わずに変換してください")
    compressed = _ZLIB in params
    binary = next((c for c in bda if _local(c.tag) == "binary"), None)
    text = (binary.text or "").strip().encode("ascii") if binary is not None else b""
    kind = "mz" if _MZ_ARRAY in params else "intensity" if _INTENSITY_ARRAY in params else None
    array_length = bda.get("arrayLength")
    length = int(array_length) if array_length else default_length
    return kind, EncodedArray(text, dtype, compressed, length)


def _parse_spectrum(element, groups):
    params = _params(element, groups)
    default_length = int(element.get("defaultArrayLength", "0") or 0)
    time_min = float("nan")
    for scan in element.iter():
        if _local(scan.tag) == "scan":
            scan_params = _params(scan, groups)
            if _SCAN_START_TIME in scan_params:
                time_min = _time_in_minutes(*scan_params[_SCAN_START_TIME])
            break
    mz_array = intensity_array = None
    for bda in element.iter():
        if _local(bda.tag) != "binaryDataArray":
            continue
        kind, array = _encoded_array(bda, groups, default_length)
        if kind == "mz":
            mz_array = array
        elif kind == "intensity":
            intensity_array = array
    polarity = 1 if _POSITIVE in params else -1 if _NEGATIVE in params else 0
    centroid = True if _CENTROID in params else False if _PROFILE in params else None
    ms_level = int(float(params[_MS_LEVEL][0])) if _MS_LEVEL in params else 1
    scan = Scan(int(element.get("index", "0") or 0), element.get("id", ""), time_min, ms_level,
                polarity, centroid, float("nan"), mz_array, intensity_array)
    if _TIC in params and params[_TIC][0]:
        scan.tic = float(params[_TIC][0])
    else:
        scan.tic = float(scan.arrays()[1].sum())
    return scan


def read_mzml(path, progress=None):
    """mzML を読む。progress(読んだスキャン数, 全スキャン数) を途中で呼ぶ(全数が分からなければ 0)。"""
    groups = {}
    scans = []
    total = 0
    try:
        context = ET.iterparse(path, events=("start", "end"))
        stack = []
        for event, element in context:
            tag = _local(element.tag)
            if event == "start":
                stack.append(tag)
                if tag == "spectrumList":
                    total = int(element.get("count", "0") or 0)
                continue
            stack.pop()
            if tag == "referenceableParamGroup":
                groups[element.get("id")] = _params(element, {})
            elif tag == "spectrum":
                scans.append(_parse_spectrum(element, groups))
                element.clear()
                if progress is not None and len(scans) % 50 == 0:
                    progress(len(scans), total)
            elif tag == "chromatogram":
                element.clear()
    except ET.ParseError as e:
        raise MzmlError(f"mzML として読めません: {e}") from e
    if not scans:
        raise MzmlError("スペクトルが1つもありません")
    if progress is not None:
        progress(len(scans), total)
    return MzmlRun(path, scans)
