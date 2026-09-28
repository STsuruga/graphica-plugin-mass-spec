"""組成式・付加イオン・同位体パターンの計算。GUI にも Graphica 本体にも依存しない。"""
import re
from dataclasses import dataclass

import numpy as np

from .isotopes import ISOTOPE_ROWS

ELECTRON_MASS = 5.48579909065e-4  # u(CODATA 2018)

_SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
    "Mc Lv Ts Og"
).split()


def _nist_value(text):
    return float(re.match(r"[\d.]+", text).group(0))


def _build_isotope_table():
    """{元素記号: [(質量数, 質量, 存在比), ...]}。天然に存在する同位体だけ。"""
    table = {}
    for z, a, mass, abundance in ISOTOPE_ROWS:
        if abundance is None:
            continue
        table.setdefault(_SYMBOLS[z - 1], []).append((a, _nist_value(mass), _nist_value(abundance)))
    # 重水素は D と書かれることが多いので、2H だけからなる元素として扱う
    h2 = next(row for row in table["H"] if row[0] == 2)
    table["D"] = [(2, h2[1], 1.0)]
    for symbol, rows in table.items():
        total = sum(r[2] for r in rows)
        table[symbol] = [(a, m, p / total) for a, m, p in sorted(rows)]
    return table


ISOTOPES = _build_isotope_table()


class FormulaError(ValueError):
    pass


_TOKEN = re.compile(r"\s*(?:([A-Z][a-z]?)|([(\[{])|([)\]}])|(\d+))")


def parse_formula(text):
    """組成式を {元素記号: 個数} にする。括弧の入れ子と D(重水素)に対応する。"""
    text = text.strip()
    if not text:
        raise FormulaError("組成式が空です")
    stack = [{}]
    pos = 0
    last = None  # 直前の元素記号か、閉じた括弧の中身(個数を掛ける対象)
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise FormulaError(f"組成式を読めません: {text!r}({pos + 1} 文字目)")
        symbol, opening, closing, number = m.groups()
        pos = m.end()
        if symbol:
            if symbol not in ISOTOPES:
                if symbol in _SYMBOLS:
                    raise FormulaError(f"{symbol} は天然の同位体組成が無いので計算できません")
                raise FormulaError(f"元素記号ではありません: {symbol}")
            stack[-1][symbol] = stack[-1].get(symbol, 0) + 1
            last = ("element", symbol)
        elif opening:
            stack.append({})
            last = None
        elif closing:
            if len(stack) == 1:
                raise FormulaError(f"括弧の対応がとれていません: {text!r}")
            group = stack.pop()
            for k, v in group.items():
                stack[-1][k] = stack[-1].get(k, 0) + v
            last = ("group", group)
        else:
            n = int(number)
            if last is None:
                raise FormulaError(f"数字の前に元素がありません: {text!r}")
            kind, value = last
            if kind == "element":
                stack[-1][value] += n - 1
            else:
                for k, v in value.items():
                    stack[-1][k] += v * (n - 1)
            last = None
    if len(stack) != 1:
        raise FormulaError(f"括弧が閉じていません: {text!r}")
    counts = {k: v for k, v in stack[0].items() if v}
    if not counts:
        raise FormulaError("組成式に元素がありません")
    return counts


def _hill_order(counts):
    keys = list(counts)
    if "C" in counts:
        head = ["C"] + (["H"] if "H" in counts else [])
        return head + sorted(k for k in keys if k not in head)
    return sorted(keys)


def format_formula(counts):
    """Hill 方式の並び(C、H、残りはアルファベット順)で書く。"""
    return "".join(k + (str(counts[k]) if counts[k] != 1 else "") for k in _hill_order(counts))


def monoisotopic_mass(counts):
    """各元素で最も存在比の大きい同位体の質量の和。"""
    return sum(n * max(ISOTOPES[el], key=lambda r: r[2])[1] for el, n in counts.items())


def average_mass(counts):
    """同位体の存在比で重み付けした平均分子量。"""
    return sum(n * sum(m * p for _a, m, p in ISOTOPES[el]) for el, n in counts.items())


# ---------------------------------------------------------------- 付加イオン

@dataclass(frozen=True)
class Adduct:
    notation: str
    multimer: int          # nM の n
    delta: dict            # 足し引きする原子 {元素記号: 個数}(負は引く)
    charge: int            # 符号付きの電荷数


_ADDUCT = re.compile(r"^\[?\s*(\d*)\s*M\s*((?:[+-]\s*\d*\s*[A-Za-z0-9()]+\s*)*)\]?\s*(\d*)\s*([+-])\s*([•.·*]?)$")
_ADDUCT_TERM = re.compile(r"([+-])\s*(\d*)\s*([A-Za-z0-9()]+)")

POSITIVE_STANDARD = ("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+•")
NEGATIVE_STANDARD = ("[M-H]-", "[M+Cl]-", "[M]-•", "[M+HCOO]-")


def _normalize_signs(text):
    return (text.replace("−", "-").replace("–", "-").replace("＋", "+")
            .replace("⁺", "+").replace("⁻", "-"))


def parse_adduct(text):
    """[M+H]+、[2M+Na]+、[M+2H]2+、[M-H-H2O]-、[M]+• のような表記を読む。"""
    s = _normalize_signs(text.strip())
    m = _ADDUCT.match(s)
    if not m:
        raise FormulaError(f"付加イオンの表記を読めません: {text!r}([M+H]+ の形で書きます)")
    multimer, terms, z_digits, z_sign, _radical = m.groups()
    delta = {}
    for sign, count, formula in _ADDUCT_TERM.findall(terms):
        factor = (int(count) if count else 1) * (1 if sign == "+" else -1)
        for el, n in parse_formula(formula).items():
            delta[el] = delta.get(el, 0) + factor * n
    charge = (int(z_digits) if z_digits else 1) * (1 if z_sign == "+" else -1)
    if charge == 0:
        raise FormulaError(f"電荷が 0 です: {text!r}")
    return Adduct(text.strip(), int(multimer) if multimer else 1,
                  {k: v for k, v in delta.items() if v}, charge)


def split_adducts(text):
    """カンマ・セミコロン・空白区切りの表記の並びを分ける。"""
    parts = re.split(r"[,;、\s]+(?=\[|\d*M)", _normalize_signs(text.strip()))
    return [p.strip().rstrip(",;、") for p in parts if p.strip()]


def ion_composition(counts, adduct):
    """付加イオンの組成。引く原子が足りなければ FormulaError。"""
    ion = {el: n * adduct.multimer for el, n in counts.items()}
    for el, n in adduct.delta.items():
        ion[el] = ion.get(el, 0) + n
    short = [el for el, n in ion.items() if n < 0]
    if short:
        raise FormulaError(f"{adduct.notation}: {', '.join(short)} が足りません")
    ion = {el: n for el, n in ion.items() if n}
    if not ion:
        raise FormulaError(f"{adduct.notation}: 原子が残りません")
    return ion


def ion_mz(neutral_mass, adduct, counts_mass_delta):
    """中性の質量と付加分の質量から m/z。電子の質量を電荷の分だけ増減する。"""
    return (neutral_mass * adduct.multimer + counts_mass_delta - adduct.charge * ELECTRON_MASS) / abs(adduct.charge)


# ---------------------------------------------------------------- 同位体パターン

def _element_distribution(symbol):
    """(最も軽い同位体からの質量数の差ごとの確率, 確率×質量)。"""
    rows = ISOTOPES[symbol]
    a0 = rows[0][0]
    size = rows[-1][0] - a0 + 1
    p = np.zeros(size)
    pm = np.zeros(size)
    for a, mass, prob in rows:
        p[a - a0] += prob
        pm[a - a0] += prob * mass
    return p, pm


def _combine(d1, d2, prune):
    p1, pm1 = d1
    p2, pm2 = d2
    p = np.convolve(p1, p2)
    pm = np.convolve(pm1, p2) + np.convolve(p1, pm2)
    keep = np.nonzero(p > p.max() * prune)[0]
    last = keep[-1] + 1 if len(keep) else 1
    return p[:last], pm[:last]


def _power(dist, n, prune):
    result = (np.array([1.0]), np.array([0.0]))
    base = dist
    while n:
        if n & 1:
            result = _combine(result, base, prune)
        n >>= 1
        if n:
            base = _combine(base, base, prune)
    return result


def isotope_pattern(counts, charge=1, min_relative=1e-3, prune=1e-12):
    """同位体パターンを、質量数の差(M, M+1, ...)ごとに束ねたピークの並びで返す。

    分解能が数万以下の装置では同じ質量数の差の微細構造は分かれないので、
    確率で重み付けした平均の質量を1本のピークにする。
    Returns: (m/z の配列, 最大を 100 とした相対強度の配列, 最も軽いピークからの質量数の差の配列)。
    電子の質量は含まない(付加イオンの m/z は ion_pattern で求める)。
    """
    total = (np.array([1.0]), np.array([0.0]))
    for el, n in counts.items():
        total = _combine(total, _power(_element_distribution(el), n, prune), prune)
    p, pm = total
    nonzero = p > 0
    mass = np.where(nonzero, pm / np.where(nonzero, p, 1.0), 0.0)
    rel = p / p.max() * 100.0
    keep = rel >= min_relative * 100.0
    offsets = np.nonzero(keep)[0]
    return mass[keep] / abs(charge), rel[keep], offsets


@dataclass
class IonPattern:
    adduct: Adduct
    composition: dict
    mz: np.ndarray            # 束ねた各ピークの m/z(電子の質量を補正済み)
    relative: np.ndarray      # 最大 100
    offsets: np.ndarray       # 最も軽いピークからの質量数の差
    monoisotopic_mz: float    # 各元素の主同位体だけからなるイオンの m/z


def ion_pattern(counts, adduct, min_relative=1e-3):
    ion = ion_composition(counts, adduct)
    z = adduct.charge
    mz, rel, offsets = isotope_pattern(ion, z, min_relative)
    correction = -z * ELECTRON_MASS / abs(z)
    return IonPattern(adduct, ion, mz + correction, rel, offsets,
                      monoisotopic_mass(ion) / abs(z) + correction)


def gaussian_profile(mz, relative, resolution, points_per_fwhm=20, span_fwhm=4.0):
    """各ピークに FWHM = m/z ÷ 分解能 のガウスを置いたプロファイル。"""
    mz = np.asarray(mz, dtype=float)
    relative = np.asarray(relative, dtype=float)
    fwhm = mz / resolution
    step = fwhm.min() / points_per_fwhm
    lo = (mz - span_fwhm * fwhm).min()
    hi = (mz + span_fwhm * fwhm).max()
    x = np.arange(lo, hi + step, step)
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    y = np.zeros_like(x)
    for center, height, s in zip(mz, relative, sigma):
        near = np.abs(x - center) <= 6 * s
        y[near] += height * np.exp(-0.5 * ((x[near] - center) / s) ** 2)
    return x, y
