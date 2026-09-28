"""ProteoWizard msconvert による .d → mzML の変換の準備とキャッシュ。実行そのものはパネルが QProcess で行う。

msconvert は利用者が各自インストールしたものを呼ぶだけで、このプラグインは同梱しない(ベンダーの DLL の再配布条件のため)。
"""
import glob
import hashlib
import os
import re
import shutil

# 変換の設定。変えたらキャッシュのキーが変わり、次に開くときに変換し直す。
# profile のまま(peakPicking なし)、64 bit の m/z(既定)、zlib 圧縮。numpress は mzml.py が読めないので使わない。
MSCONVERT_ARGS = ["--mzML", "--zlib"]
CACHE_KEEP = 30


def _version_key(path):
    numbers = re.findall(r"\d+", os.path.basename(os.path.dirname(path)))
    return [int(n) for n in numbers], os.path.getmtime(path)


def msconvert_candidates():
    """ProteoWizard の既定のインストール先にある msconvert.exe(新しい版から)。"""
    patterns = []
    local = os.environ.get("LOCALAPPDATA")
    if local:
        patterns.append(os.path.join(local, "Apps", "ProteoWizard*", "msconvert.exe"))
    for key in ("ProgramFiles", "ProgramW6432", "ProgramFiles(x86)"):
        base = os.environ.get(key)
        if base:
            patterns.append(os.path.join(base, "ProteoWizard", "ProteoWizard*", "msconvert.exe"))
            patterns.append(os.path.join(base, "ProteoWizard*", "msconvert.exe"))
    found = {os.path.normcase(os.path.abspath(p)): p for pattern in patterns for p in glob.glob(pattern)}
    return sorted(found.values(), key=_version_key, reverse=True)


def find_msconvert(configured=None):
    """設定された場所 → PATH → 既定のインストール先の順に探す。無ければ None。"""
    if configured and os.path.isfile(configured):
        return configured
    on_path = shutil.which("msconvert")
    if on_path:
        return on_path
    candidates = msconvert_candidates()
    return candidates[0] if candidates else None


def is_d_folder(path):
    return os.path.isdir(path) and path.rstrip("\\/").lower().endswith(".d")


def _raw_files(d_path):
    """変換結果に影響する中身(Bruker の測定ファイル)。"""
    names = ("analysis.baf", "analysis.tdf", "analysis.tdf_bin", "analysis.yep", "fid")
    return [os.path.join(d_path, n) for n in names if os.path.isfile(os.path.join(d_path, n))]


def cache_key(d_path, msconvert):
    """変換結果を左右する入力すべて: .d の場所と中身の大きさ・更新時刻、変換の設定、msconvert の場所と更新時刻。"""
    parts = [os.path.normcase(os.path.abspath(d_path)), " ".join(MSCONVERT_ARGS)]
    for f in _raw_files(d_path):
        st = os.stat(f)
        parts.append(f"{os.path.basename(f)}:{st.st_size}:{st.st_mtime_ns}")
    if msconvert and os.path.isfile(msconvert):
        parts.append(f"{os.path.normcase(os.path.abspath(msconvert))}:{os.stat(msconvert).st_mtime_ns}")
    return hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()[:16]


def _safe_stem(d_path):
    stem = os.path.basename(d_path.rstrip("\\/"))[:-2] or "measurement"
    return re.sub(r'[<>:"/\\|?*]', "_", stem)


def cached_mzml_path(d_path, cache_dir, msconvert):
    return os.path.join(cache_dir, f"{_safe_stem(d_path)}-{cache_key(d_path, msconvert)}.mzML")


def conversion_command(msconvert, d_path, work_dir):
    """(プログラム, 引数)。出力は作業フォルダに書き、終わってからキャッシュへ移す(途中で止まった変換を使わないため)。"""
    return msconvert, [d_path, *MSCONVERT_ARGS, "-o", work_dir, "--outfile", "converted.mzML"]


def work_dir_for(final_path):
    return final_path + ".part"


def finish_conversion(final_path):
    """作業フォルダの出力をキャッシュへ移す。出力が無ければ None。"""
    work = work_dir_for(final_path)
    produced = os.path.join(work, "converted.mzML")
    if not os.path.isfile(produced):
        shutil.rmtree(work, ignore_errors=True)
        return None
    os.replace(produced, final_path)
    shutil.rmtree(work, ignore_errors=True)
    return final_path


def prune_cache(cache_dir, keep=CACHE_KEEP):
    """古いキャッシュを消す。新しい順に keep 個を残す。"""
    files = sorted(glob.glob(os.path.join(cache_dir, "*.mzML")), key=os.path.getmtime, reverse=True)
    for path in files[keep:]:
        try:
            os.remove(path)
        except OSError:
            pass
