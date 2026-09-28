import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mass_spec import convert  # noqa: E402


def _d_folder(tmp_path, name="sample.d", content=b"abc"):
    d = tmp_path / name
    d.mkdir()
    (d / "analysis.baf").write_bytes(content)
    return str(d)


def _fake_install(tmp_path, monkeypatch, versions):
    local = tmp_path / "local"
    for v in versions:
        folder = local / "Apps" / f"ProteoWizard {v} 64-bit"
        folder.mkdir(parents=True)
        (folder / "msconvert.exe").write_bytes(b"")
    monkeypatch.setenv("LOCALAPPDATA", str(local))
    for key in ("ProgramFiles", "ProgramW6432", "ProgramFiles(x86)"):
        monkeypatch.setenv(key, str(tmp_path / "nothing"))
    monkeypatch.setattr(convert.shutil, "which", lambda name: None)


def test_find_msconvert_prefers_the_configured_path_then_the_newest_install(tmp_path, monkeypatch):
    _fake_install(tmp_path, monkeypatch, ["3.0.25100.aaa", "3.0.26267.699dd48"])
    assert "3.0.26267" in convert.find_msconvert()
    configured = tmp_path / "my" / "msconvert.exe"
    configured.parent.mkdir()
    configured.write_bytes(b"")
    assert convert.find_msconvert(str(configured)) == str(configured)
    assert "3.0.26267" in convert.find_msconvert(str(tmp_path / "missing.exe"))


def test_find_msconvert_returns_none_without_an_install(tmp_path, monkeypatch):
    _fake_install(tmp_path, monkeypatch, [])
    assert convert.find_msconvert() is None


def test_is_d_folder(tmp_path):
    assert convert.is_d_folder(_d_folder(tmp_path))
    assert not convert.is_d_folder(str(tmp_path))
    assert not convert.is_d_folder(str(tmp_path / "x.d"))


def test_cache_key_changes_with_every_input(tmp_path):
    d = _d_folder(tmp_path)
    exe = tmp_path / "msconvert.exe"
    exe.write_bytes(b"")
    key = convert.cache_key(d, str(exe))
    assert convert.cache_key(d, str(exe)) == key
    time.sleep(0.01)
    with open(os.path.join(d, "analysis.baf"), "ab") as f:
        f.write(b"more")
    changed = convert.cache_key(d, str(exe))
    assert changed != key
    other = _d_folder(tmp_path, "other.d", b"abcmore")
    assert convert.cache_key(other, str(exe)) != changed
    exe2 = tmp_path / "v2" / "msconvert.exe"
    exe2.parent.mkdir()
    exe2.write_bytes(b"")
    assert convert.cache_key(d, str(exe2)) != changed


def test_cached_path_command_and_finish(tmp_path):
    d = _d_folder(tmp_path, "odd name.d")
    cache = tmp_path / "cache"
    cache.mkdir()
    final = convert.cached_mzml_path(d, str(cache), None)
    assert os.path.dirname(final) == str(cache) and final.endswith(".mzML")
    program, args = convert.conversion_command("msconvert.exe", d, convert.work_dir_for(final))
    assert program == "msconvert.exe"
    assert args[0] == d and "--zlib" in args and args[args.index("-o") + 1] == convert.work_dir_for(final)
    assert convert.finish_conversion(final) is None
    work = convert.work_dir_for(final)
    os.makedirs(work)
    with open(os.path.join(work, "converted.mzML"), "w") as f:
        f.write("x")
    assert convert.finish_conversion(final) == final
    assert os.path.isfile(final) and not os.path.exists(work)


def test_prune_cache_keeps_the_newest(tmp_path):
    for i in range(5):
        p = tmp_path / f"{i}.mzML"
        p.write_text("x")
        os.utime(p, (1000 + i, 1000 + i))
    convert.prune_cache(str(tmp_path), keep=2)
    assert sorted(os.listdir(tmp_path)) == ["3.mzML", "4.mzML"]
