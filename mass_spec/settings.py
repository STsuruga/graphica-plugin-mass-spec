"""パネルの設定。プラグインの保存フォルダ(ctx.data_dir)の settings.json に置く。"""
import json
import logging
import os

logger = logging.getLogger(__name__)

FILE_NAME = "settings.json"

DEFAULTS = {
    "msconvert_path": "",
    "label_mode": "top",            # "top" / "percent" / "none"
    "label_top_n": 10,
    "label_percent": 5.0,
    "mz_decimals": 4,
    "intensity": "none",            # "none" / "relative" / "absolute"
    "ppm_decimals": 1,
    "subtract_background": True,
    "clip_negative": True,
    "transfer_view_only": False,
    "centroid_min_percent": 1.0,
    "subplot_target": 0,
    "formula": "",
    "preset": 0,                    # 付加イオンの標準セットの番号
    "extra_adducts": "",
    "resolution": 0.0,
    "tolerance_ppm": 50.0,
    "detect_percent": 0.5,
}


def load_settings(data_dir):
    settings = dict(DEFAULTS)
    path = os.path.join(data_dir, FILE_NAME)
    try:
        with open(path, encoding="utf-8") as f:
            stored = json.load(f)
    except FileNotFoundError:
        return settings
    except (OSError, ValueError):
        logger.warning("MS パックの設定を読めませんでした: %s", path)
        return settings
    for key, default in DEFAULTS.items():
        value = stored.get(key, default)
        # 壊れた値は既定に戻す(型が違えば使わない)
        if isinstance(default, bool):
            ok = isinstance(value, bool)
        elif isinstance(default, (int, float)):
            ok = isinstance(value, (int, float)) and not isinstance(value, bool)
        else:
            ok = isinstance(value, type(default))
        settings[key] = value if ok else default
    return settings


def save_settings(data_dir, settings):
    path = os.path.join(data_dir, FILE_NAME)
    try:
        os.makedirs(data_dir, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump({k: settings[k] for k in DEFAULTS}, f, ensure_ascii=False, indent=2)
    except OSError:
        logger.warning("MS パックの設定を保存できませんでした: %s", path)
