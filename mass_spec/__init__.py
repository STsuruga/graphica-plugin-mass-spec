"""質量分析 (MS) パック(Graphica プラグイン、P-316)。"""
from .analyzer import ANALYZER_NAME, PARAM_SCHEMA, analyze, load_baf_file, load_mzml_file
from .plot_types import STICK, draw_stick

PANEL_NAME = "MS スペクトル"
OPEN_MENU = "MS ビューアを開く"
HELP_MENU = "MS パック: 使い方"


def _open_viewer(ctx):
    """このタブの MS ビューアのウィンドウを開く。ビューアはタブごとのパネルが持っている。"""
    from .panel import viewer_for
    viewer = viewer_for(ctx)
    if viewer is None:
        ctx.show_error("このタブの MS ビューアが見つかりません。プラグイン ▸ パネル ▸ MS スペクトル から開いてください。",
                       "MS パック")
        return
    viewer.open_window()


def _show_help(ctx):
    from .panel import HELP_TEXT
    ctx.show_message(f"プラグイン ▸ {OPEN_MENU} で MS ビューアのウィンドウを開きます"
                     "(クイックアクセスツールバーにピン留めするとボタンから開けます)。\n\n" + HELP_TEXT, "MS パック")


def _create_panel(ctx):
    from .panel import MassSpecPanel
    return MassSpecPanel(ctx)


def register(api):
    api.register_plot_type(STICK, draw_stick)
    api.register_importer([".mzML"], load_mzml_file, name="mzML(質量分析)")
    api.register_importer([".baf"], load_baf_file, name="Bruker BAF(質量分析)")
    api.register_analyzer(ANALYZER_NAME, analyze, output_kind="table", param_schema=PARAM_SCHEMA)
    api.register_panel(PANEL_NAME, _create_panel, area="right")
    api.register_menu_action(OPEN_MENU, _open_viewer)
    api.register_menu_action(HELP_MENU, _show_help)
