"""質量分析 (MS) パック(Graphica プラグイン、P-316)。"""
from .analyzer import ANALYZER_NAME, PARAM_SCHEMA, analyze, load_mzml_file
from .plot_types import STICK, draw_stick

PANEL_NAME = "MS スペクトル"
HELP_MENU = "MS パック: 使い方"

HELP_TEXT = """\
MS パックは「MS スペクトル」パネルで使います(プラグインメニューからパネルを表示)。

開く: 「.d を開く…」で Bruker の .d フォルダを選ぶと、ProteoWizard の msconvert で mzML に変換して読みます
(2回目からは変換済みのものを使います)。mzML は「mzML…」か、本体のファイルを開く・ドラッグ&ドロップでも開けます。

TIC: 左ドラッグで試料の時間範囲、Shift+左ドラッグで背景の範囲。帯の端をドラッグで伸縮、帯の中で移動。
クリックでその時刻の1スキャン。

スペクトル: 左ドラッグで m/z 範囲に拡大、Ctrl+左ドラッグで矩形の拡大、Shift+左ドラッグで Δm/z を測る
(同位体の間隔なら電荷数も表示)。ピークをクリックするとラベルを固定。

共通: ホイールで拡大縮小(Shift+ホイールで縦)、中ボタンドラッグでパン、ダブルクリックで全体、Backspace で1つ前の表示、
右クリックでメニュー(プロットへの転送など)。

照合: パネルの「組成式から計算」か、実測スペクトルを選んで プラグイン ▸ 解析 ▸ 同位体パターンと照合。
"""


def _show_help(ctx):
    ctx.show_message(HELP_TEXT, "MS パック")


def _create_panel(ctx):
    from .panel import MassSpecPanel
    return MassSpecPanel(ctx)


def register(api):
    api.register_plot_type(STICK, draw_stick)
    api.register_importer([".mzML"], load_mzml_file, name="mzML(質量分析)")
    api.register_analyzer(ANALYZER_NAME, analyze, output_kind="table", param_schema=PARAM_SCHEMA)
    api.register_panel(PANEL_NAME, _create_panel, area="right")
    api.register_menu_action(HELP_MENU, _show_help)
