"""質量分析 (MS) パック(Graphica プラグイン、P-316)。"""
from .analyzer import ANALYZER_NAME, PARAM_SCHEMA, analyze, load_mzml_file
from .plot_types import STICK, draw_stick


def register(api):
    api.register_plot_type(STICK, draw_stick)
    api.register_importer([".mzML"], load_mzml_file, name="mzML(質量分析)")
    api.register_analyzer(ANALYZER_NAME, analyze, output_kind="table", param_schema=PARAM_SCHEMA)
