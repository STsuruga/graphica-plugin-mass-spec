"""MS ビューア(ドックパネル「MS スペクトル」)。

上端にメニューバー、左に開いた測定の一覧(チェックで TIC の表示、クリックで2段目の対象)、右に3段:
TIC(1段目)、選んだ時間範囲のスペクトル(2段目)、実測のコピーや計算パターンを1つずつ置く枠(3段目)。
照合・時間範囲・表示設定はメニューから別ウィンドウで開く(windows.py)。
"""
import os
import time

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import to_hex
from PySide6.QtCore import QProcess, Qt, QTimer
from PySide6.QtGui import QColor, QGuiApplication, QIcon, QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QFileDialog, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QMenu, QMenuBar, QMessageBox,
    QPushButton, QScrollArea, QSplitter, QSpinBox, QVBoxLayout, QWidget,
)

from . import convert
from .chemistry import (
    NEGATIVE_STANDARD, POSITIVE_STANDARD, FormulaError, average_mass, format_formula, gaussian_profile,
    ion_pattern, monoisotopic_mass, parse_adduct, parse_formula, split_adducts,
)
from .datasets import (
    centroid_dataset, pattern_datasets, range_text, spectrum_dataset, stick_dataset, tic_dataset,
)
from .matching import DEFAULT_RESOLUTION, match_all
from .mzml import MzmlError, read_mzml
from .plots import PanePlot, SpectrumPlot, TicPlot
from .settings import load_settings, save_settings
from .spectra import LabelFormat, average_spectrum, compact_zeros, find_peaks, subtract_background
from .windows import CalcWindow, RangeWindow, SettingsWindow

PANEL_NAME = "MS スペクトル"
MAX_PANES = 8
# ノイズの極大まで拾うとピークが数万になるので、最大の 0.1% 未満は拾わない
PEAK_FLOOR = 1e-3
RECOMPUTE_DELAY_MS = 250

# 測定ごとの TIC の色。Qt に渡すので #rrggbb にしておく
RUN_COLORS = [to_hex(c) for c in colormaps["tab10"].colors]


def _header(text):
    """各段の題名の帯(DataAnalysis の窓の題名に似せる)。色は本体のテーマのパレットから取る。"""
    label = QLabel(text)
    label.setStyleSheet("QLabel { background: palette(highlight); color: palette(highlighted-text);"
                        " font-weight: 600; padding: 2px 6px; }")
    return label


def _color_icon(color):
    pixmap = QPixmap(12, 12)
    pixmap.fill(QColor(color))
    return QIcon(pixmap)


MSCONVERT_MISSING = (
    "ProteoWizard の msconvert が見つかりません。\n\n"
    "https://proteowizard.sourceforge.io/download.html から Windows 64-bit 版"
    "(ベンダー形式を読める版)をインストールしてください。\n"
    "別の場所に入れた場合は、次の画面で msconvert.exe を選べます。")

HELP_TEXT = """\
MS ビューア(MS スペクトル パネル)

左の一覧: チェックでその測定の TIC を表示、名前をクリックでその測定を2段目の対象にする。
左の一覧: チェックでその測定の TIC を表示、名前をクリックでその測定を2段目の対象にする。
1段目 TIC: 左ドラッグで試料の時間範囲、Shift+左ドラッグで背景の範囲。帯の端をドラッグで伸縮、帯の中で移動。
クリックでその時刻の1スキャン。
2段目 スペクトル: 選んだ範囲の平均(背景を引く設定なら差し引き後)。左ドラッグでその m/z 範囲に拡大、
Ctrl+左ドラッグで矩形の拡大、Shift+左ドラッグで Δm/z を測る(同位体の間隔なら電荷数も)。ピークをクリックでラベルを固定。
右クリックで「3段目の枠にコピー」や転送。
3段目 枠: 実測のコピーか計算パターンを1つずつ。枠の数と、横軸を枠どうしで同期するかは3段目の下で切り替える。

軸の上: ホイールで拡大縮小、左ドラッグで表示範囲をずらす。
グラフの中: ホイールで拡大縮小(Shift で縦)、中ボタンドラッグでパン、ダブルクリックで全体、Backspace で1つ前。

解析 ▸ 同位体パターンの照合 で組成式から同位体パターンと付加イオンを計算し、実測と照合する(結果は枠や本体のプロットへ)。
表示 ▸ 表示設定 でラベルの本数・桁数、背景、転送の設定。
"""


class MassSpecPanel(QWidget):
    def __init__(self, ctx, parent=None):
        super().__init__(parent)
        self.ctx = ctx
        self.settings = load_settings(ctx.data_dir)
        self.runs = []
        self.ranges = {}            # run.path -> {"sample": (t0, t1), "background": (t0, t1) | None}
        self.spectrum = None        # {"mz", "y", "peaks", "name", "provenance", "run"}
        self.calc = None            # {"counts", "formula", "results", "adducts"}
        self.panes = []
        self._process = None
        self._pending = None        # 変換中の (d のパス, キャッシュの mzML, 表示名)
        self._convert_started = 0.0
        self._updating = False
        self._syncing = False
        self._recompute_timer = QTimer(self)
        self._recompute_timer.setSingleShot(True)
        self._recompute_timer.timeout.connect(self.compute_spectrum)
        self._elapsed_timer = QTimer(self)
        self._elapsed_timer.timeout.connect(self._show_conversion_progress)
        self._build()
        self.calc_window = CalcWindow(self)
        self.settings_window = SettingsWindow(self)
        self.range_window = RangeWindow(self)
        self._apply_label_settings()
        self.set_pane_count(self.settings["pane_count"])

    # ================================================================ 画面
    def _build(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 0, 4, 4)
        outer.setSpacing(4)
        self.menu_bar = QMenuBar(self)
        outer.setMenuBar(self.menu_bar)
        self._build_menus()

        self.pane_spin = QSpinBox()
        self.pane_spin.setRange(1, MAX_PANES)
        self.pane_spin.setValue(self.settings["pane_count"])
        self.pane_spin.valueChanged.connect(self.set_pane_count)
        self.sync_check = QCheckBox("横軸を同期")
        self.sync_check.setToolTip("3段目の枠どうしで m/z の範囲をそろえる")
        self.sync_check.setChecked(self.settings["sync_panes"])
        self.sync_check.toggled.connect(self._on_sync_toggled)
        self.cancel_button = QPushButton("変換を中止")
        self.cancel_button.clicked.connect(self.cancel_conversion)
        self.cancel_button.hide()
        self.last_status = ""

        # 左に測定の一覧、右に3段のグラフ
        self.main_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.main_splitter.setChildrenCollapsible(False)
        outer.addWidget(self.main_splitter, 1)
        self.run_list = QListWidget()
        self.run_list.setMinimumWidth(140)
        self.run_list.setToolTip("チェックで TIC を表示、名前をクリックでその測定のスペクトルを2段目に出す")
        self.run_list.currentRowChanged.connect(self._on_run_selected)
        self.run_list.itemChanged.connect(self._on_run_item_changed)
        self.main_splitter.addWidget(self.run_list)

        self.splitter = QSplitter(Qt.Orientation.Vertical)
        self.splitter.setChildrenCollapsible(False)
        self.main_splitter.addWidget(self.splitter)
        self.main_splitter.setStretchFactor(1, 1)
        self.main_splitter.setSizes([170, 800])

        tic_box = QWidget()
        tic_layout = QVBoxLayout(tic_box)
        tic_layout.setContentsMargins(0, 0, 0, 0)
        tic_layout.setSpacing(0)
        self.tic_header = _header("クロマトグラム")
        tic_layout.addWidget(self.tic_header)
        self.tic_plot = TicPlot()
        self.tic_plot.range_changing.connect(self._on_tic_range_changing)
        self.tic_plot.range_changed.connect(self._on_tic_range_changed)
        self.tic_plot.scan_clicked.connect(self._on_scan_clicked)
        self.tic_plot.context_requested.connect(self._tic_menu)
        tic_layout.addWidget(self.tic_plot, 1)
        self.splitter.addWidget(tic_box)

        spectrum_box = QWidget()
        spectrum_layout = QVBoxLayout(spectrum_box)
        spectrum_layout.setContentsMargins(0, 0, 0, 0)
        spectrum_layout.setSpacing(0)
        self.spectrum_header = _header("スペクトル")
        spectrum_layout.addWidget(self.spectrum_header)
        self.spectrum_plot = SpectrumPlot()
        self.spectrum_plot.hover_text.connect(self._show_readout)
        self.spectrum_plot.measured.connect(self._show_readout)
        self.spectrum_plot.context_requested.connect(self._spectrum_menu)
        spectrum_layout.addWidget(self.spectrum_plot, 1)
        self.readout_label = QLabel(" ")
        self.readout_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        spectrum_layout.addWidget(self.readout_label)
        self.splitter.addWidget(spectrum_box)

        pane_box = QWidget()
        pane_box_layout = QVBoxLayout(pane_box)
        pane_box_layout.setContentsMargins(0, 0, 0, 0)
        pane_box_layout.setSpacing(0)
        self.pane_header = _header("比較スペクトル")
        pane_box_layout.addWidget(self.pane_header)
        self.pane_area = QScrollArea()
        self.pane_area.setWidgetResizable(True)
        self.pane_container = QWidget()
        self.pane_layout = QVBoxLayout(self.pane_container)
        self.pane_layout.setContentsMargins(0, 0, 0, 0)
        self.pane_area.setWidget(self.pane_container)
        pane_box_layout.addWidget(self.pane_area, 1)
        row = QHBoxLayout()
        row.addWidget(QLabel("枠数"))
        row.addWidget(self.pane_spin)
        row.addWidget(self.sync_check)
        row.addStretch(1)
        pane_box_layout.addLayout(row)
        self.splitter.addWidget(pane_box)

        # 変換中だけ右下に出す(数秒かかるので、何も出ないと止まったように見える)
        self.progress_row = QWidget()
        progress_layout = QHBoxLayout(self.progress_row)
        progress_layout.setContentsMargins(0, 0, 0, 0)
        progress_layout.addStretch(1)
        self.progress_label = QLabel("")
        progress_layout.addWidget(self.progress_label)
        progress_layout.addWidget(self.cancel_button)
        self.progress_row.hide()
        outer.addWidget(self.progress_row)
        self.splitter.setSizes([180, 260, 360])

    def _build_menus(self):
        m = self.menu_bar.addMenu("ファイル")
        self.open_d_action = m.addAction(".d を開く…", self._choose_d_folder)
        m.addAction("mzML を開く…", self._choose_mzml)
        m.addSeparator()
        m.addAction("選んでいる測定を閉じる", self.close_current_run)
        self.cancel_action = m.addAction("変換を中止", self.cancel_conversion)
        self.cancel_action.setEnabled(False)
        m = self.menu_bar.addMenu("表示")
        m.addAction("すべてのグラフを全体表示", self.reset_all_views)
        m.addAction("3段目の枠をすべて空にする", self.clear_all_panes)
        m.addSeparator()
        m.addAction("時間範囲の詳細設定…", lambda: self._open_window(self.range_window))
        m.addAction("表示設定…", self.open_settings_window)
        m = self.menu_bar.addMenu("解析")
        m.addAction("同位体パターンの照合…", self.open_calc_window)
        m = self.menu_bar.addMenu("転送")
        m.addAction("TIC", self.transfer_tic)
        m.addAction("スペクトル(2段目)", self.transfer_spectrum)
        m.addAction("centroid(2段目)", self.transfer_centroid)
        m.addAction("ピークのラベルだけ(2段目)", lambda: self.transfer_centroid(labels_only=True))
        m.addAction("計算パターン", self.transfer_calculation)
        m = self.menu_bar.addMenu("ヘルプ")
        m.addAction("使い方", lambda: self.ctx.show_message(HELP_TEXT, "MS パック"))
        # PySide6 はメニューを Python 側で持っていないと消すことがあるので、項目を含めて保持する
        self._menu_keepalive = [menu for menu in self.menu_bar.findChildren(QMenu)]

    @staticmethod
    def _open_window(window):
        window.show()
        window.raise_()
        window.activateWindow()

    def open_calc_window(self):
        self._open_window(self.calc_window)

    def open_settings_window(self):
        self._open_window(self.settings_window)

    # ================================================================ 設定
    def _save(self):
        save_settings(self.ctx.data_dir, self.settings)

    def apply_settings(self, **values):
        """設定ウィンドウなどからの変更を反映して保存する。"""
        background_changed = any(k in values and values[k] != self.settings[k]
                                 for k in ("subtract_background", "clip_negative"))
        self.settings.update(values)
        self._save()
        self._apply_label_settings()
        for plot in [self.spectrum_plot, *self.panes]:
            plot.update_labels()
        if self.calc is not None:
            self._show_calc_table()
        if background_changed:
            self.compute_spectrum()

    def label_format(self):
        return LabelFormat(mz_decimals=self.settings["mz_decimals"], intensity=self.settings["intensity"],
                           ppm_decimals=self.settings["ppm_decimals"])

    def _apply_label_settings(self):
        for plot in [self.spectrum_plot, *self.panes]:
            plot.label_mode = self.settings["label_mode"]
            plot.label_top_n = self.settings["label_top_n"]
            plot.label_percent = self.settings["label_percent"]
            plot.label_format = self.label_format()

    def _show_readout(self, text):
        self.readout_label.setText(text)

    def _status(self, text):
        """画面には出さない(状態の行は置かない)。直近の内容はテストと調査のために持っておく。"""
        self.last_status = text

    # ================================================================ 3段目の枠
    def pane_count(self):
        return len(self.panes)

    def set_pane_count(self, count):
        count = max(1, min(MAX_PANES, int(count)))
        while len(self.panes) < count:
            pane = PanePlot(len(self.panes) + 1)
            pane.hover_text.connect(self._show_readout)
            pane.measured.connect(self._show_readout)
            pane.context_requested.connect(lambda pos, p=pane: self._pane_menu(p, pos))
            pane.x_range_changed.connect(lambda lo, hi, p=pane: self._on_pane_x_changed(p, lo, hi))
            self.pane_layout.addWidget(pane)
            self.panes.append(pane)
        while len(self.panes) > count:
            pane = self.panes.pop()
            self.pane_layout.removeWidget(pane)
            pane.deleteLater()
        self._apply_label_settings()
        if self.pane_spin.value() != count:
            self.pane_spin.setValue(count)
        self.settings["pane_count"] = count
        self._save()
        if hasattr(self, "calc_window"):
            self.calc_window.update_pane_choices(count)

    def _on_sync_toggled(self, on):
        self.settings["sync_panes"] = on
        self._save()
        filled = [p for p in self.panes if p.item is not None]
        if on and filled:
            self._on_pane_x_changed(filled[0], *filled[0].ax.get_xlim())

    def _on_pane_x_changed(self, source, lo, hi):
        if self._syncing or not self.settings["sync_panes"] or source.item is None:
            return
        self._syncing = True
        try:
            for pane in self.panes:
                if pane is not source and pane.item is not None:
                    pane.set_view((lo, hi), remember=False)
        finally:
            self._syncing = False

    def show_in_pane(self, index, item):
        if not 0 <= index < len(self.panes):
            return
        pane = self.panes[index]
        pane.show_item(item)
        others = [p for p in self.panes if p is not pane and p.item is not None]
        if self.settings["sync_panes"] and others:
            self._syncing = True
            try:
                pane.set_view(others[0].ax.get_xlim(), remember=False)
            finally:
                self._syncing = False
        self._status(f"枠 {index + 1} に表示しました: {item['name']}")

    def first_empty_pane(self):
        for i, pane in enumerate(self.panes):
            if pane.item is None:
                return i
        return None

    def copy_spectrum_to_pane(self, index=None):
        if self.spectrum is None:
            return
        if index is None:
            index = self.first_empty_pane()
            if index is None:
                self.ctx.show_message("空いている枠がありません。枠の数を増やすか、枠を空にしてください。", "MS パック")
                return
        s = self.spectrum
        mz, y = compact_zeros(s["mz"], s["y"])
        self.show_in_pane(index, {"kind": "measured", "name": s["name"], "mz": mz, "y": y, "full_mz": s["mz"],
                                  "full_y": s["y"], "peaks": list(s["peaks"]), "run": s["run"],
                                  "provenance": dict(s["provenance"])})

    def clear_all_panes(self):
        for pane in self.panes:
            pane.show_item(None)

    def reset_all_views(self):
        for plot in [self.tic_plot, self.spectrum_plot, *self.panes]:
            plot.reset_view()

    # ================================================================ 開く
    def _choose_mzml(self):
        path, _ = QFileDialog.getOpenFileName(self, "mzML を開く", "", "mzML (*.mzML *.mzml);;すべて (*)")
        if path:
            self.open_mzml(path)

    def _choose_d_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Bruker の .d フォルダを選ぶ")
        if path:
            self.open_d(path)

    def open_mzml(self, path, display_name=None):
        QGuiApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            started = time.perf_counter()
            run = read_mzml(path)
            elapsed = time.perf_counter() - started
        except (MzmlError, OSError) as e:
            self.ctx.show_error(f"{os.path.basename(path)} を読めませんでした。\n{e}", "MS パック")
            return None
        finally:
            QGuiApplication.restoreOverrideCursor()
        if display_name:
            run.display_name = display_name
        used = {getattr(r, "color", None) for r in self.runs}
        run.color = next((c for c in RUN_COLORS if c not in used), RUN_COLORS[len(self.runs) % len(RUN_COLORS)])
        self.runs.append(run)
        item = QListWidgetItem(_color_icon(run.color), self._run_name(run))
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(Qt.CheckState.Checked)
        item.setToolTip(f"{run.path}\n{self._polarity_text(run)}、{len(run.ms1_scans())} スキャン")
        self.run_list.blockSignals(True)
        self.run_list.addItem(item)
        self.run_list.blockSignals(False)
        self.run_list.setCurrentRow(len(self.runs) - 1)
        self._status(f"{self._run_name(run)} を読み込みました({elapsed:.1f} 秒)")
        return run

    @staticmethod
    def _polarity_text(run):
        return {1: "正イオン (+)", -1: "負イオン (−)"}.get(run.polarity(), "極性混在")

    def _tic_visible(self, index):
        item = self.run_list.item(index)
        return item is not None and item.checkState() == Qt.CheckState.Checked

    def redraw_tics(self, keep_view=True):
        """チェックした測定の TIC を重ね、フォーカスしている測定を太く描く。範囲の帯はフォーカスの測定のもの。"""
        current = self.current_run()
        traces = [{"times": run.times(), "tics": run.tics(), "color": run.color, "focused": run is current}
                  for i, run in enumerate(self.runs) if self._tic_visible(i)]
        self.tic_plot.set_traces(traces, keep_view=keep_view)
        if current is not None:
            self._show_ranges(self.ranges[current.path])
            sign = {1: "+", -1: "−"}.get(current.polarity(), "±")
            self.tic_header.setText(f"クロマトグラム - {self._run_name(current)}: TIC {sign}")
        else:
            self.tic_header.setText("クロマトグラム")

    def _on_run_item_changed(self, _item):
        self.redraw_tics()

    @staticmethod
    def _run_name(run):
        return getattr(run, "display_name", None) or run.name

    def open_d(self, d_path):
        if not convert.is_d_folder(d_path):
            self.ctx.show_error(f"{d_path} は .d フォルダではありません。", "MS パック")
            return
        if self._process is not None:
            self.ctx.show_error("別の測定を変換中です。終わってから開いてください。", "MS パック")
            return
        msconvert = convert.find_msconvert(self.settings["msconvert_path"])
        if msconvert is None:
            msconvert = self.ask_msconvert()
            if msconvert is None:
                return
        cache_dir = os.path.join(self.ctx.data_dir, "mzml_cache")
        os.makedirs(cache_dir, exist_ok=True)
        name = os.path.basename(d_path.rstrip("\\/"))[:-2]
        final = convert.cached_mzml_path(d_path, cache_dir, msconvert)
        if os.path.isfile(final):
            os.utime(final)
            self.open_mzml(final, display_name=name)
            return
        self._start_conversion(msconvert, d_path, final, name)

    def ask_msconvert(self, explain=True):
        if explain:
            QMessageBox.information(self, "MS パック", MSCONVERT_MISSING)
        path, _ = QFileDialog.getOpenFileName(self, "msconvert.exe を選ぶ", "", "msconvert (msconvert.exe);;すべて (*)")
        if not path:
            return None
        self.settings["msconvert_path"] = path
        self._save()
        return path

    def _start_conversion(self, msconvert, d_path, final, name):
        work = convert.work_dir_for(final)
        os.makedirs(work, exist_ok=True)
        program, args = convert.conversion_command(msconvert, d_path, work)
        process = QProcess(self)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.finished.connect(self._on_conversion_finished)
        process.errorOccurred.connect(self._on_conversion_error)
        self._process = process
        self._pending = (d_path, final, name)
        self._convert_started = time.monotonic()
        self._set_converting(True)
        process.start(program, args)
        self._elapsed_timer.start(500)
        self._show_conversion_progress()

    def _set_converting(self, converting):
        self.progress_row.setVisible(converting)
        self.cancel_button.setVisible(converting)
        self.cancel_action.setEnabled(converting)
        self.open_d_action.setEnabled(not converting)

    def _show_conversion_progress(self):
        if self._pending:
            elapsed = time.monotonic() - self._convert_started
            self.progress_label.setText(f"{self._pending[2]} を mzML に変換中… {elapsed:.0f} 秒")

    def _end_conversion(self):
        self._elapsed_timer.stop()
        process, self._process = self._process, None
        pending, self._pending = self._pending, None
        self._set_converting(False)
        output = ""
        if process is not None:
            output = bytes(process.readAll()).decode("utf-8", "replace")
            process.deleteLater()
        return pending, output

    def _on_conversion_finished(self, exit_code, _status):
        pending, output = self._end_conversion()
        if pending is None:
            return
        _d_path, final, name = pending
        result = convert.finish_conversion(final) if exit_code == 0 else None
        if result is None:
            convert.finish_conversion(final)  # 作業フォルダを片付ける
            tail = "\n".join(output.strip().splitlines()[-8:])
            self.ctx.show_error(f"{name}.d を変換できませんでした(終了コード {exit_code})。\n\n{tail}", "MS パック")
            self._status("変換に失敗しました")
            return
        convert.prune_cache(os.path.dirname(final))
        self.open_mzml(final, display_name=name)

    def _on_conversion_error(self, error):
        if error == QProcess.ProcessError.FailedToStart:
            pending, _ = self._end_conversion()
            if pending:
                convert.finish_conversion(pending[1])
            self.ctx.show_error("msconvert を起動できませんでした。ProteoWizard のインストールを確かめてください。",
                                "MS パック")
            self._status("変換に失敗しました")

    def cancel_conversion(self):
        if self._process is None:
            return
        process = self._process
        pending = self._pending
        process.finished.disconnect(self._on_conversion_finished)
        process.kill()
        process.waitForFinished(3000)
        self._end_conversion()
        if pending:
            convert.finish_conversion(pending[1])
        self._status("変換を中止しました")

    def close_current_run(self):
        i = self.run_list.currentRow()
        if i < 0:
            return
        run = self.runs.pop(i)
        self.ranges.pop(run.path, None)
        self.run_list.blockSignals(True)
        self.run_list.takeItem(i)
        self.run_list.blockSignals(False)
        if self.runs:
            self.run_list.setCurrentRow(min(i, len(self.runs) - 1))
            self._on_run_selected(self.run_list.currentRow())
        else:
            self.redraw_tics(keep_view=False)
            self.spectrum_plot.clear()
            self.spectrum = None
            self.spectrum_header.setText("スペクトル")
            self._status("測定を開いてください")

    # ================================================================ 範囲とスペクトル
    def current_run(self):
        i = self.run_list.currentRow()
        return self.runs[i] if 0 <= i < len(self.runs) else None

    def _on_run_selected(self, index):
        run = self.current_run()
        if run is None:
            return
        # 開いた直後は全範囲を試料にし、背景はなし。同じ測定に戻ったときは前の範囲を使う
        self.ranges.setdefault(run.path, {"sample": self._full_range(run), "background": None})
        self.redraw_tics(keep_view=len(self.runs) > 1)
        self.compute_spectrum(keep_view=False)

    @staticmethod
    def _full_range(run):
        times = run.times()
        return (float(times[0]), float(times[-1])) if len(times) else (0.0, 0.0)

    def reset_sample_range(self):
        run = self.current_run()
        if run is None:
            return
        self.ranges[run.path]["sample"] = self._full_range(run)
        self._show_ranges(self.ranges[run.path])
        self.compute_spectrum()

    def _show_ranges(self, ranges):
        self._updating = True
        try:
            for (lo, hi), kind in ((ranges["sample"], "sample"), (ranges["background"] or (0.0, 0.0), "background")):
                w = self.range_window
                boxes = (w.sample_from, w.sample_to) if kind == "sample" else (w.bg_from, w.bg_to)
                boxes[0].setValue(lo)
                boxes[1].setValue(hi)
            self.tic_plot.set_range("sample", *ranges["sample"])
            bg = ranges["background"]
            self.tic_plot.set_range("background", *(bg if bg else (None, None)))
        finally:
            self._updating = False

    def _on_tic_range_changing(self, kind, t0, t1):
        self._updating = True
        try:
            w = self.range_window
            boxes = (w.sample_from, w.sample_to) if kind == "sample" else (w.bg_from, w.bg_to)
            boxes[0].setValue(t0)
            boxes[1].setValue(t1)
        finally:
            self._updating = False

    def _on_tic_range_changed(self, kind, t0, t1):
        run = self.current_run()
        if run is None:
            return
        self.ranges[run.path][kind] = (t0, t1)
        self._on_tic_range_changing(kind, t0, t1)
        self.compute_spectrum()

    def _on_scan_clicked(self, t):
        run = self.current_run()
        scan = run.nearest_scan(t) if run else None
        if scan is None:
            return
        self.ranges[run.path]["sample"] = (scan.time_min, scan.time_min)
        self._show_ranges(self.ranges[run.path])
        self.compute_spectrum()

    def _on_range_spin_changed(self, *_):
        if self._updating:
            return
        run = self.current_run()
        if run is None:
            return
        ranges = self.ranges[run.path]
        ranges["sample"] = (self.range_window.sample_from.value(), self.range_window.sample_to.value())
        bg = (self.range_window.bg_from.value(), self.range_window.bg_to.value())
        ranges["background"] = bg if bg[1] > bg[0] else None
        self.tic_plot.set_range("sample", *ranges["sample"])
        self.tic_plot.set_range("background", *(ranges["background"] or (None, None)))
        self._recompute_timer.start(RECOMPUTE_DELAY_MS)

    def clear_background(self):
        run = self.current_run()
        if run is None:
            return
        self.ranges[run.path]["background"] = None
        self._show_ranges(self.ranges[run.path])
        self.compute_spectrum()

    def compute_spectrum(self, keep_view=True):
        """試料の範囲の平均(背景を引く設定なら背景の平均を引く)を作り、2段目に表示する。"""
        self._recompute_timer.stop()
        run = self.current_run()
        if run is None:
            return
        ranges = self.ranges[run.path]
        scans = run.scans_in_range(*ranges["sample"])
        if not scans:
            self.spectrum = None
            self.spectrum_plot.clear()
            self.spectrum_header.setText(f"スペクトル - {self._run_name(run)}(選んだ時間範囲にスキャンがありません)")
            self._status("選んだ時間範囲にスキャンがありません")
            return
        QGuiApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            mz, y = average_spectrum(scans)
            bg_scans = []
            if self.settings["subtract_background"] and ranges["background"]:
                bg_scans = run.scans_in_range(*ranges["background"])
            if bg_scans:
                bg_mz, bg_y = average_spectrum(bg_scans)
                y = subtract_background(mz, y, bg_mz, bg_y, clip_negative=self.settings["clip_negative"])
            top = float(np.nanmax(y)) if len(y) else 0.0
            peaks = find_peaks(mz, y, min_height=top * PEAK_FLOOR) if top > 0 else []
        finally:
            QGuiApplication.restoreOverrideCursor()
        name = f"{self._run_name(run)} {range_text(*ranges['sample'])}"
        provenance = {"plugin": "mass_spec", "source": run.path, "sample_range_min": list(ranges["sample"]),
                      "scans": len(scans)}
        if bg_scans:
            name += f" − 背景 {range_text(*ranges['background'])}"
            provenance["background_range_min"] = list(ranges["background"])
            provenance["background_scans"] = len(bg_scans)
            provenance["negative_clipped"] = self.settings["clip_negative"]
        self.spectrum = {"mz": mz, "y": y, "peaks": peaks, "name": name, "provenance": provenance, "run": run}
        shown_mz, shown_y = compact_zeros(mz, y)
        self.spectrum_plot.set_spectrum(shown_mz, shown_y, peaks, keep_view=keep_view)
        average = f"{len(scans)} スキャンの平均" if len(scans) > 1 else "1 スキャン"
        self.spectrum_header.setText(f"スペクトル - {name}({average})")
        self._status(f"{name}({len(scans)} スキャンの平均"
                     + (f"、背景 {len(bg_scans)} スキャン" if bg_scans else "") + ")")
        if self.calc is not None:
            self.run_calculation(quiet=True)

    # ================================================================ 計算
    def _adducts(self):
        preset = self.calc_window.preset_combo.currentIndex()
        adducts = list(POSITIVE_STANDARD if preset == 0 else NEGATIVE_STANDARD if preset == 1 else ())
        extra = split_adducts(self.calc_window.extra_edit.text())
        return adducts + [a for a in extra if a not in adducts]

    def run_calculation(self, quiet=False):
        w = self.calc_window
        self.settings.update(formula=w.formula_edit.text().strip(), preset=w.preset_combo.currentIndex(),
                             extra_adducts=w.extra_edit.text().strip(), resolution=w.resolution_spin.value(),
                             tolerance_ppm=w.tolerance_spin.value())
        self._save()
        try:
            counts = parse_formula(self.settings["formula"])
            adducts = self._adducts()
            if not adducts:
                raise FormulaError("付加イオンがありません。標準のセットを選ぶか、表記を入力してください")
            for a in adducts:
                parse_adduct(a)
        except FormulaError as e:
            if not quiet:
                self.ctx.show_error(str(e), "MS パック")
            return
        if self.spectrum is not None:
            results = match_all(counts, adducts, self.spectrum["mz"], self.spectrum["y"],
                                resolution=self.settings["resolution"], tolerance_ppm=self.settings["tolerance_ppm"],
                                min_detect_percent=self.settings["detect_percent"])
        else:
            results = []
            for a in adducts:
                try:
                    results.append(ion_pattern(counts, parse_adduct(a)))
                except FormulaError as e:
                    results.append(e)
        self.calc = {"counts": counts, "formula": format_formula(counts), "results": results, "adducts": adducts}
        target = f"照合の対象: {self.spectrum['name']}" if self.spectrum else "実測なし(計算だけ)"
        w.summary.setText(
            f"{format_formula(counts)}  モノアイソトピック質量 {monoisotopic_mass(counts):.{self.settings['mz_decimals']}f}"
            f"  平均分子量 {average_mass(counts):.4f}\n{target}")
        self._show_calc_table()
        self.show_overlays()

    def _overlay_colors(self):
        cycle = self.ctx.active_color_cycle() or ["#d62728"]
        return cycle[1:] + cycle[:1] if len(cycle) > 1 else cycle

    def show_overlays(self):
        overlays = []
        if self.calc is not None and self.spectrum is not None and self.calc_window.overlay_check.isChecked():
            colors = self._overlay_colors()
            for i, r in enumerate(self._found_results()):
                sx, sy = r.sticks()
                overlays.append((sx, sy, colors[i % len(colors)], r.pattern.adduct.notation, "stick"))
        self.spectrum_plot.set_overlays(overlays)

    def _found_results(self):
        if self.calc is None:
            return []
        return [r for r in self.calc["results"] if not isinstance(r, Exception) and getattr(r, "found", False)]

    def _pattern_item(self, adduct):
        """付加イオンの計算パターンを3段目の枠の中身にする。実測と照合したなら実測の強度に合わせる。"""
        if self.calc is None:
            return None
        colors = self._overlay_colors()
        for i, (notation, r) in enumerate(zip(self.calc["adducts"], self.calc["results"])):
            if notation != adduct or isinstance(r, Exception):
                continue
            color = colors[i % len(colors)]
            if hasattr(r, "found"):
                pattern, resolution = r.pattern, r.resolution
            else:
                pattern, resolution = r, self.settings["resolution"] or DEFAULT_RESOLUTION
            if getattr(r, "found", False):
                px, py = r.profile()
                sticks = r.sticks()
            else:
                px, py = gaussian_profile(pattern.mz, pattern.relative, resolution)
                sticks = (pattern.mz, pattern.relative)
            return {"kind": "calc", "name": f"{self.calc['formula']} {notation} 計算(R {resolution:.0f})",
                    "mz": px, "y": py, "sticks": sticks, "peaks": [], "color": color, "result": r,
                    "pattern": pattern, "resolution": resolution}
        return None

    def pattern_to_pane(self, adduct, index):
        item = self._pattern_item(adduct)
        if item is None:
            self.ctx.show_message(f"{adduct} の計算結果がありません。", "MS パック")
            return
        self.show_in_pane(index, item)

    def found_patterns_to_panes(self):
        if self.calc is None:
            return
        chosen = [r.pattern.adduct.notation for r in self._found_results()]
        if not chosen and self.spectrum is None:
            chosen = [a for a, r in zip(self.calc["adducts"], self.calc["results"]) if not isinstance(r, Exception)]
        if not chosen:
            self.ctx.show_message("実測に見つかった付加イオンがありません。", "MS パック")
            return
        if len(chosen) > len(self.panes):
            self.set_pane_count(min(MAX_PANES, len(chosen)))
        for i, adduct in enumerate(chosen[:len(self.panes)]):
            self.pattern_to_pane(adduct, i)

    def calc_rows(self):
        """表の行(表示用の文字列)。"""
        fmt = self.label_format()
        header = ["付加イオン", "ピーク", "計算 m/z", "実測 m/z", "誤差 (ppm)", "計算 %", "実測 %", "R", "状態"]
        rows = []
        if self.calc is None:
            return header, rows
        for adduct, r in zip(self.calc["adducts"], self.calc["results"]):
            if isinstance(r, Exception):
                rows.append([adduct, "", "", "", "", "", "", "", f"エラー: {r}"])
                continue
            if not hasattr(r, "found"):  # 実測なし: IonPattern
                mono = int(np.argmin(np.abs(r.mz - r.monoisotopic_mz)))
                for k, (m, rel) in enumerate(zip(r.mz, r.relative)):
                    if rel < 1.0:
                        continue
                    shift = int(r.offsets[k] - r.offsets[mono])
                    rows.append([adduct, "M" if shift == 0 else f"M{shift:+d}", fmt.mz(m), "", "",
                                 f"{rel:.1f}", "", "", "計算のみ"])
                continue
            for p in r.peaks:
                measured = np.isfinite(p.measured_mz)
                rows.append([
                    adduct, p.label, fmt.mz(p.calc_mz), fmt.mz(p.measured_mz) if measured else "",
                    f"{p.error_ppm:+.{fmt.ppm_decimals}f}" if measured else "", f"{p.calc_relative:.1f}",
                    f"{p.measured_relative:.1f}" if measured else "", f"{r.resolution:.0f}({r.resolution_source})",
                    r.status])
        return header, rows

    def _show_calc_table(self):
        self.calc_window.show_table(*self.calc_rows())

    def copy_table(self):
        header, rows = self.calc_rows()
        if not rows:
            return
        QGuiApplication.clipboard().setText("\n".join("\t".join(r) for r in [header, *rows]))
        self._status("照合表をコピーしました")

    # ================================================================ 転送
    def _next_color(self):
        cycle = self.ctx.active_color_cycle() or ["#1f77b4"]
        return cycle[len(self.ctx.datasets()) % len(cycle)]

    def _add(self, dataset, what):
        dataset.subplot_target = self.settings["subplot_target"]
        self.ctx.add_dataset(dataset, description=f"[MS] {what}を追加")
        self._status(f"追加しました: {dataset.name}")

    def transfer_tic(self):
        run = self.current_run()
        if run is None:
            return
        self._add(tic_dataset(run, self._next_color(), {"plugin": "mass_spec", "source": run.path}), "TIC ")

    def _transfer_measured(self, source, plot, view_only=None):
        if view_only is None:
            view_only = self.settings["transfer_view_only"]
        x_range = plot.ax.get_xlim() if view_only else None
        name = source["name"] + (f" (m/z {x_range[0]:.0f}–{x_range[1]:.0f})" if x_range else "")
        mz = source.get("full_mz", source["mz"])
        y = source.get("full_y", source["y"])
        ds = spectrum_dataset(name, mz, y, self._next_color(), source["run"].path, dict(source["provenance"]),
                              x_range=x_range)
        self._add(ds, "スペクトル")

    def transfer_spectrum(self, view_only=None):
        if self.spectrum is not None:
            self._transfer_measured(self.spectrum, self.spectrum_plot, view_only)

    def _transfer_peaks(self, source, plot, labels_only=False):
        peaks = source["peaks"]
        if self.settings["transfer_view_only"]:
            lo, hi = plot.ax.get_xlim()
            peaks = [p for p in peaks if lo <= p.mz <= hi]
        if not peaks:
            return
        mode = self.settings["label_mode"]
        ds = centroid_dataset(
            source["name"] + (" ラベル" if labels_only else " centroid"), peaks, self._next_color(),
            self.label_format(),
            label_top_n=self.settings["label_top_n"] if mode == "top" else (None if mode == "percent" else 0),
            label_min_relative=self.settings["label_percent"] if mode == "percent" else None,
            min_relative=self.settings["centroid_min_percent"], source_file=source["run"].path,
            provenance=dict(source["provenance"], kind="centroid"), labels_only=labels_only, pinned=plot.pinned)
        self._add(ds, "ピークのラベル" if labels_only else "centroid ")

    def transfer_centroid(self, labels_only=False):
        if self.spectrum is not None:
            self._transfer_peaks(self.spectrum, self.spectrum_plot, labels_only)

    def _calc_datasets(self, item):
        fmt = self.label_format()
        r = item["result"]
        provenance = {"plugin": "mass_spec", "formula": self.calc["formula"] if self.calc else "",
                      "adduct": item["pattern"].adduct.notation, "resolution": item["resolution"]}
        if hasattr(r, "found") and r.found:
            return pattern_datasets(r, self.calc["formula"], item["color"], fmt, provenance=provenance)
        sx, sy = item["sticks"]
        labels = [fmt.mz(m) if rel >= 5.0 else "" for m, rel in zip(sx, item["pattern"].relative)]
        return [stick_dataset(item["name"], sx, sy, labels, item["color"], provenance=provenance)]

    def transfer_calculation(self):
        if self.calc is None:
            return
        results = self._found_results() if self.spectrum is not None else [
            r for r in self.calc["results"] if not isinstance(r, Exception)]
        if not results:
            self.ctx.show_message("実測に見つかった付加イオンがないので、転送するものがありません。", "MS パック")
            return
        for r in results:
            notation = r.pattern.adduct.notation if hasattr(r, "found") else r.adduct.notation
            for ds in self._calc_datasets(self._pattern_item(notation)):
                self._add(ds, "計算パターン")

    def transfer_pane(self, pane):
        item = pane.item
        if item is None:
            return
        if item["kind"] == "measured":
            self._transfer_measured(item, pane)
        else:
            for ds in self._calc_datasets(item):
                self._add(ds, "計算パターン")

    # ================================================================ 右クリック
    def _tic_menu(self, pos):
        menu = QMenu(self)
        menu.addAction("TIC をプロットに転送", self.transfer_tic)
        menu.addAction("背景の範囲を解除", self.clear_background)
        menu.addAction("全体を表示", self.tic_plot.reset_view)
        menu.exec(pos)

    def _spectrum_menu(self, pos):
        menu = QMenu(self)
        copy_menu = menu.addMenu("3段目の枠にコピー")
        copy_menu.addAction("空いている枠", self.copy_spectrum_to_pane)
        for i in range(len(self.panes)):
            copy_menu.addAction(f"枠 {i + 1}", lambda i=i: self.copy_spectrum_to_pane(i))
        menu.addSeparator()
        menu.addAction("表示中の範囲をプロットに転送", lambda: self.transfer_spectrum(view_only=True))
        menu.addAction("全範囲をプロットに転送", lambda: self.transfer_spectrum(view_only=False))
        menu.addAction("centroid を転送", self.transfer_centroid)
        menu.addAction("ラベルだけを転送", lambda: self.transfer_centroid(labels_only=True))
        menu.addAction("計算パターンを転送", self.transfer_calculation)
        menu.addSeparator()
        menu.addAction("表示中のピーク一覧をコピー", lambda: self.copy_visible_peaks(self.spectrum_plot))
        menu.addAction("全体を表示", self.spectrum_plot.reset_view)
        menu.exec(pos)

    def _pane_menu(self, pane, pos):
        menu = QMenu(self)
        if pane.item is not None:
            menu.addAction("プロットに転送", lambda: self.transfer_pane(pane))
            if pane.item["kind"] == "measured":
                menu.addAction("centroid を転送", lambda: self._transfer_peaks(pane.item, pane))
            menu.addAction("表示中のピーク一覧をコピー", lambda: self.copy_visible_peaks(pane))
            menu.addAction("枠を空にする", lambda: pane.show_item(None))
            menu.addAction("全体を表示", pane.reset_view)
        else:
            menu.addAction("2段目のスペクトルをここにコピー",
                           lambda: self.copy_spectrum_to_pane(self.panes.index(pane)))
        menu.exec(pos)

    def copy_visible_peaks(self, plot):
        x0, x1 = plot.ax.get_xlim()
        fmt = self.label_format()
        lines = ["m/z\t強度"] + [f"{fmt.mz(p.mz)}\t{p.height:.6g}" for p in plot.peaks if x0 <= p.mz <= x1]
        QGuiApplication.clipboard().setText("\n".join(lines))
        self._status(f"ピーク {len(lines) - 1} 本をコピーしました")

    # ================================================================ 後片付け
    def closeEvent(self, event):
        if self._process is not None:
            self.cancel_conversion()
        super().closeEvent(event)
