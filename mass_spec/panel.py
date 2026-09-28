"""ドックパネル「MS スペクトル」。測定を開き、TIC から範囲を選んで平均し、ラベル・照合・本体への転送を行う。"""
import os
import time

import numpy as np
from PySide6.QtCore import QProcess, Qt, QTimer
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMenu, QMessageBox, QPushButton, QScrollArea, QSpinBox,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from . import convert
from .analyzer import PRESET_NEGATIVE, PRESET_NONE, PRESET_POSITIVE
from .chemistry import (
    NEGATIVE_STANDARD, POSITIVE_STANDARD, FormulaError, average_mass, format_formula, ion_pattern,
    monoisotopic_mass, parse_adduct, parse_formula, split_adducts,
)
from .datasets import centroid_dataset, pattern_datasets, range_text, spectrum_dataset, stick_dataset, tic_dataset
from .matching import DEFAULT_RESOLUTION, match_all
from .mzml import MzmlError, read_mzml
from .plots import SpectrumPlot, TicPlot
from .settings import load_settings, save_settings
from .spectra import LabelFormat, average_spectrum, compact_zeros, find_peaks, subtract_background

PANEL_NAME = "MS スペクトル"
PRESETS = [PRESET_POSITIVE, PRESET_NEGATIVE, PRESET_NONE]
LABEL_MODES = [("top", "強い順に N 本"), ("percent", "相対強度 % 以上"), ("none", "表示しない")]
INTENSITY_MODES = [("none", "なし"), ("relative", "相対 %"), ("absolute", "絶対値")]
# ノイズの極大まで拾うとピークが数万になるので、最大の 0.1% 未満は拾わない
PEAK_FLOOR = 1e-3
RECOMPUTE_DELAY_MS = 250

MSCONVERT_MISSING = (
    "ProteoWizard の msconvert が見つかりません。\n\n"
    "https://proteowizard.sourceforge.io/download.html から Windows 64-bit 版"
    "(ベンダー形式を読める版)をインストールしてください。\n"
    "別の場所に入れた場合は、次の画面で msconvert.exe を選べます。")


def _spin(minimum, maximum, value, decimals=3, step=0.1, width=80):
    box = QDoubleSpinBox()
    box.setDecimals(decimals)
    box.setRange(minimum, maximum)
    box.setSingleStep(step)
    box.setValue(value)
    box.setMaximumWidth(width)
    box.setKeyboardTracking(False)
    return box


class MassSpecPanel(QWidget):
    def __init__(self, ctx, parent=None):
        super().__init__(parent)
        self.ctx = ctx
        self.settings = load_settings(ctx.data_dir)
        self.runs = []
        self.ranges = {}            # run.path -> {"sample": (t0, t1), "background": (t0, t1) | None}
        self.spectrum = None        # {"mz", "y", "peaks", "name", "provenance"}
        self.calc = None            # {"counts", "formula", "results"}
        self._process = None
        self._pending = None        # 変換中の (d のパス, キャッシュの mzML)
        self._convert_started = 0.0
        self._updating = False
        self._recompute_timer = QTimer(self)
        self._recompute_timer.setSingleShot(True)
        self._recompute_timer.timeout.connect(self.compute_spectrum)
        self._elapsed_timer = QTimer(self)
        self._elapsed_timer.timeout.connect(self._show_conversion_progress)
        self._build()
        self._apply_label_settings()

    # ================================================================ 画面
    def _build(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll)
        body = QWidget()
        scroll.setWidget(body)
        layout = QVBoxLayout(body)

        row = QHBoxLayout()
        self.run_combo = QComboBox()
        self.run_combo.setMinimumContentsLength(12)
        self.run_combo.currentIndexChanged.connect(self._on_run_selected)
        self.open_d_button = QPushButton(".d を開く…")
        self.open_d_button.clicked.connect(self._choose_d_folder)
        self.open_mzml_button = QPushButton("mzML…")
        self.open_mzml_button.clicked.connect(self._choose_mzml)
        self.close_button = QPushButton("閉じる")
        self.close_button.clicked.connect(self.close_current_run)
        self.cancel_button = QPushButton("変換を中止")
        self.cancel_button.clicked.connect(self.cancel_conversion)
        self.cancel_button.hide()
        row.addWidget(self.run_combo, 1)
        for w in (self.open_d_button, self.open_mzml_button, self.close_button, self.cancel_button):
            row.addWidget(w)
        layout.addLayout(row)
        self.status_label = QLabel("測定を開いてください(.d は ProteoWizard の msconvert で自動的に変換します)")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        self.tic_plot = TicPlot()
        self.tic_plot.range_changing.connect(self._on_tic_range_changing)
        self.tic_plot.range_changed.connect(self._on_tic_range_changed)
        self.tic_plot.scan_clicked.connect(self._on_scan_clicked)
        self.tic_plot.context_requested.connect(self._tic_menu)
        layout.addWidget(self.tic_plot)

        row = QHBoxLayout()
        self.sample_from = _spin(0, 1e4, 0)
        self.sample_to = _spin(0, 1e4, 0)
        self.bg_from = _spin(0, 1e4, 0)
        self.bg_to = _spin(0, 1e4, 0)
        for box in (self.sample_from, self.sample_to, self.bg_from, self.bg_to):
            box.valueChanged.connect(self._on_range_spin_changed)
        self.subtract_check = QCheckBox("背景を引く")
        self.subtract_check.setChecked(self.settings["subtract_background"])
        self.subtract_check.toggled.connect(self._on_background_option)
        self.clip_check = QCheckBox("負は 0")
        self.clip_check.setChecked(self.settings["clip_negative"])
        self.clip_check.toggled.connect(self._on_background_option)
        row.addWidget(QLabel("試料"))
        row.addWidget(self.sample_from)
        row.addWidget(QLabel("–"))
        row.addWidget(self.sample_to)
        row.addWidget(QLabel("min  背景"))
        row.addWidget(self.bg_from)
        row.addWidget(QLabel("–"))
        row.addWidget(self.bg_to)
        row.addStretch(1)
        layout.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(self.subtract_check)
        row.addWidget(self.clip_check)
        self.clear_bg_button = QPushButton("背景の範囲を解除")
        self.clear_bg_button.clicked.connect(self.clear_background)
        row.addWidget(self.clear_bg_button)
        row.addStretch(1)
        layout.addLayout(row)

        self.spectrum_plot = SpectrumPlot()
        self.spectrum_plot.hover_text.connect(self._show_readout)
        self.spectrum_plot.measured.connect(self._show_readout)
        self.spectrum_plot.context_requested.connect(self._spectrum_menu)
        layout.addWidget(self.spectrum_plot)
        self.readout_label = QLabel(" ")
        self.readout_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.readout_label)

        row = QHBoxLayout()
        self.label_mode_combo = QComboBox()
        for _key, text in LABEL_MODES:
            self.label_mode_combo.addItem(text)
        self.label_mode_combo.setCurrentIndex([k for k, _ in LABEL_MODES].index(self.settings["label_mode"]))
        self.label_n_spin = QSpinBox()
        self.label_n_spin.setRange(1, 200)
        self.label_n_spin.setValue(self.settings["label_top_n"])
        self.label_percent_spin = _spin(0.0, 100.0, self.settings["label_percent"], decimals=1, step=1.0, width=70)
        self.mz_decimals_spin = QSpinBox()
        self.mz_decimals_spin.setRange(0, 6)
        self.mz_decimals_spin.setValue(self.settings["mz_decimals"])
        self.intensity_combo = QComboBox()
        for _key, text in INTENSITY_MODES:
            self.intensity_combo.addItem(text)
        self.intensity_combo.setCurrentIndex([k for k, _ in INTENSITY_MODES].index(self.settings["intensity"]))
        for w in (self.label_mode_combo, self.intensity_combo):
            w.currentIndexChanged.connect(self._on_label_settings_changed)
        for w in (self.label_n_spin, self.mz_decimals_spin, self.label_percent_spin):
            w.valueChanged.connect(self._on_label_settings_changed)
        row.addWidget(QLabel("ラベル"))
        row.addWidget(self.label_mode_combo)
        row.addWidget(self.label_n_spin)
        row.addWidget(self.label_percent_spin)
        row.addWidget(QLabel("m/z 桁"))
        row.addWidget(self.mz_decimals_spin)
        row.addWidget(QLabel("強度"))
        row.addWidget(self.intensity_combo)
        row.addStretch(1)
        layout.addLayout(row)

        layout.addWidget(self._build_calc_group())
        layout.addWidget(self._build_transfer_group())
        layout.addStretch(1)

    def _build_calc_group(self):
        group = QGroupBox("組成式から計算")
        form = QFormLayout(group)
        self.formula_edit = QLineEdit(self.settings["formula"])
        self.formula_edit.setPlaceholderText("例: C6H12O6")
        self.formula_edit.returnPressed.connect(self.run_calculation)
        form.addRow("組成式", self.formula_edit)
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(PRESETS)
        self.preset_combo.setCurrentIndex(min(max(self.settings["preset"], 0), len(PRESETS) - 1))
        form.addRow("付加イオン", self.preset_combo)
        self.extra_edit = QLineEdit(self.settings["extra_adducts"])
        self.extra_edit.setPlaceholderText("追加: [2M+Na]+, [M+2H]2+")
        self.extra_edit.returnPressed.connect(self.run_calculation)
        form.addRow("", self.extra_edit)
        row = QHBoxLayout()
        self.resolution_spin = _spin(0, 1e7, self.settings["resolution"], decimals=0, step=1000, width=90)
        self.resolution_spin.setSpecialValueText("実測から")
        self.tolerance_spin = _spin(1, 5000, self.settings["tolerance_ppm"], decimals=1, step=5, width=70)
        self.ppm_decimals_spin = QSpinBox()
        self.ppm_decimals_spin.setRange(0, 3)
        self.ppm_decimals_spin.setValue(self.settings["ppm_decimals"])
        self.ppm_decimals_spin.valueChanged.connect(self._on_label_settings_changed)
        row.addWidget(self.resolution_spin)
        row.addWidget(QLabel("±ppm"))
        row.addWidget(self.tolerance_spin)
        row.addWidget(QLabel("ppm 桁"))
        row.addWidget(self.ppm_decimals_spin)
        row.addStretch(1)
        form.addRow("分解能 R", row)
        row = QHBoxLayout()
        self.calc_button = QPushButton("計算して重ねる")
        self.calc_button.clicked.connect(self.run_calculation)
        self.copy_table_button = QPushButton("表をコピー")
        self.copy_table_button.clicked.connect(self.copy_table)
        row.addWidget(self.calc_button)
        row.addWidget(self.copy_table_button)
        row.addStretch(1)
        form.addRow(row)
        self.calc_summary = QLabel("")
        self.calc_summary.setWordWrap(True)
        self.calc_summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        form.addRow(self.calc_summary)
        self.table = QTableWidget(0, 0)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setMinimumHeight(140)
        form.addRow(self.table)
        return group

    def _build_transfer_group(self):
        group = QGroupBox("プロットに転送")
        layout = QVBoxLayout(group)
        row = QHBoxLayout()
        self.view_only_check = QCheckBox("表示中の m/z 範囲だけ")
        self.view_only_check.setChecked(self.settings["transfer_view_only"])
        self.view_only_check.toggled.connect(self._on_transfer_option)
        self.subplot_spin = QSpinBox()
        self.subplot_spin.setRange(1, 99)
        self.subplot_spin.setValue(self.settings["subplot_target"] + 1)
        self.subplot_spin.valueChanged.connect(self._on_transfer_option)
        row.addWidget(self.view_only_check)
        row.addWidget(QLabel("サブプロット"))
        row.addWidget(self.subplot_spin)
        row.addStretch(1)
        layout.addLayout(row)
        row = QHBoxLayout()
        self.send_tic_button = QPushButton("TIC")
        self.send_tic_button.clicked.connect(self.transfer_tic)
        self.send_spectrum_button = QPushButton("スペクトル")
        self.send_spectrum_button.clicked.connect(self.transfer_spectrum)
        self.send_centroid_button = QPushButton("centroid")
        self.send_centroid_button.clicked.connect(self.transfer_centroid)
        self.send_labels_button = QPushButton("ラベルだけ")
        self.send_labels_button.setToolTip("ピークのラベルを、棒を描かない Stick として送る(profile に重ねる用)")
        self.send_labels_button.clicked.connect(lambda: self.transfer_centroid(labels_only=True))
        self.send_calc_button = QPushButton("計算パターン")
        self.send_calc_button.clicked.connect(self.transfer_calculation)
        for w in (self.send_tic_button, self.send_spectrum_button, self.send_centroid_button,
                  self.send_labels_button, self.send_calc_button):
            row.addWidget(w)
        row.addStretch(1)
        layout.addLayout(row)
        return group

    # ================================================================ 設定
    def _save(self):
        save_settings(self.ctx.data_dir, self.settings)

    def label_format(self):
        return LabelFormat(mz_decimals=self.settings["mz_decimals"], intensity=self.settings["intensity"],
                           ppm_decimals=self.settings["ppm_decimals"])

    def _apply_label_settings(self):
        plot = self.spectrum_plot
        plot.label_mode = self.settings["label_mode"]
        plot.label_top_n = self.settings["label_top_n"]
        plot.label_percent = self.settings["label_percent"]
        plot.label_format = self.label_format()
        self.label_n_spin.setVisible(plot.label_mode == "top")
        self.label_percent_spin.setVisible(plot.label_mode == "percent")

    def _on_label_settings_changed(self, *_):
        self.settings.update(
            label_mode=LABEL_MODES[self.label_mode_combo.currentIndex()][0],
            label_top_n=self.label_n_spin.value(),
            label_percent=self.label_percent_spin.value(),
            mz_decimals=self.mz_decimals_spin.value(),
            intensity=INTENSITY_MODES[self.intensity_combo.currentIndex()][0],
            ppm_decimals=self.ppm_decimals_spin.value(),
        )
        self._apply_label_settings()
        self.spectrum_plot.update_labels()
        if self.calc is not None:
            self._show_calc_table()
        self._save()

    def _on_transfer_option(self, *_):
        self.settings["transfer_view_only"] = self.view_only_check.isChecked()
        self.settings["subplot_target"] = self.subplot_spin.value() - 1
        self._save()

    def _show_readout(self, text):
        self.readout_label.setText(text)

    def _status(self, text):
        self.status_label.setText(text)

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
        self.runs.append(run)
        polarity = {1: "+", -1: "−"}.get(run.polarity(), "±")
        self.run_combo.addItem(f"{self._run_name(run)}({polarity}、{len(run.ms1_scans())} スキャン)")
        self.run_combo.setCurrentIndex(len(self.runs) - 1)
        self._status(f"{self._run_name(run)} を読み込みました({elapsed:.1f} 秒)")
        return run

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
            msconvert = self._ask_msconvert()
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

    def _ask_msconvert(self):
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
        self.cancel_button.setVisible(converting)
        self.open_d_button.setEnabled(not converting)

    def _show_conversion_progress(self):
        if self._pending:
            elapsed = time.monotonic() - self._convert_started
            self._status(f"{self._pending[2]} を mzML に変換中… {elapsed:.0f} 秒")

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
        d_path, final, name = pending
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
        i = self.run_combo.currentIndex()
        if i < 0:
            return
        run = self.runs.pop(i)
        self.ranges.pop(run.path, None)
        self.run_combo.removeItem(i)
        if not self.runs:
            self.tic_plot.set_data([], [])
            self.spectrum_plot.clear()
            self.spectrum = None
            self._status("測定を開いてください")

    # ================================================================ 範囲とスペクトル
    def current_run(self):
        i = self.run_combo.currentIndex()
        return self.runs[i] if 0 <= i < len(self.runs) else None

    def _on_run_selected(self, index):
        run = self.current_run()
        if run is None:
            return
        times = run.times()
        self.tic_plot.set_data(times, run.tics())
        ranges = self.ranges.setdefault(run.path, {"sample": (float(times[0]), float(times[-1])) if len(times) else (0, 0),
                                                   "background": None})
        self._show_ranges(ranges)
        self.compute_spectrum(keep_view=False)

    def _show_ranges(self, ranges):
        self._updating = True
        try:
            for (lo, hi), kind in ((ranges["sample"], "sample"), (ranges["background"] or (0.0, 0.0), "background")):
                boxes = (self.sample_from, self.sample_to) if kind == "sample" else (self.bg_from, self.bg_to)
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
            boxes = (self.sample_from, self.sample_to) if kind == "sample" else (self.bg_from, self.bg_to)
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
        ranges["sample"] = (self.sample_from.value(), self.sample_to.value())
        bg = (self.bg_from.value(), self.bg_to.value())
        ranges["background"] = bg if bg[1] > bg[0] else None
        self.tic_plot.set_range("sample", *ranges["sample"])
        self.tic_plot.set_range("background", *(ranges["background"] or (None, None)))
        self._recompute_timer.start(RECOMPUTE_DELAY_MS)

    def _on_background_option(self, *_):
        self.settings["subtract_background"] = self.subtract_check.isChecked()
        self.settings["clip_negative"] = self.clip_check.isChecked()
        self._save()
        self.compute_spectrum()

    def clear_background(self):
        run = self.current_run()
        if run is None:
            return
        self.ranges[run.path]["background"] = None
        self._show_ranges(self.ranges[run.path])
        self.compute_spectrum()

    def compute_spectrum(self, keep_view=True):
        """試料の範囲の平均(背景を引く設定なら背景の平均を引く)を作り、パネルに表示する。"""
        self._recompute_timer.stop()
        run = self.current_run()
        if run is None:
            return
        ranges = self.ranges[run.path]
        scans = run.scans_in_range(*ranges["sample"])
        if not scans:
            self.spectrum = None
            self.spectrum_plot.clear()
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
        sample_text = range_text(*ranges["sample"])
        name = f"{self._run_name(run)} {sample_text}"
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
        self._status(f"{name}({len(scans)} スキャンの平均"
                     + (f"、背景 {len(bg_scans)} スキャン" if bg_scans else "") + ")")
        if self.calc is not None:
            self.run_calculation(quiet=True)

    # ================================================================ 計算
    def _adducts(self):
        preset = self.preset_combo.currentIndex()
        adducts = list(POSITIVE_STANDARD if preset == 0 else NEGATIVE_STANDARD if preset == 1 else ())
        extra = split_adducts(self.extra_edit.text())
        return adducts + [a for a in extra if a not in adducts]

    def run_calculation(self, quiet=False):
        self.settings.update(formula=self.formula_edit.text().strip(), preset=self.preset_combo.currentIndex(),
                             extra_adducts=self.extra_edit.text().strip(), resolution=self.resolution_spin.value(),
                             tolerance_ppm=self.tolerance_spin.value())
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
        self.calc_summary.setText(
            f"{format_formula(counts)}  モノアイソトピック質量 {monoisotopic_mass(counts):.{self.settings['mz_decimals']}f}"
            f"  平均分子量 {average_mass(counts):.4f}")
        self._show_calc_table()
        self._show_overlays()

    def _overlay_colors(self):
        cycle = self.ctx.active_color_cycle() or ["#d62728"]
        return cycle[1:] + cycle[:1] if len(cycle) > 1 else cycle

    def _show_overlays(self):
        overlays = []
        if self.calc is not None and self.spectrum is not None:
            colors = self._overlay_colors()
            found = [r for r in self.calc["results"] if not isinstance(r, Exception) and getattr(r, "found", False)]
            for i, r in enumerate(found):
                sx, sy = r.sticks()
                overlays.append((sx, sy, colors[i % len(colors)], r.pattern.adduct.notation, "stick"))
        self.spectrum_plot.set_overlays(overlays)

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
        header, rows = self.calc_rows()
        self.table.clear()
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                self.table.setItem(i, j, QTableWidgetItem(value))
        self.table.resizeColumnsToContents()

    def copy_table(self):
        header, rows = self.calc_rows()
        if not rows:
            return
        text = "\n".join("\t".join(r) for r in [header, *rows])
        QGuiApplication.clipboard().setText(text)
        self._status("照合表をコピーしました")

    # ================================================================ 転送
    def _next_color(self):
        cycle = self.ctx.active_color_cycle() or ["#1f77b4"]
        return cycle[len(self.ctx.datasets()) % len(cycle)]

    def _add(self, dataset, what):
        dataset.subplot_target = self.settings["subplot_target"]
        self.ctx.add_dataset(dataset, description=f"[MS] {what}を追加")
        self._status(f"追加しました: {dataset.name}")

    def _view_range(self):
        return self.spectrum_plot.ax.get_xlim() if self.settings["transfer_view_only"] else None

    def transfer_tic(self):
        run = self.current_run()
        if run is None:
            return
        self._add(tic_dataset(run, self._next_color(), {"plugin": "mass_spec", "source": run.path}), "TIC ")

    def transfer_spectrum(self, view_only=None):
        if self.spectrum is None:
            return
        s = self.spectrum
        x_range = self.spectrum_plot.ax.get_xlim() if view_only else (None if view_only is False else self._view_range())
        name = s["name"] + (f" (m/z {x_range[0]:.0f}–{x_range[1]:.0f})" if x_range else "")
        ds = spectrum_dataset(name, s["mz"], s["y"], self._next_color(), s["run"].path, dict(s["provenance"]),
                              x_range=x_range)
        self._add(ds, "スペクトル")

    def transfer_centroid(self, labels_only=False):
        if self.spectrum is None or not self.spectrum["peaks"]:
            return
        s = self.spectrum
        peaks = s["peaks"]
        x_range = self._view_range()
        if x_range:
            peaks = [p for p in peaks if x_range[0] <= p.mz <= x_range[1]]
        if not peaks:
            return
        mode = self.settings["label_mode"]
        ds = centroid_dataset(
            s["name"] + (" ラベル" if labels_only else " centroid"), peaks, self._next_color(), self.label_format(),
            label_top_n=self.settings["label_top_n"] if mode == "top" else (None if mode == "percent" else 0),
            label_min_relative=self.settings["label_percent"] if mode == "percent" else None,
            min_relative=self.settings["centroid_min_percent"], source_file=s["run"].path,
            provenance=dict(s["provenance"], kind="centroid"), labels_only=labels_only,
            pinned=self.spectrum_plot.pinned)
        self._add(ds, "ピークのラベル" if labels_only else "centroid ")

    def transfer_calculation(self):
        if self.calc is None:
            return
        fmt = self.label_format()
        colors = self._overlay_colors()
        added = 0
        for i, r in enumerate([r for r in self.calc["results"] if not isinstance(r, Exception)]):
            color = colors[i % len(colors)]
            notation = r.pattern.adduct.notation if hasattr(r, "found") else r.adduct.notation
            provenance = {"plugin": "mass_spec", "formula": self.calc["formula"], "adduct": notation}
            if hasattr(r, "found"):
                if not r.found:
                    continue
                for ds in pattern_datasets(r, self.calc["formula"], color, fmt,
                                           provenance=dict(provenance, resolution=r.resolution)):
                    self._add(ds, "計算パターン")
            else:
                resolution = self.settings["resolution"] or DEFAULT_RESOLUTION
                labels = [fmt.mz(m) if rel >= 5.0 else "" for m, rel in zip(r.mz, r.relative)]
                ds = stick_dataset(f"{self.calc['formula']} {r.adduct.notation} 計算", r.mz, r.relative, labels,
                                   color, provenance=dict(provenance, resolution=resolution))
                self._add(ds, "計算パターン")
            added += 1
        if not added:
            self.ctx.show_message("実測に見つかった付加イオンがないので、転送するものがありません。", "MS パック")

    # ================================================================ 右クリック
    def _tic_menu(self, pos):
        menu = QMenu(self)
        menu.addAction("TIC をプロットに転送", self.transfer_tic)
        menu.addAction("背景の範囲を解除", self.clear_background)
        menu.addAction("全体を表示", self.tic_plot.reset_view)
        menu.exec(pos)

    def _spectrum_menu(self, pos):
        menu = QMenu(self)
        menu.addAction("表示中の範囲を転送", lambda: self.transfer_spectrum(view_only=True))
        menu.addAction("全範囲を転送", lambda: self.transfer_spectrum(view_only=False))
        menu.addAction("centroid を転送", self.transfer_centroid)
        menu.addAction("ラベルだけを転送", lambda: self.transfer_centroid(labels_only=True))
        menu.addAction("計算パターンを転送", self.transfer_calculation)
        menu.addSeparator()
        toggle = menu.addAction("ラベルを表示")
        toggle.setCheckable(True)
        toggle.setChecked(self.settings["label_mode"] != "none")
        toggle.toggled.connect(self._toggle_labels)
        menu.addAction("表示中のピーク一覧をコピー", self.copy_visible_peaks)
        menu.addAction("全体を表示", self.spectrum_plot.reset_view)
        menu.exec(pos)

    def _toggle_labels(self, on):
        self.label_mode_combo.setCurrentIndex(0 if on else len(LABEL_MODES) - 1)

    def copy_visible_peaks(self):
        if self.spectrum is None:
            return
        x0, x1 = self.spectrum_plot.ax.get_xlim()
        fmt = self.label_format()
        lines = ["m/z\t強度"] + [f"{fmt.mz(p.mz)}\t{p.height:.6g}" for p in self.spectrum["peaks"] if x0 <= p.mz <= x1]
        QGuiApplication.clipboard().setText("\n".join(lines))
        self._status(f"ピーク {len(lines) - 1} 本をコピーしました")

    # ================================================================ 後片付け
    def closeEvent(self, event):
        if self._process is not None:
            self.cancel_conversion()
        super().closeEvent(event)
