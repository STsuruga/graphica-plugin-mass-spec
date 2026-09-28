"""MS ビューアのメニューから開く別ウィンドウ: 「同位体パターンの照合」「時間範囲の詳細設定」「表示設定」。

どちらもビューアの子ウィンドウで、ビューアが持つ窓口 ctx を通して本体とやりとりする(ビューアと一緒に閉じる)。
"""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from .analyzer import PRESET_NEGATIVE, PRESET_NONE, PRESET_POSITIVE

PRESETS = [PRESET_POSITIVE, PRESET_NEGATIVE, PRESET_NONE]
LABEL_MODES = [("top", "強い順に N 本"), ("percent", "相対強度 % 以上"), ("none", "表示しない")]
INTENSITY_MODES = [("none", "なし"), ("relative", "相対 %"), ("absolute", "絶対値")]


def spin(minimum, maximum, value, decimals=3, step=0.1, width=90):
    box = QDoubleSpinBox()
    box.setDecimals(decimals)
    box.setRange(minimum, maximum)
    box.setSingleStep(step)
    box.setValue(value)
    box.setMaximumWidth(width)
    box.setKeyboardTracking(False)
    return box


class CalcWindow(QWidget):
    """組成式・付加イオン・分解能を入れて計算し、照合表を出す。結果は3段目の枠・2段目・本体のプロットへ出せる。"""

    def __init__(self, viewer):
        super().__init__(viewer, Qt.WindowType.Window)
        self.viewer = viewer
        s = viewer.settings
        self.setWindowTitle("同位体パターンの照合 - MS パック")
        self.resize(720, 520)
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.formula_edit = QLineEdit(s["formula"])
        self.formula_edit.setPlaceholderText("例: C6H12O6")
        self.formula_edit.returnPressed.connect(viewer.run_calculation)
        form.addRow("組成式", self.formula_edit)
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(PRESETS)
        self.preset_combo.setCurrentIndex(min(max(s["preset"], 0), len(PRESETS) - 1))
        form.addRow("付加イオン", self.preset_combo)
        self.extra_edit = QLineEdit(s["extra_adducts"])
        self.extra_edit.setPlaceholderText("追加: [2M+Na]+, [M+2H]2+")
        self.extra_edit.returnPressed.connect(viewer.run_calculation)
        form.addRow("", self.extra_edit)
        row = QHBoxLayout()
        self.resolution_spin = spin(0, 1e7, s["resolution"], decimals=0, step=1000)
        self.resolution_spin.setSpecialValueText("実測から")
        self.tolerance_spin = spin(1, 5000, s["tolerance_ppm"], decimals=1, step=5, width=70)
        row.addWidget(self.resolution_spin)
        row.addWidget(QLabel("探す幅 ±ppm"))
        row.addWidget(self.tolerance_spin)
        row.addStretch(1)
        form.addRow("分解能 R", row)
        layout.addLayout(form)

        row = QHBoxLayout()
        self.calc_button = QPushButton("計算")
        self.calc_button.setDefault(True)
        self.calc_button.clicked.connect(viewer.run_calculation)
        self.overlay_check = QCheckBox("2段目に重ねる")
        self.overlay_check.setChecked(True)
        self.overlay_check.toggled.connect(lambda _on: viewer.show_overlays())
        row.addWidget(self.calc_button)
        row.addWidget(self.overlay_check)
        row.addStretch(1)
        layout.addLayout(row)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.summary)
        self.table = QTableWidget(0, 0)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        layout.addWidget(self.table, 1)

        row = QHBoxLayout()
        self.pane_combo = QComboBox()
        self.to_pane_button = QPushButton("選んだ付加イオンを枠に表示")
        self.to_pane_button.clicked.connect(self._send_selected_to_pane)
        self.found_to_panes_button = QPushButton("検出されたものを枠 1 から順に")
        self.found_to_panes_button.clicked.connect(viewer.found_patterns_to_panes)
        row.addWidget(QLabel("3段目の枠"))
        row.addWidget(self.pane_combo)
        row.addWidget(self.to_pane_button)
        row.addWidget(self.found_to_panes_button)
        row.addStretch(1)
        layout.addLayout(row)
        row = QHBoxLayout()
        self.transfer_button = QPushButton("計算パターンをプロットに転送")
        self.transfer_button.clicked.connect(viewer.transfer_calculation)
        self.copy_button = QPushButton("表をコピー")
        self.copy_button.clicked.connect(viewer.copy_table)
        row.addWidget(self.transfer_button)
        row.addWidget(self.copy_button)
        row.addStretch(1)
        layout.addLayout(row)
        self.update_pane_choices(viewer.pane_count())

    def update_pane_choices(self, count):
        current = self.pane_combo.currentIndex()
        self.pane_combo.clear()
        self.pane_combo.addItems([f"枠 {i + 1}" for i in range(count)])
        self.pane_combo.setCurrentIndex(min(max(current, 0), count - 1))

    def selected_adduct(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            return None
        item = self.table.item(rows[0], 0)
        return item.text() if item else None

    def _send_selected_to_pane(self):
        adduct = self.selected_adduct()
        if adduct is None:
            self.viewer.ctx.show_message("表で付加イオンの行を選んでください。", "MS パック")
            return
        self.viewer.pattern_to_pane(adduct, self.pane_combo.currentIndex())

    def show_table(self, header, rows):
        self.table.clear()
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                self.table.setItem(i, j, QTableWidgetItem(value))
        self.table.resizeColumnsToContents()


class RangeWindow(QWidget):
    """試料と背景の時間範囲を数値で直す。ふだんは TIC のドラッグで決め、ここは細かく合わせるときだけ使う。"""

    def __init__(self, viewer):
        super().__init__(viewer, Qt.WindowType.Window)
        self.viewer = viewer
        self.setWindowTitle("時間範囲の詳細設定 - MS パック")
        form = QFormLayout(self)
        self.sample_from = spin(0, 1e4, 0, width=110)
        self.sample_to = spin(0, 1e4, 0, width=110)
        self.bg_from = spin(0, 1e4, 0, width=110)
        self.bg_to = spin(0, 1e4, 0, width=110)
        for box in (self.sample_from, self.sample_to, self.bg_from, self.bg_to):
            box.setSuffix(" min")
            box.valueChanged.connect(viewer._on_range_spin_changed)
        for label, lo, hi in (("試料", self.sample_from, self.sample_to), ("背景", self.bg_from, self.bg_to)):
            row = QHBoxLayout()
            row.addWidget(lo)
            row.addWidget(QLabel("–"))
            row.addWidget(hi)
            row.addStretch(1)
            form.addRow(label, row)
        row = QHBoxLayout()
        full = QPushButton("試料を全範囲にする")
        full.clicked.connect(viewer.reset_sample_range)
        clear = QPushButton("背景を解除")
        clear.clicked.connect(viewer.clear_background)
        row.addWidget(full)
        row.addWidget(clear)
        row.addStretch(1)
        form.addRow(row)
        note = QLabel("背景の終わりが始まりより後でないときは、背景なしとして扱います。")
        note.setWordWrap(True)
        form.addRow(note)


class SettingsWindow(QWidget):
    """ラベル・桁数・背景・転送・msconvert の設定。変えるとすぐ反映して保存する。"""

    def __init__(self, viewer):
        super().__init__(viewer, Qt.WindowType.Window)
        self.viewer = viewer
        s = viewer.settings
        self.setWindowTitle("表示設定 - MS パック")
        layout = QVBoxLayout(self)

        group = QGroupBox("ピークのラベル")
        form = QFormLayout(group)
        self.label_mode_combo = QComboBox()
        for _key, text in LABEL_MODES:
            self.label_mode_combo.addItem(text)
        self.label_mode_combo.setCurrentIndex([k for k, _ in LABEL_MODES].index(s["label_mode"]))
        form.addRow("付け方", self.label_mode_combo)
        self.label_n_spin = QSpinBox()
        self.label_n_spin.setRange(1, 200)
        self.label_n_spin.setValue(s["label_top_n"])
        form.addRow("本数", self.label_n_spin)
        self.label_percent_spin = spin(0.0, 100.0, s["label_percent"], decimals=1, step=1.0)
        form.addRow("相対強度の下限 (%)", self.label_percent_spin)
        self.mz_decimals_spin = QSpinBox()
        self.mz_decimals_spin.setRange(0, 6)
        self.mz_decimals_spin.setValue(s["mz_decimals"])
        form.addRow("m/z の小数点以下の桁数", self.mz_decimals_spin)
        self.intensity_combo = QComboBox()
        for _key, text in INTENSITY_MODES:
            self.intensity_combo.addItem(text)
        self.intensity_combo.setCurrentIndex([k for k, _ in INTENSITY_MODES].index(s["intensity"]))
        form.addRow("強度の表示", self.intensity_combo)
        self.ppm_decimals_spin = QSpinBox()
        self.ppm_decimals_spin.setRange(0, 3)
        self.ppm_decimals_spin.setValue(s["ppm_decimals"])
        form.addRow("ppm の小数点以下の桁数", self.ppm_decimals_spin)
        layout.addWidget(group)

        group = QGroupBox("背景")
        form = QFormLayout(group)
        self.subtract_check = QCheckBox("背景の範囲の平均を引く")
        self.subtract_check.setChecked(s["subtract_background"])
        self.clip_check = QCheckBox("引いて負になった点は 0 にする")
        self.clip_check.setChecked(s["clip_negative"])
        form.addRow(self.subtract_check)
        form.addRow(self.clip_check)
        layout.addWidget(group)

        group = QGroupBox("プロットへの転送")
        form = QFormLayout(group)
        self.view_only_check = QCheckBox("スペクトルは表示中の m/z 範囲だけ送る")
        self.view_only_check.setChecked(s["transfer_view_only"])
        form.addRow(self.view_only_check)
        self.subplot_spin = QSpinBox()
        self.subplot_spin.setRange(1, 99)
        self.subplot_spin.setValue(s["subplot_target"] + 1)
        form.addRow("送り先のサブプロット", self.subplot_spin)
        self.centroid_min_spin = spin(0.0, 100.0, s["centroid_min_percent"], decimals=2, step=0.5)
        form.addRow("centroid に入れる下限 (%)", self.centroid_min_spin)
        layout.addWidget(group)


        group = QGroupBox("msconvert(ProteoWizard)")
        form = QFormLayout(group)
        self.msconvert_label = QLabel(s["msconvert_path"] or "自動で探す")
        self.msconvert_label.setWordWrap(True)
        self.msconvert_button = QPushButton("msconvert.exe を選ぶ…")
        self.msconvert_button.clicked.connect(self._choose_msconvert)
        self.msconvert_reset = QPushButton("自動に戻す")
        self.msconvert_reset.clicked.connect(self._reset_msconvert)
        row = QHBoxLayout()
        row.addWidget(self.msconvert_button)
        row.addWidget(self.msconvert_reset)
        row.addStretch(1)
        form.addRow(self.msconvert_label)
        form.addRow(row)
        layout.addWidget(group)
        layout.addStretch(1)

        for w in (self.label_mode_combo, self.intensity_combo):
            w.currentIndexChanged.connect(self._changed)
        for w in (self.label_n_spin, self.label_percent_spin, self.mz_decimals_spin, self.ppm_decimals_spin,
                  self.subplot_spin, self.centroid_min_spin):
            w.valueChanged.connect(self._changed)
        for w in (self.subtract_check, self.clip_check, self.view_only_check):
            w.toggled.connect(self._changed)
        self._update_visibility()

    def _update_visibility(self):
        mode = LABEL_MODES[self.label_mode_combo.currentIndex()][0]
        self.label_n_spin.setEnabled(mode == "top")
        self.label_percent_spin.setEnabled(mode == "percent")

    def _changed(self, *_):
        self._update_visibility()
        self.viewer.apply_settings(
            label_mode=LABEL_MODES[self.label_mode_combo.currentIndex()][0],
            label_top_n=self.label_n_spin.value(),
            label_percent=self.label_percent_spin.value(),
            mz_decimals=self.mz_decimals_spin.value(),
            intensity=INTENSITY_MODES[self.intensity_combo.currentIndex()][0],
            ppm_decimals=self.ppm_decimals_spin.value(),
            subtract_background=self.subtract_check.isChecked(),
            clip_negative=self.clip_check.isChecked(),
            transfer_view_only=self.view_only_check.isChecked(),
            subplot_target=self.subplot_spin.value() - 1,
            centroid_min_percent=self.centroid_min_spin.value(),
        )

    def _choose_msconvert(self):
        path = self.viewer.ask_msconvert(explain=False)
        if path:
            self.msconvert_label.setText(path)

    def _reset_msconvert(self):
        self.viewer.apply_settings(msconvert_path="")
        self.msconvert_label.setText("自動で探す")
