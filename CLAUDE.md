# CLAUDE.md

Graphica(https://github.com/STsuruga/Graphica)用プラグイン P-316 質量分析 (MS) パックのリポジトリ。

## 最初に読むもの
- 開発ハブ(Artifact): https://claude.ai/artifact/GZ3LTLJjbxj1LQsAhZFg2o
  Artifact ツールの action: "read" で読む。共通ルール・API 早見表・この項目の仕様と状態がある。
- 実装計画(Artifact、非公開): https://claude.ai/artifact/ULochMNYRvtNCsPj28i75k
  ソースは docs/implementation_plan.html。直したらこのファイルを同じ URL へ再公開する。
- 本体の CLAUDE.md「Plugin API」節と docs/plugin_development.md(本体リポジトリ)

## 作業の範囲
- このリポジトリ専用。ほかのプラグインのリポジトリと Graphica 本体は変更しない
  (本体への書き込みは docs/dev/PLUGIN_DEVELOPMENT_PROGRESS.md の表への1行だけ)。
- 新しいプラグインを始めるときは、新しいチャットでハブの引継ぎプロンプトから始める。

## このプラグイン
- 項目: P-316。使うフック: register_importer, register_panel, register_analyzer, register_plot_type, register_menu_action
- v1.0 の範囲(2026-09-28 にユーザーと決めた内容。詳細は実装計画):
  - 読み込み: パネルで .d フォルダを選ぶと、利用者が入れた ProteoWizard の msconvert を QProcess で実行して
    mzML に変換し、ctx.data_dir にキャッシュして読む。mzML を直接開くこともできる。mzML は標準ライブラリと numpy で読む。
  - ドックパネル(タブごと)に開いた測定の一覧、TIC、平均スペクトル。TIC のドラッグ(+数値欄)で試料と背景の
    時間範囲を選び、背景を差し引ける。マウス操作(拡大・パン・Δm/z の測定・範囲の伸縮と移動)は計画の表のとおり。
  - ピークに m/z のラベル(桁数を設定可)。本体へは点のラベル(ラベル列)で出す。
  - TIC・平均スペクトル・centroid・計算パターンを本体のプロットへ転送(ctx.add_dataset)。
  - importer(.mzML): スキャンが1つならスペクトル、複数なら TIC。
  - analyzer「同位体パターンと照合」: 組成式・付加イオン(標準セットか表記の入力)・分解能(0 なら実測から推定)。
    実測に見つかった付加イオンだけ重ね、表(ppm 誤差・相対強度)は全部。パネルにも「組成式から計算」欄。
  - plot type「Stick」(棒表示。線種「なし」なら棒を描かず、ラベルだけにできる)。
- 見送ったこと・次の版に回したこと: .d の直接読み込み(baf2sql DLL)、m/z からの組成式の推定、フラグメント解析、
  LC-MS の netCDF(P-108 で扱う)。
- ライセンス: ProteoWizard と Bruker の DLL は同梱しない(利用者が各自インストールし、条件に同意する)。
  .baf を自前で解析しない。
- ユーザーの測定データ(未発表)はリポジトリに入れない。ファイル名もコミットしない。テストは合成した mzML だけ。
- isotopes.py は P-805(graphica-plugin-element-constants)の生成物のコピー。手で編集しない。

## コマンド
python -m venv .venv                  # 初回だけ。Python 3.11 以上
.venv\Scripts\activate               # macOS / Linux は source .venv/bin/activate
pip install "graphica-plot>=2.0,<3"   # 本体の未リリースの変更で試すときは pip install -e <PlotterApp>/Graphica_project
pip install -r requirements-dev.txt
pytest
graphica                               # 本体を起動(zip は 編集 ▸ 環境設定 ▸ プラグイン から入れる)
python scripts/build_zip.py --all      # dist/mass_spec-<version>.zip

## ルール
- Graphica 本体のコードは変更しない。足りない拡張点は本体の Issue(ユーザーの了承を得て作成)とハブの note に記録する。
- 依存は本体同梱のパッケージのみ(PySide6 6.11 / matplotlib 3.11 / numpy / pandas / scipy / openpyxl / xlrd)。
- matplotlib の色は組で返ることがあるので、Qt に渡す前に matplotlib.colors.to_hex で #rrggbb にする。
- プラグイン内は相対 import。他プラグインは import できない。
- 受け取った Dataset は書き換えない。新しい Dataset は name / df / x_col_name / y_col_name 必須。
- 本体から import するのは graphica.plugin / graphica.plugin.testing だけ。本体の操作は窓口 ctx(PluginContext)で行う。
- ctx の呼び出しは GUI スレッドだけ(msconvert の完了通知は QProcess のシグナルで GUI スレッドに戻る)。
- リリースしたら plugin.json の version とタグを揃え、ハブの db(collection "plugins", doc_id "P-316")を更新する。

## 現状
- 2026-09-28: 仕様確定。リポジトリと仮想環境を作成。chemistry.py(組成式・付加イオン・同位体パターン・ガウス)と
  テストまで完了。次は ProteoWizard の導入を待って実データ(mzML の中身・点数・時間)を確認し、mzml.py / convert.py。
