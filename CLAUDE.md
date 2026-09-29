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
- v1.0 の範囲(2026-09-28〜29 にユーザーと決めた内容。詳細は実装計画):
  - 読み込み: .d(analysis.baf)を baf.py が直接読む(Win/Mac、Bruker のライブラリも ProteoWizard も不要)。
    圧縮方式が違う測定だけ、ProteoWizard の msconvert があれば QProcess で mzML に変換して ctx.data_dir にキャッシュ(convert.py)。
    mzML は mzml.py(標準ライブラリと numpy)。importer は .mzML と .baf(スキャンが1つならスペクトル、複数なら TIC)。
  - MS ビューアは独立したウィンドウだけ(ドックでは使わない、ユーザー判断)。プラグイン ▸ MS ビューアを開く で開き、
    パネル(register_panel)はタブごとの ctx を持つためだけに登録し、表示されたらドックを隠してウィンドウを出す。
    ウィンドウは WA_QuitOnClose を切る(開いたままでも本体を閉じたら終わる)。macOS ではメニューバーを setNativeMenuBar(False)。
  - 画面: メニューバー(ファイル・表示・解析・転送・ヘルプ)、左に測定の一覧(チェックで TIC、クリックで2段目の対象)、
    右に3段(TIC / 選んだ範囲の平均スペクトル / 比較の枠。枠は1つに1件、数と横軸の同期は3段目の下)。
    開いた直後は範囲なし。時間範囲の数値・表示設定・照合は別ウィンドウ。状態の行は出さない(変換中だけ右下に表示)。
  - ラベル: 表示範囲でいちばん強いピークに対する % 以上(既定 20%)。同位体の間隔で並ぶピークは1つのシグナルとし、
    代表(いちばん強い)だけ必ず出す。ほかはラベル・線と重ならず、枠からはみ出さないときだけ。本体へは点のラベル(代表だけ)。
  - 照合: 組成式・付加イオン(標準セットか表記。自動で極性に合わせない、ユーザー判断)・計算パターンの半値全幅
    (m/z、既定 0.1。実測からの推定はしない、ユーザー判断)。結果は表、3段目の空いている枠(なければ1つ増やす)、本体へ転送。
  - plot type「Stick」(線の太さが 0 なら棒を描かず、ラベルだけ)。
- 見送ったこと・次の版に回したこと: m/z からの組成式の推定、フラグメント解析、LC-MS の netCDF(P-108 で扱う)。
- ライセンス: ProteoWizard と Bruker のライブラリは同梱しない。.baf はデータファイルだけを読み、Bruker のプログラムは
  使わない・改変しない・逆コンパイルしない。Bruker の使用許諾(TopSpin の例)はソフトウェアのリバースエンジニアリングを禁じ、
  データと結果は顧客のものとしている。これを読んだうえで、ユーザーの判断で .baf の読み込みを公開した(2026-09-29)。
- ユーザーの測定データ(未発表)はリポジトリに入れない。ファイル名もコミットしない。テストは合成した mzML と BAF だけ。
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

## 実データで分かったこと(ユーザーの micrOTOF の測定3件、2026-09-28)
- msconvert(ProteoWizard 3.0.26267、`%LOCALAPPDATA%\Apps\ProteoWizard … 64-bit\`)の変換は1件 5〜6 秒、mzML は 26〜94 MB。
- profile だけ。20〜77 スキャン(約 1 秒ごと)、1スキャン常に 150,912 点(0 も省略されない。0 が 47〜99.7%)。
  m/z は 64 bit・強度は 32 bit・zlib。点の間隔は約 20 ppm。
- スキャン間の m/z の格子のずれは最大 0.36 ppm → 点の番号ごとに平均してよい(spectra.GRID_TOLERANCE_PPM)。
- 分解能は負イオンで 1.1〜1.7 万、正イオンの弱いピークで約 6 千。半値全幅に 2〜6 点しか乗らないので頂点はガウスの3点補間。
- 装置の TIC(mzML の total ion current)は profile の総和と値が違う。表示は装置の値。
- パネルで開いて全スキャンを平均し、ピークを拾うまで 0.3〜1.1 秒(GUI スレッド)。配列はそのつど復号(全部持つと 186 MB)。
- テストの合成データ(tests/synthetic.py)はこの形に合わせてある。

## 現状
- 2026-09-29: v1.0.0 をリリース(https://github.com/STsuruga/graphica-plugin-mass-spec/releases/tag/v1.0.0)。
  Windows(仮想環境の graphica)と Mac の実機でユーザーが動作を確認。テスト 146 件、CI(windows-latest)通過。
- 次: 未定。候補は m/z からの組成式の推定、フラグメント解析、ほかの Bruker 機種の .baf(圧縮方式が違えば msconvert に任せている)。
