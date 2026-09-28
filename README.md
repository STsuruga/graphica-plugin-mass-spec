# graphica-plugin-mass-spec

[Graphica](https://github.com/STsuruga/Graphica) 用の質量分析 (MS) パック(P-316)。

## できること

- Bruker の `.d` フォルダを開く(ProteoWizard の msconvert で自動的に mzML に変換し、変換結果はキャッシュ)。mzML も直接開ける
- 「MS スペクトル」パネルで TIC を表示し、ドラッグで選んだ時間範囲の平均スペクトルを作る(背景の範囲の平均を差し引ける)
- ピークに m/z のラベル(桁数・本数・強度の表示を設定可)
- 組成式からモノアイソトピック質量・平均分子量・同位体パターン・付加イオンの m/z を計算し、実測に重ねて ppm 誤差と相対強度を比べる
- TIC・スペクトル・centroid・ピークのラベル・計算パターンを Graphica のプロットへ転送
- 棒(Stick)表示の種類

## 必要なもの

- Graphica 2.x
- `.d` を開く場合: [ProteoWizard](https://proteowizard.sourceforge.io/download.html) の Windows 64-bit 版(ベンダー形式を読める版)。
  **このプラグインには含まれていません。** 各自でインストールし、ダウンロード時に示されるベンダーのライセンスに同意してください。
  既定の場所(`%LOCALAPPDATA%\Apps\ProteoWizard …` または `Program Files`)に入れれば自動で見つかります。
  別の場所なら、初めて `.d` を開くときに `msconvert.exe` を選びます。

## 使い方

プラグインメニューから「MS スペクトル」パネルを表示します。操作の一覧はプラグインメニューの「MS パック: 使い方」でも見られます。

| グラフ | 操作 | 動作 |
|---|---|---|
| TIC | 左ドラッグ / Shift+左ドラッグ | 試料 / 背景の時間範囲 |
| TIC | 帯の端・帯の中をドラッグ | 範囲を伸縮・移動 |
| TIC | クリック | その時刻の1スキャン |
| スペクトル | 左ドラッグ | その m/z 範囲に拡大 |
| スペクトル | Ctrl+左ドラッグ | 矩形で拡大 |
| スペクトル | Shift+左ドラッグ | Δm/z と ppm を測る(同位体の間隔なら電荷数も) |
| スペクトル | ピークをクリック | ラベルを固定 / 解除 |
| 共通 | ホイール / Shift+ホイール | 横 / 縦の拡大縮小 |
| 共通 | 中ボタンドラッグ | パン |
| 共通 | ダブルクリック / Backspace | 全体 / 1つ前の表示 |
| 共通 | 右クリック | メニュー(転送など) |

照合は、パネルの「組成式から計算」か、実測スペクトルのデータセットを選んで プラグイン ▸ 解析 ▸ 同位体パターンと照合 で行います。
付加イオンは標準のセット(正: [M+H]+ [M+Na]+ [M+K]+ [M+NH4]+ [M]+•、負: [M−H]− [M+Cl]− [M]−• [M+HCOO]−)に、
`[2M+Na]+, [M+2H]2+, [M-H-H2O]-` のような表記を足せます。計算パターンの線幅は、半値全幅(m/z の単位、既定 0.1)で指定します。

本体のプロットに送ったピークのラベルは、本体の「点のラベル」機能で表示されます(点数が本体の上限、既定 1000 点を超えるデータセットには描かれないので、ラベル付きで送るのはピークだけです)。
棒を描かずラベルだけにしたいときは「ラベルだけ」で送るか、本体で Stick のデータセットの線の太さを 0 にします。

## msconvert の設定

プラグインは `msconvert <.d> --mzML --zlib` で変換します(profile のまま、m/z は 64 bit、zlib 圧縮)。自分で変換する場合も同じ設定にしてください。
**numpress 圧縮と peakPicking は使わないでください**(numpress は読めません。peakPicking をすると profile が失われます)。

Bruker micrOTOF(APCI、m/z 50〜3000)の測定での目安: 変換は1件 5〜6 秒、mzML は 26〜94 MB、
パネルで開いて全スキャンを平均するまで 0.3〜1.1 秒(20〜77 スキャン、1スキャン 150,912 点)。

## 開発環境

```
python -m venv .venv
.venv\Scripts\activate
pip install "graphica-plot>=2.0,<3"
pip install -r requirements-dev.txt
pytest
```

テストは合成したデータと合成した mzML だけを使います。

## zip のビルド

```
python scripts/build_zip.py --all      # dist/mass_spec-<version>.zip
```

Graphica の 編集 ▸ 環境設定 ▸「プラグイン」タブ ▸ プラグインをインストール... から入れます。

## ライセンス

MIT。同位体のデータは NIST(米国政府の著作物)。詳しくは [LICENSE](LICENSE)。
