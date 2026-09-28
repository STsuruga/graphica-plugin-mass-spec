# graphica-plugin-mass-spec

[Graphica](https://github.com/STsuruga/Graphica) 用の質量分析 (MS) パック(P-316)。開発中です。

## できること(v1.0 の予定)

- Bruker の `.d` フォルダ(ProteoWizard の msconvert で自動的に mzML に変換)と mzML の読み込み
- TIC の表示と、ドラッグで選んだ時間範囲の平均スペクトル(背景の差し引き付き)
- ピークへの m/z のラベル(桁数を設定可)、Graphica のプロットへの転送
- 組成式からモノアイソトピック質量・平均分子量・同位体パターン・付加イオンの m/z を計算し、実測に重ねて ppm 誤差と相対強度を比べる
- 棒(Stick)表示

## 必要なもの

- Graphica 2.x
- `.d` を開く場合: [ProteoWizard](https://proteowizard.sourceforge.io/download.html)(Windows 64-bit、ベンダー形式を読める版)。
  このプラグインには含まれていません。各自でインストールし、ダウンロード時に示されるベンダーのライセンスに同意してください。

## 開発環境

```
python -m venv .venv
.venv\Scripts\activate
pip install "graphica-plot>=2.0,<3"
pip install -r requirements-dev.txt
pytest
```

## zip のビルド

```
python scripts/build_zip.py --all      # dist/mass_spec-<version>.zip
```

Graphica の 編集 ▸ 環境設定 ▸「プラグイン」タブ ▸ プラグインをインストール... から入れます。

## ライセンス

MIT。同位体のデータは NIST(米国政府の著作物)。詳しくは [LICENSE](LICENSE)。
