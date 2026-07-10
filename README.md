# Counting Galaxies in a Deep Sky Image

*[日本語はこちら / Japanese version below](#日本語版)*

<p align="center">
  <img src="figures/tiles.png" alt="The sky image with the analysed regions highlighted" width="600"/>
</p>

## What is this?

This project takes a real photograph of deep space (taken with the 4-metre telescope at Kitt Peak Observatory, USA) and answers a simple-sounding question: **how many galaxies are in this picture, and how bright is each one?**

That's harder than it sounds. The image contains thousands of faint smudges of light sitting on top of background noise, plus bright stars that "bleed" streaks across the camera sensor. So I wrote a program that:

1. **Measures the background** — fits a curve to the sky brightness so we know what "empty sky" looks like
2. **Finds the galaxies** — hunts for spots that are clearly brighter than empty sky, and works out how big each one is by expanding a ring around it until the brightness drops back to background level
3. **Measures each galaxy's brightness** — adds up the light inside each galaxy, and subtracts the background light that would be there anyway
4. **Counts them by brightness** — builds a catalogue of every galaxy found, then plots how the number of galaxies grows as you count fainter and fainter ones

The final result is compared against a classic theoretical prediction (that the count should grow at a specific rate as you go fainter). My measured rate came out noticeably shallower than the theory — mostly because the faintest galaxies get harder and harder to detect, so the count flattens out. Quantifying that difference, and why it happens, was the point of the project.

<p align="center">
  <img src="figures/number_counts.png" alt="Galaxy counts vs brightness, with the fitted line and theoretical prediction" width="600"/>
</p>

## The files

| File | What it does |
|---|---|
| `src.py` | The main pipeline. Reads the image, detects galaxies in 8 hand-picked regions (chosen to avoid the damaged/noisy parts of the image), measures their brightness, saves a catalogue, and makes the plots. |
| `threshold_sweep.py` | A testing tool. The detector needs a "how bright counts as a galaxy?" setting — this script tries a range of settings on a crowded patch of the image and measures when neighbouring galaxies start getting wrongly merged into one. It's how I justified the setting used in `src.py`. |
| `visualise_tiles.py` | Draws the 8 analysed regions on top of the image (the picture at the top of this page), so you can see exactly which parts were used and which were avoided. |
| `catalogue.csv` | Example output from `src.py` — the list of detected galaxies with their positions and brightnesses. |

## Figures

| | |
|---|---|
| <img src="figures/gaussian_fit.png" width="380"/> | **Background fit** — a histogram of pixel brightness. The big peak is empty sky; the small tail on the right is the galaxies. |
| <img src="figures/detections_overlay.png" width="380"/> | **Detections** — every galaxy the program found, marked on the original image. |
| <img src="figures/number_counts.png" width="380"/> | **The result** — how the galaxy count grows with faintness, with my fitted line (cyan) vs. the theoretical prediction (red). |

## Running it

You'll need Python 3 and a few libraries:

```
pip install -r requirements.txt
```

You'll also need the image itself (`mosaic.fits`, ~23 MB) placed in the same folder — it's the deep survey image provided as part of Imperial College London's third-year lab course, so it isn't included in this repository.

Then:

```
python src.py              # full pipeline: detect, measure, catalogue, plots
python visualise_tiles.py  # show which regions of the image were analysed
python threshold_sweep.py --fits mosaic.fits --x1 390 --x2 1300 --y1 1050 --y2 1700 \
    --seed_sigma 3.0 --ring_sigmas 3.0,2.5,2.0,1.5,1.0 --output_dir sweep_results
```

## Report

The full write-up (method, results, error analysis) is in [`report.pdf`](report.pdf).

---

# 日本語版

<a name="日本語版"></a>

## これは何？

このプロジェクトは、実際の深宇宙の画像(アメリカ・キットピーク天文台の4m望遠鏡で撮影)を使って、シンプルに聞こえる問いに答えるものです。**この画像には銀河がいくつ写っていて、それぞれどのくらいの明るさなのか？**

これは見た目より難しい問題です。画像には、背景ノイズの上に何千もの淡い光の点が写っており、さらに明るい星がセンサー上に筋状の「にじみ」を作っています。そこで、次の処理を行うプログラムを書きました。

1. **背景を測定する** — 空の明るさに曲線をフィットさせ、「何もない空」がどう見えるかを把握する
2. **銀河を見つける** — 空より明らかに明るい点を探し、その周りにリングを広げていき、明るさが背景レベルに戻るまでの範囲からサイズを求める
3. **各銀河の明るさを測る** — 銀河の範囲内の光を合計し、もともとそこにあるはずの背景光を差し引く
4. **明るさごとに数える** — 検出した銀河のカタログを作成し、暗い銀河まで含めていくと数がどう増えるかをプロットする

最終結果は、古典的な理論予測(暗くなるにつれて数が特定の割合で増えるという予測)と比較しました。測定された増加率は理論よりも明らかに緩やかでした。これは主に、暗い銀河ほど検出が難しくなり、数え上げが頭打ちになるためです。その差を定量化し、原因を考察することがこのプロジェクトの目的でした。

## ファイル構成

| ファイル | 内容 |
|---|---|
| `src.py` | メインのパイプライン。画像を読み込み、(画像の損傷部分やノイズの多い部分を避けて選んだ)8つの領域で銀河を検出し、明るさを測定し、カタログを保存してプロットを作成します。 |
| `threshold_sweep.py` | 検証用ツール。検出器には「どのくらい明るければ銀河とみなすか」という設定が必要です。このスクリプトは、銀河が密集した領域で設定値を変えながら、隣り合う銀河が誤って1つに統合され始めるポイントを測定します。`src.py` で使用した設定値の根拠となっています。 |
| `visualise_tiles.py` | 解析に使用した8つの領域を画像上に描画します(ページ上部の画像)。どの部分を使い、どの部分を避けたかが一目で分かります。 |
| `catalogue.csv` | `src.py` の出力例 — 検出された銀河の位置と明るさの一覧です。 |

## 実行方法

Python 3 といくつかのライブラリが必要です:

```
pip install -r requirements.txt
```

また、画像ファイル(`mosaic.fits`、約23MB)を同じフォルダに置く必要があります。これはインペリアル・カレッジ・ロンドンの3年次実験で提供されたデータのため、このリポジトリには含まれていません。

```
python src.py              # フルパイプライン: 検出、測定、カタログ作成、プロット
python visualise_tiles.py  # 解析に使用した領域の表示
```

## レポート

詳細な手法・結果・誤差解析は [`report.pdf`](report.pdf) をご覧ください。