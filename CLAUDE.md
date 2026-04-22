# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

大阪大学の講義「データサイエンスA 行列テンソル分解・教師なし学習応用」のハンズオン資料。R言語で書かれた2つのハンズオンから構成される。Google Colaboratory上での実行を想定。

## Structure

各ハンズオンは `.R`（ソースコード）と `.ipynb`（Jupyter Notebook、R kernel）のペアで構成される。`.R` が正本で、`.ipynb` は Colab 配布用。

- **faceimages** (行列分解): ORL顔画像データセットに対するPCAとNMF
- **serology** (テンソル分解): COVID-19血清学データ（Zenodo）に対するTucker分解（HOSVD）

## Key R Packages

- `rTensor`: テンソルデータ構造とテンソル分解（Tucker分解、ORL顔画像データの読み込み）
- `nnTensor`: NMF（非負値行列因子分解）
- `TeachingDemos`: `subplot()` によるスコアプロット上への画像埋め込み（faceimages）
- `ggplot2` / `reshape2` / `RColorBrewer`: ヒートマップ可視化（serology）

## Notebook Execution

ipynbの実行・出力保存には R kernel (IRkernel) が必要。conda環境での実行例:

```bash
conda activate r-notebook
jupyter nbconvert --to notebook --execute --inplace faceimages.ipynb
jupyter nbconvert --to notebook --execute --inplace serology.ipynb
```

serology.ipynbはZenodoからデータをダウンロードするため、実行にはインターネット接続と `options(timeout = 600)` が必要。

## Conventions

- `.R` と `.ipynb` の内容は同期させること。`.R` を編集したら `.ipynb` も対応箇所を更新する。ipynbのコードセルは `.R` と完全一致させる。
- 各ファイルは「パッケージインストール → パッケージ読み込み → 定数設定 → 関数定義 → 分析 → sessionInfo()」の順で構成される。
- パッケージインストールは `install.packages()` で CRAN から source install する形式。
- プロット設定は `options(repr.plot.width=..., repr.plot.height=...)` で Jupyter 出力サイズを制御。
- コメントおよびMarkdownセルは日本語で記述する。
