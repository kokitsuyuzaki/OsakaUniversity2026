# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

大阪大学の講義「データサイエンスA 行列テンソル分解・教師なし学習応用」のハンズオン資料。R言語で書かれた2つのハンズオンから構成される。Google Colaboratory上での実行を想定。

## Structure

各ハンズオンは `.R`（ソースコード）と `.ipynb`（Jupyter Notebook、R kernel）のペアで構成される。`.R` が正本で、`.ipynb` は Colab 配布用。

- **faceimages**: ORL顔画像データセットに対する行列分解（PCA, SVD, NMF, ICA）
- **serology**: COVID-19血清学データ（Zenodo）に対するTucker分解（HOSVD）

## Key R Packages

- `rTensor`: テンソルデータ構造とテンソル分解（Tucker, CP等）
- `nnTensor`: NMF（非負値行列因子分解）
- `iTensor`: ICA（独立成分分析）
- `TeachingDemos`: `subplot()` によるスコアプロット上への画像埋め込み（faceimages）
- `ggplot2` / `reshape2`: ヒートマップ可視化（serology）

## Conventions

- `.R` と `.ipynb` の内容は同期させること。`.R` を編集したら `.ipynb` も対応箇所を更新する。
- 各ファイルは「Package Download → Package Loading → (Function Definition) → 分析 → sessionInfo()」の順で構成される。
- パッケージインストールは `install.packages()` で CRAN から source install する形式。
- プロット設定は `options(repr.plot.width=..., repr.plot.height=...)` で Jupyter 出力サイズを制御。
