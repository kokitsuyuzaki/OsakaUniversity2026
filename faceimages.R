# パッケージのインストール
list.of.packages <- c("rTensor", "nnTensor",
    "RColorBrewer", "TeachingDemos")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages,
        repos="https://cloud.r-project.org/", type="source")
}

# パッケージの読み込み
library("rTensor")
library("nnTensor")
library("RColorBrewer")
library("TeachingDemos")

# 定数の設定
N_DISPLAY <- 10        # 表示する被験者数
N_RANK <- 20           # 行列分解のランク（基底数）
IMG_DIM <- c(92, 112)  # 顔画像の画素数 (幅 x 高さ)

# 関数定義
# 顔画像を2x5グリッドで表示する関数
plot_faces <- function(tensor, slice, colvec){
    layout(rbind(1:5, 6:10))
    for(i in seq(N_DISPLAY)){
        image(tensor@data[,,i,slice], col=colvec, main=i)
    }
}

# 行列分解の係数（基底画像）を表示する関数
plot_coef <- function(coef_matrix, colvec, n=N_RANK){
    nrow <- ceiling(n / 5)
    layout(matrix(seq(nrow * 5), nrow=nrow, byrow=TRUE))
    for(i in seq(n)){
        tmp <- as.matrix(coef_matrix[,i])
        dim(tmp) <- IMG_DIM
        image(tmp, col=rev(colvec), main=i)
    }
}

# スコア（各被験者の重み）を2次元散布図で表示する関数
plot_score <- function(score_matrix, facedata, colvec){
    plot(score_matrix[,1:2], pch=16, cex=2,
        axes=FALSE, xlab="", ylab="")
    abline(v=0)
    abline(h=0)
    n <- dim(facedata@data)[3]
    for(i in seq(n)){
        subplot(image(facedata@data[,,i,1], col=colvec,
            main=i, axes=FALSE),
            x=score_matrix[i,1], y=score_matrix[i,2],
            size=c(0.45, 0.45))
    }
}

# カラーパレットの設定
colvec <- rev(brewer.pal(9, "Greys"))

# ORL顔画像データのダウンロード
# AT&T Laboratoriesが収集した顔画像データセット
# 40人の被験者 × 10種類の撮影条件（角度・表情・照明）
# 各画像は92×112ピクセルのグレースケール → 4階テンソル (92×112×40×10)
options(timeout = 600)
facedata <- load_orl()
str(facedata)

# 顔画像の表示（角度1、最初の10人）
# 同一の撮影条件でも、被験者ごとに顔の形状・特徴が異なることを確認
options(repr.plot.width=12, repr.plot.height=6)
plot_faces(facedata, slice=1, colvec)

# 顔画像の表示（角度2、最初の10人）
# 撮影条件が変わると、同じ被験者でも見え方が変化する
options(repr.plot.width=12, repr.plot.height=6)
plot_faces(facedata, slice=2, colvec)

# テンソルの3モード展開により行列化（全40人、角度1のみ使用）
# 92×112=10304次元の画素ベクトルが40本並んだ行列 (10304×40) を作成
# この高次元データを行列分解により低次元に要約する
matdata <- cs_unfold(facedata[,,,1], m=3)@data

# 平均顔の計算と表示
# 40人の画素値を平均した「平均顔」
# PCAではこの平均顔を引いた残差（個人差）を分析する
options(repr.plot.width=6, repr.plot.height=6)
average_face <- rowMeans(matdata)
dim(average_face) <- IMG_DIM
image(average_face, col=colvec, main="Average Face")

# =============================================================================
# PCA（主成分分析）
# データを中心化（平均を引く）した上で、分散が最大となる方向を順に求める
# 各主成分は顔全体にわたるグローバルなパターン（固有顔, Eigenface）になる
# =============================================================================
res_pca <- prcomp(t(matdata), center=TRUE, scale=FALSE)

## PCAの係数（固有顔, Eigenface）の表示（上位N_RANK個）
## 正負の値を取り、顔全体に広がるパターンになる
## 第1主成分は全体的な明暗、第2主成分以降は顔の左右差や上下差などを捉える
## 表示は黒=値が大きい、白=値が小さい（元画像と白黒反転）
options(repr.plot.width=12, repr.plot.height=12)
plot_coef(res_pca$rotation, colvec)

## PCAのスコアの表示（全40人）
## 各被験者が主成分空間上のどこに位置するかを示す
## 近い位置にある被験者は、顔の特徴が似ていることを意味する
options(repr.plot.width=12, repr.plot.height=12)
plot_score(res_pca$x, facedata, colvec)

# =============================================================================
# NMF（非負値行列因子分解）
# 元の行列を2つの非負行列の積 X ≈ U V^T に分解する
# 非負制約により、基底画像は「足し合わせて顔を構成するパーツ」になる
# PCAと異なり、引き算を使えないため、局所的な特徴を捉えやすい
# =============================================================================
res_nmf <- NMF(matdata, J=N_RANK)

## NMFの係数（基底画像, パーツ）の表示（N_RANK個）
## PCAの固有顔と比較すると、顔の局所的な領域（額・目・鼻・口・輪郭など）に
## 対応するパーツベースの表現が得られる
## 表示は黒=値が大きい（そのパーツが強い）、白=値が小さい
options(repr.plot.width=12, repr.plot.height=12)
plot_coef(res_nmf$U, colvec)

## NMFのスコアの表示（全40人）
## 各被験者がどのパーツをどの程度持っているかの重み
## スコアも非負のため、各パーツの「存在量」として解釈できる
options(repr.plot.width=12, repr.plot.height=12)
plot_score(res_nmf$V, facedata, colvec)

# セッション情報
sessionInfo()
