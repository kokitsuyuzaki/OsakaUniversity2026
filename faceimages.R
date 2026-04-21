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
N_SUBJECTS <- 10       # 使用する被験者数
IMG_DIM <- c(92, 112)  # 顔画像の画素数 (幅 x 高さ)

# 関数定義
# 顔画像を2x5グリッドで表示する関数
plot_faces <- function(tensor, slice, colvec){
    layout(rbind(1:5, 6:10))
    for(i in seq(N_SUBJECTS)){
        image(tensor@data[,,i,slice], col=colvec, main=i)
    }
}

# 行列分解の係数（基底画像）を表示する関数
plot_coef <- function(coef_matrix, colvec){
    layout(rbind(1:5, 6:10))
    for(i in seq(N_SUBJECTS)){
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
    for(i in seq(N_SUBJECTS)){
        subplot(image(facedata@data[,,i,1], col=colvec,
            main=i, axes=FALSE),
            x=score_matrix[i,1], y=score_matrix[i,2],
            size=c(0.45, 0.45))
    }
}

# カラーパレットの設定
colvec <- rev(brewer.pal(9, "Greys"))

# ORL顔画像データのダウンロード（10人×40枚×10角度の4階テンソル）
options(timeout = 600)
facedata <- load_orl()
str(facedata)

# 顔画像の表示（角度1）
options(repr.plot.width=12, repr.plot.height=6)
plot_faces(facedata, slice=1, colvec)

# 顔画像の表示（角度2）
options(repr.plot.width=12, repr.plot.height=6)
plot_faces(facedata, slice=2, colvec)

# テンソルの3モード展開により行列化（10人分、角度1のみ使用）
matdata <- cs_unfold(facedata[,,1:N_SUBJECTS,1], m=3)@data

# 平均顔の計算と表示
options(repr.plot.width=6, repr.plot.height=6)
average_face <- rowMeans(matdata)
dim(average_face) <- IMG_DIM
image(average_face, col=colvec, main="Average Face")

# PCA（主成分分析）
res_pca <- prcomp(t(matdata), center=TRUE, scale=FALSE)

## PCAの係数（固有顔）の表示
options(repr.plot.width=12, repr.plot.height=6)
plot_coef(res_pca$rotation, colvec)

## PCAのスコアの表示
options(repr.plot.width=12, repr.plot.height=12)
plot_score(res_pca$x, facedata, colvec)

# NMF（非負値行列因子分解）
res_nmf <- NMF(matdata, J=N_SUBJECTS)

## NMFの係数（基底画像）の表示
options(repr.plot.width=12, repr.plot.height=6)
plot_coef(res_nmf$U, colvec)

## NMFのスコアの表示
options(repr.plot.width=12, repr.plot.height=12)
plot_score(res_nmf$V, facedata, colvec)

# セッション情報
sessionInfo()
