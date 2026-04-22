# パッケージのインストール
list.of.packages <- c("rTensor", "ggplot2", "reshape2", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages,
        repos="https://cloud.r-project.org/", type="source")
}

# パッケージの読み込み
library("rTensor")
library("ggplot2")
library("reshape2")
library("RColorBrewer")

# 定数の設定
# 元論文ではCP分解（ランク2）を使用し、2成分で分散の74%を説明した
# ここではTucker分解（HOSVD）を用いるが、同様に2成分で主要パターンを抽出する
J <- 2  # Tucker分解のランク（成分数）
HEATMAP_COLORS <- brewer.pal(11, "PiYG")
HEATMAP_THEME <- theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 30))

# 関数定義
# Zenodoから血清学データをダウンロードし、3階テンソル（受容体×抗原×検体）に整形する関数
# データ概要:
#   COVID-19患者の血清学的アッセイ（蛍光強度）を測定したデータ
#   - 受容体/抗体 (11種): IgG1, IgG2, IgG3, IgA1, IgA2, IgM, FcRalpha, FcR2A, FcR2B, FcR3A, FcR3B
#   - 抗原 (6種): S, RBD, N, S1, S2, S1.Trimer（SARS-CoV-2由来）
#   - 検体 (438サンプル): 一人の患者から複数時点で採取されるため、サンプル数 > 患者数
#   患者の重症度分類: Negative(陰性), Mild(軽症), Moderate(中等症), Severe(重症), Deceased(死亡)
load_serology <- function(){
    # Download from Zenodo Server
    td <- tempdir()
    tempfile <- paste0(td, "/meyer-lab-systemsSerology-fd9ef61.zip")
    download.file("https://zenodo.org/record/5184449/files/meyer-lab/systemsSerology-v1.0.zip?download=1", tempfile)
    # Preprocessing
    unzip(tempfile, exdir=td)
    csvfile <- paste0(td, "/meyer-lab-systemsSerology-fd9ef61/syserol/data/ZoharCovData.csv")
    data <- read.csv(csvfile, row.names=1, header=TRUE)
    serology <- t(t(as.matrix(data[, 23:140])) - unlist(data[443, 23:140]))
    serology <- serology[1:438, ]
    # Data frame -> Array
    receptor <- c("IgG1", "IgG2", "IgG3", "IgA1", "IgA2", "IgM", "FcRalpha", "FcR2A", "FcR2B", "FcR3A", "FcR3B")
    antigen <- c("S", "RBD", "N", "S1", "S2", "S1.Trimer")
    arr <- array(0, dim=c(11, 6, 438))
    dimnames(arr) <- list(
        receptor = receptor,
        antigen = antigen,
        samples = rownames(serology))
    for(i in receptor){
        for(j in antigen){
            arr[i, j, ] <- serology[, paste0(i, "_", j)]
        }
    }
    # Log Transformation
    arr[which(arr < 0)] <- 10
    arr <- log10(arr)
    # Centering
    for(k in seq_len(dim(arr)[3])){
        arr[,,k] <- arr[,,k] - mean(arr[,,k])
    }
    # Array -> Tensor
    covid19 <- as.tensor(arr)
    # Group Label
    group <- cbind(rownames(data)[1:438], data$group[1:438])
    colnames(group) <- c("Sample", "Group")
    # Output
    list(covid19=covid19, group=group)
}

# Tucker分解の因子行列をヒートマップ用のデータフレームに変換する関数
prepare_pattern <- function(U, dim_names, col_names, group=NULL){
    rownames(U) <- dim_names
    colnames(U) <- paste("Component", seq(ncol(U)))
    df <- melt(U)
    if(!is.null(group)){
        df <- merge(df, group, by.x="Var1", by.y="Sample")
        colnames(df) <- c(col_names[1], "Component", "Value", "Group")
        df$Group <- factor(df$Group,
            level=c("Negative", "Mild", "Moderate", "Severe", "Deceased"))
    } else {
        colnames(df) <- c(col_names[1], "Component", "Value")
    }
    df[[col_names[1]]] <- factor(df[[col_names[1]]],
        level=rev(unique(df[[col_names[1]]])))
    df
}

# ヒートマップを描画する関数
plot_heatmap <- function(df, y_col, facet_col=NULL, extra_theme=NULL){
    g <- ggplot(df, aes(x=Component, y=.data[[y_col]], fill=Value)) +
        geom_tile() +
        scale_fill_gradientn("value", colours=HEATMAP_COLORS) +
        HEATMAP_THEME
    if(!is.null(facet_col)){
        g <- g + facet_wrap(as.formula(paste0("~", facet_col)), ncol=2)
    }
    if(!is.null(extra_theme)){
        g <- g + extra_theme
    }
    g
}

# プロットサイズの設定
options(repr.plot.width=6, repr.plot.height=12)

# COVID-19血清学データのダウンロードと前処理
# 3階テンソル (11受容体 × 6抗原 × 438検体) を構築する
# 前処理: 陰性コントロールの減算 → 負値の補正 → log10変換 → 検体ごとの中心化
options(timeout = 600)
serology_data <- load_serology()
covid19 <- serology_data$covid19
group <- serology_data$group

# =============================================================================
# Tucker分解（HOSVD: 高次特異値分解）
# テンソルを各モードの因子行列とコアテンソルに分解する
# 元論文ではCP分解を使用しているため結果は完全には一致しないが、
# 各モード（受容体・抗原・検体）の主要パターンを抽出できる点は共通している
# =============================================================================
res_tucker <- hosvd(covid19, ranks=rep(J, length=3))

# 受容体パターンのヒートマップ（モード1の因子行列）
# 各受容体/抗体が2つの成分にどの程度寄与しているかを示す
# 元論文の結論:
#   Component 1（急性型）: IgG3, IgM, IgA や各種Fcタンパク質が強く関与
#   Component 2（長期型）: IgG1に特異的な応答パターン
# HOSVDでも類似の傾向が見られるか確認する
df_r <- prepare_pattern(res_tucker$U[[1]],
    dimnames(covid19@data)$receptor, c("Receptor"))
plot_heatmap(df_r, "Receptor")

# 抗原パターンのヒートマップ（モード2の因子行列）
# 各抗原（SARS-CoV-2のスパイクタンパク質やその部位）が
# 2つの成分にどの程度寄与しているかを示す
df_a <- prepare_pattern(res_tucker$U[[2]],
    dimnames(covid19@data)$antigen, c("Antigen"))
plot_heatmap(df_a, "Antigen")

# 検体パターンのヒートマップ（モード3の因子行列、重症度グループ別）
# 各検体の2成分上のスコアを、重症度グループごとに並べて表示する
# 重症度分類: Negative(陰性), Mild(軽症), Moderate(中等症), Severe(重症), Deceased(死亡)
# グループ間でスコアの分布に違いがあれば、免疫応答パターンが重症度と関連することを示唆する
df_s <- prepare_pattern(res_tucker$U[[3]],
    dimnames(covid19@data)$samples, c("Sample"), group=group)
plot_heatmap(df_s, "Sample", facet_col="Group",
    extra_theme=theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_blank()))

# セッション情報
sessionInfo()
