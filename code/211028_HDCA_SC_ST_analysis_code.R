# Spinal Cord ST Analysis


## Required packages
library(Matrix)
library(magrittr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(STutility)
library(Rcpp)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(clustree)
library(harmony)
library(niceRplots)


## Set WD
visium.dir <- "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/visium_data/"
export_path = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/r_results/"

## Load ST data
infoTable <- data.frame(section.name, samples, spotfiles, imgs, json, stringsAsFactors = FALSE)

se <- InputFromTable(infotable = infoTable, 
                     min.gene.count = 100, 
                     min.gene.spots = 5,
                     min.spot.count = 500,
                     platform="Visium")


## QC
se$percent.mt <- PercentageFeatureSet(se, pattern = "^MT")
se$percent.ribo <- PercentageFeatureSet(se, pattern = "RP")
VlnPlot(se, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0, cols = ("darkred"))


## Filtering
enids <- read.table(file = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/visium_data/GRCh38_genes.tsv", header = T, stringsAsFactors = T)
enids <- data.frame(apply(enids, 2, as.character), stringsAsFactors = F)
rownames(enids) <- enids$gene_id

se <- SubsetSTData(se, expression = nFeature_RNA > 200)

enids <- subset(enids, gene_biotype %in% c("protein_coding", "lncRNA"))
keep.genes <- intersect(rownames(se), enids$gene_name)
se <- se[keep.genes, ]

mt.genes <- grep(pattern = "^MT-", x = rownames(se), value = T)
keep.genes2 <- setdiff(rownames(se), mt.genes)
se <- se[keep.genes2, ]

hb.genes <- grep(pattern = "^HB", x = rownames(se), value = T)
keep.genes3 <- setdiff(rownames(se), hb.genes)
se <- se[keep.genes3, ]


malat1 <- grep(pattern = "MALAT1", x = rownames(se), value = T)
keep.genes4 <- setdiff(rownames(se), malat1)
se <- se[keep.genes4, ]


ribo.genes <- grep(pattern = "^RP", x = rownames(se), value = T)
keep.genes5 <- setdiff(rownames(se), ribo.genes)
se_filtered <- se[keep.genes5, ]


## Manual selection of sections
se_filtered <- LoadImages(se_filtered, time.resolve = FALSE)
se_filtered <- ManualAnnotation(se_filtered)
spot.remove <- colnames(se_filtered)[!se_filtered$labels %in% "Default"]
se_SC_sel <- SubsetSTData(se_filtered, spots = spot.remove)


## Crop Images
coords <- GetStaffli(se_SC_sel)@meta.data[, c("pixel_x", "pixel_y", "sample")]
coords$group <- se_SC_sel$labels

crop.windows <- Reduce(c, lapply(unique(coords$sample), function(s) {
  

  coords.sample <- subset(coords, sample %in% s)
  
  crop.geoms <- setNames(lapply(unique(coords.sample$group), function(grp) {
    coords.subset <- subset(coords.sample, group %in% grp)
    minxy <- apply(coords.subset[, c("pixel_x", "pixel_y")], 2, range)
    minxy[1, ] <- minxy[1, ] - 50 # Adjust padding along x axis
    minxy[2, ] <- minxy[2, ] + 50 # Adjust padding aong y axis
    minxy <- round(minxy)
    wh <- apply(minxy, 2, diff)
    geom <- paste0(wh[1], "x", wh[2], "+", minxy[1, 1], "+", minxy[1, 2])
  }), nm = rep(s, length(unique(coords.sample$group))))
}))


imgs <- GetStaffli(se_SC_sel)@rasterlists[["raw"]]
sf <- 400/2000

par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
for (i in seq_along(imgs)) {
  plot(as.cimg(imgs[[i]]), axes = F)
  for (j in which(names(crop.windows) == paste0(i))) {
    crpw <- crop.windows[[j]]
    hw <- as.numeric(unlist(strsplit(crpw, "x|\\+")))*sf
    rect(hw[3], hw[4], hw[3] + hw[1], hw[4] + hw[2])
  }
}


se.cropped.list <- lapply(seq_along(crop.windows), function(i) {
  se_cropped <- CropImages(se_SC_sel, crop.geometry.list = crop.windows[i])
  top_lbl <- names(table(se.cropped$labels) %>% sort(decreasing = TRUE))[1]
  se_cropped <- SubsetSTData(se_cropped, expression = labels %in% top_lbl)
})

se_cropped <- MergeSTData(se.cropped.list[[1]], se.cropped.list[2:length(se.cropped.list)])

msk.fun <- function(im) {
  im <- as.cimg(im[, , , 2])
  im <- imager::renorm(im, min = 0, max = 1)
  im <- 1 - im
  im <- im^0.2
  im <- imager::isoblur(im, 5)
  out <- im > otsu(im)
  out <- out %>%
    EBImage::as.Image() %>%
    EBImage::watershed()
  ids <- sort(table(out), decreasing = TRUE)
  id.remove <- names(ids)[-c(1, 2)]
  out <- rmObjects(out, index = id.remove)
  out <- EBImage::dilate(out, kern = makeBrush(size = 21, shape = "gaussian"))
  out <- as.cimg(out) %>% as.pixset()
  return(out)
}


se_cropped <- MaskImages(se_cropped, custom.msk.fkn = msk.fun, verbose = TRUE)


## Align and orient images
se_cropped <- ManualAlignImages(se_cropped, type = "masked", edges = FALSE, maxnum = 1e4)

## Save object
saveRDS(se_cropped, file = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/r_results/se_SC_cropped_aligned_76sections.rds")


## Normalization
se_norm <- SCTransform(se.cropped, vars.to.regress = "section.name")



## PCA
se_norm <- RunPCA(se_norm)


## Integration
se_harmony <- RunHarmony(se_norm, reduction = "pca", group.by.vars = "section.name", assay.use = "SCT")


## Colour palette 
pal <- c(RColorBrewer::brewer.pal(7,"Set1"), RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(8,"Pastel2") ,scales::hue_pal()(8))


## Clustering
se_harmony <- FindNeighbors(se_harmony, reduction = "harmony", dims = 1:30)
se_harmony <- FindClusters(se_harmony, verbose= F, resolution = 0.7)
se_harmony <- RunUMAP(se_harmony, reduction = "harmony", dims = 1:30, spread = 1, min.dist = 0.1, n.neighbors = 20)
DimPlot(se_harmony, dims = 1:2, label = T, cols = pal)
ST.FeaturePlot(se_harmony, features = "seurat_clusters", pt.size = 2, ncol = 4, cols = pal)

## Find marker genes
detable <- FindAllMarkers(se_harmony, only.pos = T,assay = "RNA", max.cells.per.ident = 1000,
                          min.pct = 0.05)
detable <- detable[ detable$p_val < 0.05,  ]
detable$pct.diff <- detable$pct.1 - detable$pct.2
detable$log.pct.diff <- log2(detable$pct.1 / (detable$pct.2+0.01) )

detable %>% group_by(cluster)  %>% top_n(-10, p_val) %>% top_n(5, log.pct.diff) -> top5

ord <- factor(sapply(unique(as.character(top5$gene)),function(x){getcluster(se_harmony, x, "seurat_clusters")}))



plot_dots(se_harmony, unique(as.character(top5$gene))[order(as.numeric( as.character(ord) ))], clustering = "seurat_clusters", show_grid = T,main = "top cluster markers",cex.main=1,font.main=1, pal = c("grey90","grey70", "darkblue"),cex.col = 1,srt = 90,cex.row = 1.1)

## Spatial plot - genes (one example)
ST.FeaturePlot(se_harmony, features = c("HIST1H4C", "HES6", "PTTG1", "S100B", "TOP2A", "MAD2L1"),
               indices = 7,
               ncol= 4, pt.size = 2,
               grid.ncol = 2,
               value.scale = 0:3)



