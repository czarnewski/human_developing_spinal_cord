# Spinal Cord Stereoscope Analysis

## Set WD
stsc.dir <- '/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/stsc/210524_stsc/'


##Import, format and add stsc output as an assay to Seurat object
stsc.files <- list.files(stsc.dir, recursive = TRUE, full.names = TRUE, pattern = '^W')
stsc.output <- lapply(stsc.files, read.delim, row.names = 1)
stsc.output <- lapply(stsc.output, t)

ids <- sapply(stsc.output, function(x) {
  unlist(strsplit(colnames(x)[1], "_"))[2]
})

stsc.files <- setNames(stsc.files, nm = ids)
stsc.files <- stsc.files[paste0(1:24)]

stsc.output2 <- do.call(cbind, stsc.output)
emptyMat <- matrix(0, nrow = nrow(stsc.output2), ncol = ncol(se_harmony))
colnames(emptyMat) <- colnames(se_harmony)
rownames(emptyMat) <- rownames(stsc.output2)
emptyMat[, intersect(colnames(se_harmony), colnames(stsc.output2))] <- stsc.output2[, intersect(colnames(se_harmony), colnames(stsc.output2))]


# Create assay object
stereoscope <- CreateAssayObject(data = emptyMat)
se_harmony[["stereoscope"]] <- stereoscope
DefaultAssay(se_harmony) <- "stereoscope"


## Spatial plot - cell types (one example)
ST.FeaturePlot(se_harmony, features = "ENs", ncol = 6)

## Add PCW
se_harmony$age <- "unknown"

se_harmony$age[grep("V19N18-116", se_harmony$section.name)] <- "W12"
se_harmony$age[grep("V10T03-302", se_harmony$section.name)] <- "W5"
se_harmony$age[grep("V10T03-290_A1", se_harmony$section.name)] <- "W5"
se_harmony$age[grep("V10T03-290_B1", se_harmony$section.name)] <- "W5"
se_harmony$age[grep("V10T03-290_C1", se_harmony$section.name)] <- "W12"
se_harmony$age[grep("V10T03-290_D1", se_harmony$section.name)] <- "W12"
se_harmony$age[grep("V10T03-288", se_harmony$section.name)] <- "W9"
se_harmony$age[grep("V10F24-107", se_harmony$section.name)] <- "W8"
se_harmony$age[grep("V10F24-106", se_harmony$section.name)] <- "W12"


## Subset per PCW
W5_stsc <- SubsetSTData(se_harmony, age == "W5")
DefaultAssay(W5_stsc) <- "stereoscope"

## Remove Blood cells
features.remove <- rownames(W5_stsc)[!rownames(W5_stsc) %in% "Blood"]
w5_stsc_noblood <- SubsetSTData(W5_stsc, features = features.remove)


## Cord plot per PCW
df <- t(GetAssayData(w5_stsc_noblood[["stereoscope"]], slot = "data"))
df <- as.data.frame(t(GetAssayData(w5_stsc_noblood[["stereoscope"]], slot = "data")))


cts <- colnames(df)
nspt <- nrow(df)
ntcs <- length(cts)
vals <- as.data.frame(matrix(0,ntcs^2 / 2 - ntcs / 2 -ntcs,3))
colnames(vals) <- c("from","to","value")



k = 0

for (ii in cts){
  for (jj in cts) {
    if (ii <jj){
      k = k +1
      prod <-df[[ii]] %*% df[[jj]] / nspt
      vals[k,1] <- ii
      vals[k,2] <- jj
      vals[k,3] <- prod
    }
  }
}

chordDiagram(vals, annotationTrack = c("grid", "name"), link.visible = T, grid.col = colours_diagram)




# Dorsal-Ventral axis patterning

## Load libraries
library(Seurat)
library(ggplot2)
library(gridExtra)
library(axis.projection)


## Manual DV annotation
se2 <- readRDS(file = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/r_results/se_SC_cropped_aligned_76sections.rds")
se_annot <- ManualAnnotation(se2)
DV_annot <- as.data.frame(se_annot$labels)
DV_annot <- DV_annot %>% rename(axis_annot = "se_annot$labels")
se_new <- AddMetaData(se2, DV_annot)

## Normalization, PCA, Harmony
se_norm <- SCTransform(se_new, vars.to.regress = "section.name")
se_norm <- RunPCA(se_norm)
se_harmony <- RunHarmony(se_norm, reduction = "pca", group.by.vars = "section.name", assay.use = "SCT")
se_harmony <- FindNeighbors(se_harmony, reduction = "harmony", dims = 1:30)
se_harmony <- FindClusters(se_harmony, verbose= F, resolution = 0.6)
se_harmony <- RunUMAP(se_harmony, reduction = "harmony", dims = 1:30, spread = 1, min.dist = 0.1, n.neighbors = 20)

## Save object
saveRDS(se_harmony, file = "/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/r_results/se_SC_annotated_dorsal_ventral_harmony.rds")

## Load obejct
se_axis <- readRDS("/home/st-analysis_home/zaneta.andrusivova/projects/hdca_spinalcord/r_results/se_SC_annotated_dorsal_ventral_harmony.rds")
head(se_axis@meta.data[c("labels","axis_annot")])


se_axis@meta.data["DV"] <- gsub("[0-9]*_","",se_axis@meta.data$axis_annot,perl = T)
Seurat::DefaultAssay(object = se_axis) <- "SCT"


## Add PCW
se_axis$age <- "unknown"

se_axis$age[grep("V19N18-116", se_axis$section.name)] <- "W12"
se_axis$age[grep("V10T03-302", se_axis$section.name)] <- "W5"
se_axis$age[grep("V10T03-290_A1", se_axis$section.name)] <- "W5"
se_axis$age[grep("V10T03-290_B1", se_axis$section.name)] <- "W5"
se_axis$age[grep("V10T03-290_C1", se_axis$section.name)] <- "W12"
se_axis$age[grep("V10T03-290_D1", se_axis$section.name)] <- "W12"
se_axis$age[grep("V10T03-288", se_axis$section.name)] <- "W9"
se_axis$age[grep("V10F24-107", se_axis$section.name)] <- "W8"
se_axis$age[grep("V10F24-106", se_axis$section.name)] <- "W12"


## DV plots
ab.column <- "DV"
a.label <- "D"
b.label <- "V"
split.on <- "labels"


se_axis <- axis.projection(se_axis,
                           ab.column = ab.column,
                           a.label = a.label,
                           b.label = b.label,
                           split.on = split.on,
                           include.orthogonal = T)



features <- c("PBX3", "RUSC2", "ZIC2","FGFBP3")
grobs <- lapply(features,
                function(feature){plot.axis.projection(se_axis,
                                                       feature,
                                                       ab.column,
                                                       a.label,
                                                       b.label,
                                                       show.data =F,
                                                       normalize = T,
                                                       title = feature,
                                                       split.on = "age")})


## Plot features
grid.arrange(grobs = grobs, ncol = 2)


## Subset based on PCW and anatom.region
W8_axis <- SubsetSTData(se_axis, age == "W8")
W8_subset <- SubsetSTData(W8_axis, expression = labels %in% c("section_11", "section_9", "section_15", "section_14", "section_18", "section_23", "section_26"))

W8_subset$anatom.region[grep("section_11", W8_subset$labels)] <- "Thoracic"
W8_subset$anatom.region[grep("section_14", W8_subset$labels)] <- "Thoracic"
W8_subset$anatom.region[grep("section_9", W8_subset$labels)] <- "Lumbar"
W8_subset$anatom.region[grep("section_18", W8_subset$labels)] <- "Lumbar"
W8_subset$anatom.region[grep("section_23", W8_subset$labels)] <- "Lumbar"
W8_subset$anatom.region[grep("section_15", W8_subset$labels)] <- "Cervical"
W8_subset$anatom.region[grep("section_26", W8_subset$labels)] <- "Cervical"


## Plots
g1 <- plot.axis.projection(W8_subset,
                           "LHX2",                          
                           ab.column,
                           a.label,
                           b.label,
                           show.data =F,
                           normalize = T,
                           color.color = scale_color_hue(),
                           fill.color = scale_color_hue(),
                           split.on = "anatom.region", 
                           title = "LHX2")


g1 + scale_color_manual(values = c("darkgrey"))







