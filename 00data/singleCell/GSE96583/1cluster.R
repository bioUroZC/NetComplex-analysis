

rm(list = ls())
library(Seurat)
library(Matrix)
library(dplyr)



dir_path <- "/proj/c.zihao/work2/1pathway/2pertu/GSE96583/"

subfolders <- list.dirs(dir_path, recursive = FALSE, full.names = TRUE)
subfolders <- subfolders[grepl("GSM", basename(subfolders))]
subfolders

samples1 <- data.frame(
  type = rep("10x_dir", length(subfolders)),
  path = subfolders,
  stringsAsFactors = FALSE)

samples1$sample_id <- basename(samples1$path)

head(samples1)

samples <- samples1




out_dir <- "/proj/c.zihao/work2/1pathway/2pertu/GSE96583/out/"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


objs <- list()

for (i in seq_len(nrow(samples))) {
  print(i)
  sample <- samples$sample_id[i]
  print(sample)
  type <- samples$type[i]
  path <- samples$path[i]
  counts <- Read10X(data.dir = path)
  obj <- CreateSeuratObject(counts = counts, 
                            project = sample, 
                            min.cells = 3, 
                            min.features = 200)
  obj <- RenameCells(obj, add.cell.id = sample)
  obj$sample_id <- sample
  objs[[sample]] <- obj
}



mt_pattern <- "^MT-"

for (nm in names(objs)) {
  cat("QC + Normalize:", nm, "\n")
  obj <- objs[[nm]]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj <- subset(
    obj,
    subset = nFeature_RNA >= 200 &
      nFeature_RNA <= 6000 &
      nCount_RNA   >= 800 &
      percent.mt   <= 25
  )
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  sample_dir <- file.path(out_dir, nm)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  mat_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  mat_data   <- GetAssayData(obj, assay = "RNA", layer = "data")
  writeMM(mat_counts, file.path(sample_dir, "counts.mtx"))
  writeMM(mat_data,   file.path(sample_dir, "logNorm.mtx"))
  writeLines(rownames(mat_counts), file.path(sample_dir, "genes.txt"))
  writeLines(colnames(mat_counts), file.path(sample_dir, "cells.txt"))
  write.table(
    as.matrix(mat_data),
    file = file.path(sample_dir, "logNorm.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
  write.table(
    as.matrix(mat_counts),
    file = file.path(sample_dir, "counts.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
  cat("  Exported to:", sample_dir, "\n")
  cat("  counts dim:", dim(mat_counts), " nnz:", length(mat_counts@x), "\n")
  cat("  data   dim:", dim(mat_data),   " nnz:", length(mat_data@x), "\n")
  objs[[nm]] <- obj
}

sapply(objs, ncol)

head(rownames(objs[[1]]))




sce.anchors <- FindIntegrationAnchors(
  object.list = objs,
  dims = 1:30,
  reduction = "cca"
)

sce.CCA <- IntegrateData(
  anchorset = sce.anchors,
  dims = 1:30,
  new.assay.name = "CCA"
)

DefaultAssay(sce.CCA) <- "CCA"

sce.CCA <- ScaleData(
  sce.CCA,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)

sce.CCA <- RunPCA(
  sce.CCA,
  npcs = 30,
  verbose = FALSE
)



pdf(
  file = file.path(out_dir, "ElbowPlot_PCA.pdf"),
  width = 6,
  height = 5
)
ElbowPlot(sce.CCA, ndims = 30)
dev.off()



seleNumber <- 13
set.seed(123)
sce.CCA <- FindNeighbors(sce.CCA, dims = 1:seleNumber)
sce.CCA <- FindClusters(sce.CCA, resolution = 0.8)
sce.CCA <- RunUMAP(sce.CCA, dims = 1:seleNumber)

pdf(
  file = file.path(out_dir, "UMAP_clusters.pdf"),
  width = 6,
  height = 5
)
DimPlot(
  sce.CCA,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.3
)
dev.off()


Idents(sce.CCA) <- "seurat_clusters"

cell_cluster <- data.frame(
  cell = colnames(sce.CCA),
  cluster = as.character(Idents(sce.CCA)),
  sample_id = sce.CCA$sample_id,
  stringsAsFactors = FALSE
)

write.csv(
  cell_cluster,
  file = file.path(out_dir, "cell_cluster_sample.csv"),
  row.names = FALSE
)



DefaultAssay(sce.CCA) <- "RNA"

sce.CCA[["RNA"]] <- JoinLayers(sce.CCA[["RNA"]])

Idents(sce.CCA) <- "seurat_clusters"

markers <- FindAllMarkers(
  sce.CCA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


head(markers)

saveRDS(
  markers,
  file = file.path(out_dir, "markers.rds")
)


top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50)%>%
  arrange(cluster, desc(avg_log2FC))


top_markers <- as.data.frame(top_markers)
top_markers

write.csv(
  top_markers,
  file = file.path(out_dir, "top_markers_by_cluster.csv"),
  row.names = FALSE
)

saveRDS(sce.CCA, file = file.path(out_dir, "sce_CCA.rds"))

pdf(file.path(out_dir, "UMAP_by_sample.pdf"), width = 6, height = 5)
DimPlot(sce.CCA, group.by = "sample_id", pt.size = 0.3)
dev.off()

