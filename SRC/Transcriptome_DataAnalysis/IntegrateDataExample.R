source("Main.R")
library(SeuratData)


library("future")
library(future.apply)

plan("multiprocess", workers = 12)
options(future.globals.maxSize=+Inf)

InstallData("panc8")
data(panc8)
print("11111111111111111")
bm40k.list <- SplitObject(panc8, split.by = "orig.ident")
print("22222222222")

bm40k.list <- bm40k.list[1:5]
print("33333333333")

bm40k.list <- lapply(X = bm40k.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
print("444444444444")

features <- SelectIntegrationFeatures(object.list = bm40k.list)

print("5555555555555")
print(features)

bm40k.list <- lapply(X = bm40k.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 30)
})

print("66666666666666")

anchors <- FindIntegrationAnchors(object.list = bm40k.list, reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:30)

print("7777777777")
saveRDS(anchors, paste0(dataDir, "/RDSFiles/seuratIntegratedAnchorsExample.rds"))

bm40k.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
print("888888888888")

saveRDS(bm40k.integrated, paste0(dataDir, "/RDSFiles/seuratIntegratedSamplesExample.rds"))

