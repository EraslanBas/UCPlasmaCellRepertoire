source("Main.R")
reticulate::use_python("/home/beraslan/anaconda3/bin/python", required = T)


samps <-  readRDS(paste0(dataDir, "/RDSFiles/allSamplesReplicatesMerged.rds"))


sAt <- unique(samples[,c("SampleAttributes", "CanonicalSampleName")])

perTypeIntegratedSamples <- list()
allObjectMarkers <- list()

for(elem in unique(sAt$SampleAttributes)){
    sAtTemp <- sAt[sAt$SampleAttributes == elem,]
    sampsTemp <- samps[sAtTemp$CanonicalSampleName]
    
    for(i in 1:length(sampsTemp)){
      serObj <- sampsTemp[[i]]
      print(unique(serObj$stim))
      print(dim(serObj@assays$RNA))
      serObj <- serObj[rownames(serObj[["RNA"]]) %ni% IgGenes,]
      
      allCellNames <- colnames(serObj@assays$RNA)
      sampleIdentifier <- sapply(allCellNames, function(x){strsplit(x,"_")[[1]][2]})
      serObj$sampleName <- sampleIdentifier
      sampsTemp[[i]] <- serObj
    }
    
    features <- SelectIntegrationFeatures(object.list = sampsTemp)
    
    print(str(features))
    print(features)
    
    sampsTemp <- lapply(X = sampsTemp, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE, npcs = 40)
      x <- RunUMAP(x, dims = 1:40)
      x <- RunTSNE(x, dims = 1:40)
    })
    
   
    anchors <- FindIntegrationAnchors(object.list = sampsTemp,
                                      reduction = "cca",
                                      dims = 1:40,
                                      normalization.method="LogNormalize",
                                      scale = F)
    
    sampsTemp.integrated <- IntegrateData(anchorset = anchors, dims = 1:40)
    
    sampsTemp.integrated <- ScaleData(sampsTemp.integrated, verbose = FALSE)
    sampsTemp.integrated <- RunPCA(object = sampsTemp.integrated, verbose = FALSE)
    sampsTemp.integrated <- RunUMAP(sampsTemp.integrated, dims = 1:40)
    sampsTemp.integrated <- RunTSNE(sampsTemp.integrated, dims = 1:40)
    
    sampsTemp.integrated <- FindNeighbors(sampsTemp.integrated, reduction = "pca", dims = 1:40)
    sampsTemp.integrated <- FindClusters(sampsTemp.integrated, resolution = 0.2)
    obj.markersTemp <- FindAllMarkers(sampsTemp.integrated,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)
    
    perTypeIntegratedSamples <- lappend(perTypeIntegratedSamples, sampsTemp.integrated)
    allObjectMarkers <- lappend(allObjectMarkers, obj.markersTemp)
}

names(perTypeIntegratedSamples) <- unique(sAt$SampleAttributes)
names(allObjectMarkers) <- unique(sAt$SampleAttributes)

saveRDS(perTypeIntegratedSamples, paste0(dataDir, "/RDSFiles/perTypeIntegratedSamples_12292019.rds"))
saveRDS(allObjectMarkers, paste0(dataDir, "/RDSFiles/perTypeSamplesObjectMarkers_12292019.rds"))

