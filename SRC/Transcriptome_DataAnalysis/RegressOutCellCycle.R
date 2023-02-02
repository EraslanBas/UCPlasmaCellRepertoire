source("Main.R")
samps <-  readRDS(paste0(dataDir, "/RDSFiles/allSamplesReplicatesMerged.rds"))
IG_genes <- read.csv(paste0(dataDir,"IG_Genes.txt"), sep = '\t', stringsAsFactors = F)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

markerGenes <- list()

for(i in 1:length(samps)){
  serObj <- samps[[i]]
  print(serObj$orig.ident)
  print(serObj$stim)
  print(dim(serObj@assays$RNA))
  serObj <- serObj[rownames(serObj[["RNA"]]) %ni% IgGenes,]
  serObj <- ScaleData(serObj, verbose = FALSE)
  serObj <- RunPCA(serObj, npcs = 100, verbose = FALSE)
 
  serObj <- CellCycleScoring(serObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  RidgePlot(serObj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
  
  hh <- RunPCA(serObj, features = c(s.genes, g2m.genes), npcs=40)
  DimPlot(hh)
  
  serObj <- FindNeighbors( serObj, dims = 1:100)
  serObj <- FindClusters( serObj, resolution = 0.3)

  serObj <- RunTSNE( serObj, dims = 1:100)

  allCellNames <- colnames(serObj@assays$RNA)
  sampleIdentifier <- sapply(allCellNames, function(x){strsplit(x,"_")[[1]][2]})
  serObj$sampleName <- sampleIdentifier

  if(length(unique(serObj$seurat_clusters)) > 1){
    obj.markers <- FindAllMarkers(serObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    obj.markers <- obj.markers[obj.markers$p_val_adj < 0.1,]
    top10 <- obj.markers  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

    markerGenes <- lappend(markerGenes, top10)
  }else{
    markerGenes <- lappend(markerGenes, "None")
  }

  samps[[i]] <- serObj
}

saveRDS(samps, paste0(dataDir, "/RDSFiles/cellCycleRegressedTrans.rds"))
saveRDS(markerGenes, paste0(dataDir, "/RDSFiles/cellCycleRegressedMarkerGenes.rds"))
