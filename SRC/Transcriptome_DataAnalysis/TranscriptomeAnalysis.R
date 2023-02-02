source("Main.R")


samps <-  readRDS(paste0(dataDir, "/RDSFiles/allSamplesReplicatesMerged.rds"))
IG_genes <- read.csv(paste0(dataDir,"IG_Genes.txt"), sep = '\t', stringsAsFactors = F)

df=data.frame( sampleNames=sampleNames,
               sampleAttributes=sampleAttributesPrint,
               donorID=as.factor(donorID), stringsAsFactors = F)


for(i in 1:length(samps)){
  serObj <- samps[[i]]
  print(serObj$orig.ident)
  print(serObj$stim)
  print(dim(serObj@assays$RNA))
  samps[[i]] <- serObj[rownames(serObj[["RNA"]]) %ni% c(as.character(IG_genes$Gene.name),
                                                        "CH17-262H11.1", "AC233755.1","CH17-212P11.4", "AC233755.2", "CH17-224D4.1"),]
}



immune.anchors <- FindIntegrationAnchors(object.list = samps, k.filter = NA, k.score=60)
saveRDS(immune.anchors, paste0(dataDir, "/RDSFiles/integratedAnchors_1.rds"))
immune.anchors <- readRDS(paste0(dataDir, "/RDSFiles/integratedAnchors_1.rds"))
immune.combined <- IntegrateData(anchorset = immune.anchors,
                                 dims = 1:100)

saveRDS(immune.combined, paste0(dataDir, "/RDSFiles/immuneCombined_0.rds"))

#immune.combined <- readRDS(paste0(dataDir, "/RDSFiles/immuneCombined.rds")) 
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 100, verbose = FALSE)
Idents(immune.combined) <- "stim"

immune.combined <- FindNeighbors( immune.combined, dims = 1:100)
immune.combined <- FindClusters( immune.combined, resolution = 0.2)

saveRDS(immune.combined, paste0(dataDir, "/RDSFiles/immuneCombined_1.rds"))
immune.combined <- RunTSNE( immune.combined, dims = 1:100)
immune.combined <- RunUMAP( immune.combined, dims = 1:100)

saveRDS(immune.combined, paste0(dataDir, "/RDSFiles/immuneCombined_2.rds"))