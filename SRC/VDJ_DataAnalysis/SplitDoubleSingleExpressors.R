source("Main.R")


for(sampleName in sampleNames){
  print(sampleName)
  allContigs <- read.csv(paste0(dataDir,"/VDJ_Files/",sampleName,"/filtered_contig_annotations.csv"), stringsAsFactors = F)
  allContigs <- allContigs[allContigs$chain %in% c("IGH","IGL","IGK"),]
  
  allFasta <- seqinr::read.fasta(paste0(dataDir,"/VDJ_Files/",sampleName,"/filtered_contig.fasta"), as.string = T)
  
  allContigs <- as.data.table(allContigs)
  allContigs[,numberOfRowsPerCell := .(.N),by=barcode]
  allContigsSingle <- allContigs[allContigs$numberOfRowsPerCell ==2,]
  allContigsSingle[,numberOfChainTypes := length(unique(chain)),by=barcode]
  allContigsSingle <- allContigsSingle[allContigsSingle$numberOfChainTypes==2,]
  
  allContigsSingle <- allContigsSingle[order(allContigsSingle$barcode, allContigsSingle$chain),]
  
  allContigsSingle <- allContigsSingle[, ChainPairs := do.call(paste, as.list(chain)), by = .(barcode)]
  allContigsSingle$ChainPairs <- sapply(allContigsSingle$ChainPairs, function(x){str_replace_all(x, pattern=" ", replacement="_")})
  allContigsSingle <- allContigsSingle[allContigsSingle$ChainPairs %in% c("IGH_IGK", "IGH_IGL"),]
  
  allContigsMulti <- allContigs[allContigs$numberOfRowsPerCell > 2,] 
  
  fastaSingle <- allFasta[names(allFasta) %in% allContigsSingle$contig_id]
  fastaMulti <- allFasta[names(allFasta) %in% allContigsMulti$contig_id]
  
  write.csv(allContigsSingle, 
            paste0(dataDir,"/VDJ_Files/",sampleName,"/filtered_single_contig_annotations.csv"))
  write.csv(allContigsMulti, 
            paste0(dataDir,"/VDJ_Files/",sampleName,"/filtered_multi_contig_annotations.csv"))
  
  write.fasta(sequences = getSequence(fastaSingle, as.string = TRUE),
              names = names(fastaSingle),
              file.out = paste0(dataDir,"/VDJ_Files/",sampleName,"/filtered_single_contig.fasta"))
 
  write.fasta(sequences = getSequence(fastaMulti, as.string = TRUE),
              names = names(fastaMulti),
              file.out = paste0(dataDir,"/VDJ_Files/",sampleName,"/filtered_multi_contig.fasta"))
}
  