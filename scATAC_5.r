library(ArchR)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(dplyr)
library(chromVARmotifs)
source('/share2/pub/zhouyj/zhouyj/organoid/data/function/heatmap_zscore.r')
source('/share2/pub/zhouyj/zhouyj/organoid/data/function/RenameGenesSeurat.r')

referance_path <- "/share2/pub/zhouyj/zhouyj/organoid/data/referance/analysis/eye/drived_data/eye.rds"
addArchRGenome('hg38')
inputFiles <- c('/share2/pub/zhouyj/zhouyj/organoid/data/eye/E_MTAB_12714/GB2_ATAC_S1_2/outs/fragments.tsv.gz',
               '/share2/pub/zhouyj/zhouyj/organoid/data/eye/E_MTAB_12714/GB3_ATAC_S3_5/outs/fragments.tsv.gz',
               '/share2/pub/zhouyj/zhouyj/organoid/data/eye/E_MTAB_12714/GB4_ATAC_S6_7/outs/fragments.tsv.gz')
names(inputFiles) <-  c('GB2_ATAC_S1_2','GB3_ATAC_S3_5','GB4_ATAC_S6_7')

out_path <- '/share2/pub/zhouyj/zhouyj/organoid/data/eye/analysis/HRO-5'
setwd(glue::glue(out_path,'/drived_data'))

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Hsapiens.UCSC.hg38)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, OrgDb = org.Hs.eg.db)

geneAnnotation <- createGeneAnnotation(
  TSS = geneAnnotation$TSS,
  exons = geneAnnotation$exons,
  genes = geneAnnotation$genes
)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, # 这个参数不需要过高，后续可以调整
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation()
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)

projHeme2 <- filterDoublets(projHeme1)

projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

###添加z-score提取代码

###

seRNA <- readRDS(referance_path)
seRNA <- RenameGenesSeurat(obj=seRNA,newnames=seRNA@assays[["RNA"]]@meta.features[["gene_symbols"]])
rownames(seRNA[["RNA"]]@meta.features) <- rownames(seRNA[["RNA"]])

projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "annotation2",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

p1 <- plotEmbedding(
    projHeme2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    embedding = 'UMAP'
)

projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
#    groupList = groupList,
    groupRNA = "annotation2",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM,1,which.max)]
projHeme3$annotation2 <- mapLabels(projHeme3$Clusters,newLabels = labelNew,oldLabels = labelOld)
p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "annotation2")

projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "annotation2")

X_pca <- projHeme4@embeddings@listData[["UMAP"]]@listData[["df"]]
X_pca2 <- data.matrix(X_pca) %>% as.data.frame
rownames(X_pca2) <- NULL
readr::write_tsv(X_pca2,'./umap_embedding.tsv',col_names = F)

meta_data <- cbind(rownames(projHeme4), projHeme4$annotation2)
meta_data <- as.data.frame(meta_data)
colnames(meta_data) <- c("Cell","annotation2")
readr::write_tsv(meta_data,'./meta.tsv')
readr::write_tsv(as.data.frame(table(meta_data$annotation2)),'./cell_count.tsv',col_names = F)


pathToMacs2 <- findMacs2()
projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "annotation2", 
    pathToMacs2 = pathToMacs2
)

projHeme5 <- addPeakMatrix(projHeme4)

markersPeaks <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "PeakMatrix", 
    groupBy = "annotation2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
peakList <- getMarkers(markersPeaks,cutOff = 'FDR > -1')
for (i in names(peakList)){
peakList[[i]][,'cell_type'] <- i
    }
peakList_full<-Reduce(rbind,peakList)
peakList_full <- as.data.frame(peakList_full)
colnames(peakList_full)<- c('seqnames','idx','start','end','avg_log2FC','p_val_adj','MeanDiff','cell_type')
readr::write_tsv(peakList_full,'./all.peaks.celltype.tsv')

markersGS <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "annotation2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmap_zscore <- heatmap_zscore(markersGS)
readr::write_tsv(as.data.frame(heatmap_zscore),'./heatmap_zscore.tsv')

markerList <- getMarkers(markersGS,cutOff = 'FDR > -1')

for (i in names(markerList)){
markerList[[i]][,'cell_type'] <- i
    }
markerList_full<-Reduce(rbind,markerList)
readr::write_tsv(as.data.frame(markerList_full),'./all.markers.celltype.tsv')

m_top50 <- lapply(unique(markerList_full$cell_type),function(Cell_type){
        markers_select <- markerList_full[grep(Cell_type,markerList_full$cell_type),]
        return(head(markers_select,n = 50))
    })
m_top50_full <- Reduce(rbind,m_top50)
readr::write_tsv(as.data.frame(m_top50_full),'./top50.tsv')

projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
motif_sub_List <- lapply(seq_len(ncol(motifsUp)), function(x){
  sub_motif <- motifsUp[, x]
  df <- data.frame(TF = rownames(sub_motif), mlog10Padj = assay(sub_motif)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  df$cell_type <- colnames(sub_motif)
  return(df)
})

motif_enrich <- Reduce(rbind,motif_sub_List)
readr::write_tsv(motif_enrich,'./motif_enrich.tsv')

if("Motif" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}

projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-ProjHeme5", load = FALSE)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = FALSE)
readr::write_tsv(as.data.frame(plotVarDev),'./motif_deviations.tsv')
















