library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)

match_tb <- readr::read_tsv("/share2/pub/zhouyj/zhouyj/organoid/data/table/ENSG2ENTERZ_10X.tsv")
match_features <- read.table("/share2/pub/zhouyj/zhouyj/organoid/data/table/features.tsv",header = F)

match_tb$symbols <- match_features$V2[match(match_tb$ensembl_gene_id,match_features$V1)]

enterz_to_symbol <- function(enterz_all){
        enterz_all <- strsplit(enterz_all,"/")[[1]]
        symbols_all <- match_tb$symbols[match(enterz_all,match_tb$entrezgene_id)]
        return(paste(symbols_all,sep="",collapse="/"))
}

markers <- readr::read_tsv(glue::glue("./all.markers.celltype.tsv"))
mingene = 2
markers$ENSG <- match_features$V1[match(markers$name,match_features$V2)]
markers$enterz <- match_tb$entrezgene_id[match(markers$ENSG,match_tb$ensembl_gene_id)]

enrich_go <- lapply(unique(markers$cell_type),function(cell_type){
        markers_select <- markers[grep(cell_type,markers$cell_type),]
        down_genes <- markers_select %>% filter(Log2FC <= -0.25) %>% filter(FDR < 0.05) %>% pull(name)
        if (length(down_genes)>mingene){
        down_go = clusterProfiler::enrichGO(gene=down_genes, 
                        ont ="ALL", 
                        keyType = "SYMBOL", 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        qvalueCutoff = 0.05, 
                        OrgDb = "org.Hs.eg.db")
        head(down_go)
        down_go_tb <- down_go@result %>% mutate(nlog10p_adj = -log10(p.adjust)) %>% mutate(Type = "DOWN_GO")
        }
        up_genes <- markers_select %>% filter(Log2FC >= 0.25) %>% filter(FDR < 0.05) %>% pull(name)
        if (length(up_genes)>mingene){
        up_go = clusterProfiler::enrichGO(gene=up_genes, 
                        ont ="ALL", 
                        keyType = "SYMBOL", 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        qvalueCutoff = 0.05, 
                        OrgDb = "org.Hs.eg.db")
        up_go_tb <- up_go@result %>% mutate(nlog10p_adj = -log10(p.adjust)) %>% mutate(Type = "UP_GO")
        }
        if(length(down_genes)>mingene & length(up_genes)>mingene ){
        return(rbind(up_go_tb,down_go_tb) %>% mutate(cell_type=cell_type))
        }else if(length(down_genes)>mingene){
        return(down_go_tb%>% mutate(cell_type=cell_type))
        }else if(length(up_genes)>mingene){
        return(up_go_tb %>% mutate(cell_type=cell_type))
        }
})
enrich_go <- Reduce(rbind,enrich_go)
#enrich_go$geneID <- sapply(enrich_go$geneID,name_to_symbol)        
enrich_go %>% readr::write_tsv(glue::glue("./enrich_go.tsv"))

enrich_kegg <- lapply(unique(markers$cell_type),function(cell_type){
        markers_select <- markers[grep(cell_type,markers$cell_type),]
        down_genes <- markers_select %>% filter(Log2FC <= -0.25) %>% filter(FDR < 0.05) %>% pull(enterz)
        if (length(down_genes)>mingene){
        down_kegg = clusterProfiler::enrichKEGG(gene=down_genes, 
                organism = "hsa",
                minGSSize = 10, 
                maxGSSize = 500, 
                pvalueCutoff = 0.05,
                use_internal_data =T
                )
        if(is.null(down_kegg)){
          length(down_genes)<-mingene-1
        }
        else{
          down_kegg_tb <- down_kegg@result %>% mutate(nlog10p_adj = -log10(p.adjust)) %>% mutate(Type = "DOWN_KEGG")
        }
        }
        up_genes <- markers_select %>% filter(Log2FC >= 0.25) %>% filter(FDR < 0.05) %>% pull(enterz)
        if (length(up_genes)>mingene){
        up_kegg = clusterProfiler::enrichKEGG(gene=up_genes, 
                organism = "hsa",
                minGSSize = 10, 
                maxGSSize = 500, 
                pvalueCutoff = 0.05,
                use_internal_data =T 
                )
        if(is.null(up_kegg)){
          length(up_genes)<-mingene-1
        }
        else{
          up_kegg_tb <- up_kegg@result %>% mutate(nlog10p_adj = -log10(p.adjust)) %>% mutate(Type = "UP_KEGG")
        }
        }
        if(length(down_genes)>mingene & length(up_genes)>mingene ){
        return(rbind(up_kegg_tb,down_kegg_tb) %>% mutate(cell_type=cell_type))
        }else if(length(down_genes)>mingene){
        return(down_kegg_tb%>% mutate(cell_type=cell_type))
        }else if(length(up_genes)>mingene){
        return(up_kegg_tb %>% mutate(cell_type=cell_type))
        }
})
enrich_kegg <- Reduce(rbind,enrich_kegg) %>% filter(pvalue < 0.05)
enrich_kegg$geneID <- sapply(enrich_kegg$geneID,enterz_to_symbol)  
enrich_kegg %>% readr::write_tsv(glue::glue("./enrich_kegg.tsv"))