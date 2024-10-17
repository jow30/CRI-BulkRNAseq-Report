if(!dir.exists(file.path(params$report_out_dir, "4.Functional_analysis"))) dir.create(file.path(params$report_out_dir, "4.Functional_analysis"))

if (params$species == "human") {
  gsorDb <- "org.Hs.eg.db"
  keggOrg <- 'hsa'
  msigdbr_species <- "Homo sapiens"
} else if (params$species == "mouse") {
  gsorDb <- "org.Mm.eg.db"
  keggOrg <- 'mmu'
  msigdbr_species <- "Mus musculus"
} else if (params$species == "rat") {
  gsorDb <- "org.Rn.eg.db"
  keggOrg <- 'rno'
  msigdbr_species <- "Rattus norvegicus"
} else {
  stop("Species not supported. ")
}

if(params$ora_go){
  cat("<h3>Gene Ontology (GO)</h3>")
  
  enrichGO_res <- list()
  go_cat <- c('BP','CC','MF')
  
  if (params$ora_all) {
    # run GO enrichment analysis
    enrichGO_res$all <- lapply(go_cat, function(x) clusterProfiler::compareCluster(geneCluster = deGenes, fun = 'enrichGO', keyType = 'SYMBOL', OrgDb = gsorDb, ont = x, readable = TRUE)) %>% setNames(go_cat)
    # plot the enriched GO terms
    lapply(go_cat, function(x) gs_dotplot(enrichGO_res$all[[x]], params$ora_fdr_thres, paste0("Top 10 enriched GO ", x, " terms per cluster"))) %>% invisible()
    # save the enriched GO terms
    lapply(go_cat, function(x) lapply(my_contrasts, function(y) write_csv(subset(enrichGO_res$all[[x]]@compareClusterResult, Cluster==y & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_GO_", x, "_", y, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), ".csv")))) %>% invisible()
  } else {
    enrichGO_res$down <- lapply(go_cat, function(x) clusterProfiler::compareCluster(geneCluster = downGenes, fun = 'enrichGO', keyType = 'SYMBOL', OrgDb = gsorDb, ont = x, readable = TRUE)) %>% setNames(go_cat)
    enrichGO_res$up <- lapply(go_cat, function(x) clusterProfiler::compareCluster(geneCluster = upGenes, fun = 'enrichGO', keyType = 'SYMBOL', OrgDb = gsorDb, ont = x, readable = TRUE)) %>% setNames(go_cat)
    
    lapply(go_cat, function(x) lapply(c("up", "down"), function(y) gs_dotplot(enrichGO_res[[y]][[x]], params$ora_fdr_thres, paste0("Top 10 enriched GO ", x, " terms per contrast (", y, "-regulated)")))) %>% invisible()
    
    lapply(go_cat, function(x) lapply(my_contrasts, function(y) write_csv(subset(enrichGO_res$down[[x]]@compareClusterResult, Cluster==y & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_GO_", x, "_", y, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_down.csv")))) %>% invisible()
    lapply(go_cat, function(x) lapply(my_contrasts, function(y) write_csv(subset(enrichGO_res$up[[x]]@compareClusterResult, Cluster==y & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_GO_", x, "_", y, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_up.csv")))) %>% invisible()
  }
}

if(params$ora_kegg){
  cat("<h3>KEGG</h3>")
  
  enrichKEGG_res <- list()
  
  if (params$ora_all) {
    # map gene names to entrez ids
    deEntrez <- lapply(deGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    # run KEGG enrichment analysis
    enrichKEGG_res$all <- clusterProfiler::compareCluster(geneCluster = deEntrez, fun = 'enrichKEGG', organism = keggOrg)
    # plot the enriched KEGG pathways
    gs_dotplot(enrichKEGG_res$all, params$ora_fdr_thres, "Top 10 enriched KEGG pathways per contrast")
    # map entrez ids back to gene names
    enrichKEGG_res$all@compareClusterResult$geneID <- lapply(enrichKEGG_res$all@compareClusterResult$geneID, function(x) unlist(strsplit(x, "/")) %>% bitr(fromType = "ENTREZID", toType = "SYMBOL", OrgDb = gsorDb) %>% dplyr::select(SYMBOL) %>% unlist() %>% paste(collapse = "/")) %>% unlist()
    # save the enriched KEGG pathways
    lapply(my_contrasts, function(x) write_csv(subset(enrichKEGG_res$all@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_KEGG_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), ".csv"))) %>% invisible()
  }else{
    downEntrez <- lapply(downGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    upEntrez <- lapply(upGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    
    enrichKEGG_res$down <- clusterProfiler::compareCluster(geneCluster = downEntrez, fun = 'enrichKEGG', organism = keggOrg)
    enrichKEGG_res$up <- clusterProfiler::compareCluster(geneCluster = upEntrez, fun = 'enrichKEGG', organism = keggOrg)
    
    lapply(c("up", "down"), function(x) gs_dotplot(enrichKEGG_res[[x]], params$ora_fdr_thres, paste0("Top 10 enriched KEGG pathways per contrast (", x, "-regulated)"))) %>% invisible()
    
    enrichKEGG_res$down@compareClusterResult$geneID <- lapply(enrichKEGG_res$down@compareClusterResult$geneID, function(x) unlist(strsplit(x, "/")) %>% bitr(fromType = "ENTREZID", toType = "SYMBOL", OrgDb = gsorDb) %>% dplyr::select(SYMBOL) %>% unlist() %>% paste(collapse = "/")) %>% unlist()
    enrichKEGG_res$up@compareClusterResult$geneID <- lapply(enrichKEGG_res$up@compareClusterResult$geneID, function(x) unlist(strsplit(x, "/")) %>% bitr(fromType = "ENTREZID", toType = "SYMBOL", OrgDb = gsorDb) %>% dplyr::select(SYMBOL) %>% unlist() %>% paste(collapse = "/")) %>% unlist()
    
    lapply(my_contrasts, function(x) write_csv(subset(enrichKEGG_res$down@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_KEGG_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_down.csv"))) %>% invisible()
    lapply(my_contrasts, function(x) write_csv(subset(enrichKEGG_res$up@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_KEGG_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_up.csv"))) %>% invisible()
  }
}

if(params$ora_reactome){
  cat("<h3>Reactome</h3>")
  
  enrichReact_res <- list()
  
  if (params$ora_all) {
    # map gene names to entrez ids
    deEntrez <- lapply(deGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    # run Reactome enrichment analysis
    enrichReact_res$all <- clusterProfiler::compareCluster(geneCluster = deEntrez, fun = 'enrichPathway', organism = params$species, readable = TRUE)
    # plot the enriched Reactome pathways
    gs_dotplot(enrichReact_res$all, params$ora_fdr_thres, "Top 10 enriched Reactome pathways per contrast")
    # save the enriched Reactome pathways
    lapply(my_contrasts, function(x) write_csv(subset(enrichReact_res$all@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_Reactome_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), ".csv"))) %>% invisible()
  }else{
    downEntrez <- lapply(downGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    upEntrez <- lapply(upGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    
    enrichReact_res$down <- clusterProfiler::compareCluster(geneCluster = downEntrez, fun = 'enrichPathway', organism = params$species, readable = TRUE)
    enrichReact_res$up <- clusterProfiler::compareCluster(geneCluster = upEntrez, fun = 'enrichPathway', organism = params$species, readable = TRUE)
    
    lapply(c("up", "down"), function(x) gs_dotplot(enrichReact_res[[x]], params$ora_fdr_thres, paste0("Top 10 enriched Reactome pathways per contrast (", x, "-regulated)"))) %>% invisible()
    
    lapply(my_contrasts, function(x) write_csv(subset(enrichReact_res$down@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_Reactome_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_down.csv"))) %>% invisible()
    lapply(my_contrasts, function(x) write_csv(subset(enrichReact_res$up@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_Reactome_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_up.csv"))) %>% invisible()
  }
}

if(params$ora_msigdb){
  cat("<h3>MSigDB</h3>")
  
  # view all available gene sets: print(msigdbr_collections(), n=24)
  # see the full explanation of gene sets at https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
  # fetch gene sets and genes
  gs <- msigdbr::msigdbr(species = msigdbr_species, category = params$ora_msigdbr_category, subcategory = params$ora_msigdbr_subcategory) %>% dplyr::select(gs_name, entrez_gene)

  enrichMSigDB_res <- list()
  
  if (params$ora_all) {
    # map gene names to entrez ids
    deEntrez <- lapply(deGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    # run MSigDB enrichment analysis
    enrichMSigDB_res$all <- clusterProfiler::compareCluster(geneCluster = deEntrez, fun = 'enricher', TERM2GENE = gs)
    # plot the enriched MSigDB pathways
    gs_dotplot(enrichMSigDB_res$all, params$ora_fdr_thres, paste("Top 10 enriched MSigDB", paste(c(params$ora_msigdbr_category, params$ora_msigdbr_subcategory), collapse = " "), "pathways per contrast"))
    # map entrez ids back to gene names
    enrichMSigDB_res$all@compareClusterResult$geneID <- lapply(enrichMSigDB_res$all@compareClusterResult$geneID, function(x) unlist(strsplit(x, "/")) %>% bitr(fromType = "ENTREZID", toType = "SYMBOL", OrgDb = gsorDb) %>% dplyr::select(SYMBOL) %>% unlist() %>% paste(collapse = "/")) %>% unlist()
    # save the enriched MSigDB pathways
    lapply(my_contrasts, function(x) write_csv(subset(enrichMSigDB_res$all@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_MSigDB_", paste(c(params$ora_msigdbr_category, params$ora_msigdbr_subcategory), collapse = "_"),  "_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), ".csv"))) %>% invisible()
  }else{
    downEntrez <- lapply(downGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    upEntrez <- lapply(upGenes, function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = gsorDb)$ENTREZID)
    
    enrichMSigDB_res$down <- clusterProfiler::compareCluster(geneCluster = downEntrez, fun = 'enricher', TERM2GENE = gs)
    enrichMSigDB_res$up <- clusterProfiler::compareCluster(geneCluster = upEntrez, fun = 'enricher', TERM2GENE = gs)
    
    lapply(c("up", "down"), function(x) gs_dotplot(enrichMSigDB_res[[x]], params$ora_fdr_thres, paste0("Top 10 enriched MSigDB ", paste(c(params$ora_msigdbr_category, params$ora_msigdbr_subcategory), collapse = " "), " pathways per contrast (", x, "-regulated)"))) %>% invisible()
    
    lapply(my_contrasts, function(x) write_csv(subset(enrichMSigDB_res$down@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_MSigDB_", paste(c(params$ora_msigdbr_category, params$ora_msigdbr_subcategory), collapse = "_"),  "_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_down.csv"))) %>% invisible()
    lapply(my_contrasts, function(x) write_csv(subset(enrichMSigDB_res$up@compareClusterResult, Cluster==x & p.adjust<params$ora_fdr_thres), file = paste0(file.path(params$report_out_dir, "4.Functional_analysis"), "/ORA_MSigDB_", paste(c(params$ora_msigdbr_category, params$ora_msigdbr_subcategory), collapse = "_"),  "_", x, "_padj", sub("0\\.", "", as.character(params$ora_fdr_thres)), "_up.csv"))) %>% invisible()
  }
}

if(params$gsea_msigdb){
  cat("<h2>4.2 GSEA</h2>")
  cat("Gene set enrichment analysis or [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) uses a ranked list of all protein-coding genes, in this case based on log2(FC) for a specified comparison. For a gene set of interest, GSEA calculates an enrichment score by running through the ranked gene list and examining whether the gene set is over-represented at the top or bottom of the gene list, which can be interpreted as enrichment in a particular phenotypic condition/group for the comparison of interest. <br><br>")
  cat("In this analysis, we analyzed KEGG pathways. Please contact our staffs to analyze more pathways. <br><br>")
  cat("Note that the GSEA plots shown below only present the top 5 pathways. And the bubble plots only present the top 10. To plot the pathway of your interest, please further contact our staffs. <br><br>")

  lapply(my_contrasts, function(x) run_GSEA(x, res[[x]], msigdbr_species, params$gsea_msigdbr_category, params$gsea_msigdbr_subcategory, params$gsea_fdr_thres, file.path(params$report_out_dir, "4.Functional_analysis")))
}

