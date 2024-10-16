library(rmarkdown)
library(knitr)
library(markdown)
library(tidyverse)
library(DT) 
library(rtracklayer)
library(tximport)
library(edgeR)
library(DESeq2)
library(SummarizedExperiment)
library(patchwork)
library(limma)
library(plotly)
library(htmltools)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr) # access to msigdb collections directly within R
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots
library(ReactomePA)
library(pathview)

get_sampleinfo <- function(file){
  sampleinfo <- read.delim(file)
  sampleinfo <- sampleinfo[order(sampleinfo$group),]
  rownames(sampleinfo) <- sampleinfo$sample
  return(sampleinfo)
}

remove_samples <- function(sample_to_remove, sampleinfo){
  if (length(sample_to_remove) > 0) {
    if (length(setdiff(sample_to_remove, sampleinfo$sample)) > 0) {
      stop("Samples to remove are not found. Please check whether the samples to remove exist in the metadata file.")
    } else {
      sampleinfo <- sampleinfo[!sampleinfo$sample %in% sample_to_remove,]
      cat("<span style='color:red;'>The following samples were removed from the analysis:", paste(sample_to_remove, collapse = ", "), "</span><br>")
      return(sampleinfo)
    }
  }else{
    return(sampleinfo)
  }
}

plot_cpm_density <- function(myDGEList, subtitle, log2.cpm.cutoff){
  # Get the log2 cpm matrix and change to tibble (in tidyverse)
  log2.cpm <- edgeR::cpm(myDGEList, log=TRUE) %>% as_tibble(rownames = "gene_name")
  # Pivot cpm dataframe
  log2.cpm.pivot <- tidyr::pivot_longer(log2.cpm,
                                        cols = -1, # select sample columns to pivot into longer format
                                        names_to = "samples", # name of the new variable (column) storing sample names
                                        values_to = "expression") # name of the new variable (column) storing sample values
  p1 <- ggplot(log2.cpm.pivot, aes(x=expression,group=samples,color=samples)) +
    geom_density(show.legend = FALSE) +
    labs(y="density", x = "log2 CPM") +
    theme_bw()+
    geom_vline(xintercept = log2.cpm.cutoff, color="gray", linetype="dashed")
  p2 <- ggplot(log2.cpm.pivot) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    geom_hline(yintercept = log2.cpm.cutoff, color="gray", linetype="dashed") +
    stat_summary(fun = "median",
                 geom = "point",
                 shape = 95,
                 size = 10,
                 color = "black",
                 show.legend = FALSE) +
    labs(y="log2 CPM", x = "sample") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0))
  p<- p1 + p2 +
    plot_layout(widths = c(1, 2)) +
    plot_annotation(title = "Log2 Counts per Million (CPM)",
                    subtitle = subtitle,
                    caption = "Grey dash lines indicate the log2(cpm) cutoff. Solid lines indicate the median.") &
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

plot_pca <- function(pca_attr, pca_df, color_by, shape_by, top_var){
  sample <- factor(pca_attr$sample)
  if(color_by %in% colnames(pca_attr)) {
    color_attr <- factor(pca_attr[[color_by]])
  } else{ stop("Invalid sample attribute provided for colors in PCA plot. ") }
  
  if(is.numeric(top_var) && top_var>0){
    most_variable_genes <- rownames(pca_df)[order(rowVars(pca_df), decreasing=TRUE)] %>% head(top_var)
    pca.res <- prcomp(t(pca_df[most_variable_genes,]), scale.=F, retx=T)
    cat("The top", top_var, "most variable genes were used in PCA. <br>")
  }else{ 
    pca.res <- prcomp(t(pca_df), scale.=F, retx=T) 
    cat("All pre-processed genes were used in PCA. <br>")
  }
  
  pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
  pc.per<-round(pc.var/sum(pc.var)*100, 1)
  pca.res.df <- as_tibble(pca.res$x)
  
  p1 <- ggplot(pca.res.df) +
    aes(x = PC1, y = PC2, label = sample, color = color_attr) +
    geom_point(size = 2) +
    ggtitle("PCA Plot") +
    xlab(paste0("PC1 (", pc.per[1], "%", ")")) +
    ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
    labs(color = color_by) +
    coord_fixed() +
    theme_bw()
  
  if (!is.null(shape_by)) {
    if (shape_by %in% colnames(pca_attr)) {
      shape_attr <- factor(pca_attr[[shape_by]])
      p1 <- ggplot(pca.res.df) +
        aes(x = PC1, y = PC2, label = sample, color = color_attr, shape = shape_attr) +
        geom_point(size = 2) +
        ggtitle("PCA Plot") +
        xlab(paste0("PC1 (", pc.per[1], "%", ")")) +
        ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
        labs(color = color_by, shape = shape_by) +
        coord_fixed() +
        theme_bw()
    } else{
      stop("Invalid sample attribute provided for shapes in PCA plot. ")
    }
  }
  
  scree.df <- data.frame(PC=as.character(1:length(pc.per)),Eigen=pc.per)
  scree.df$PC <- factor(scree.df$PC, levels = scree.df$PC)
  p2 <- ggplot(scree.df, aes(PC,Eigen)) +
    geom_bar(stat="identity", fill = "gray")+
    theme_bw()+
    ggtitle("Scree Plot") +
    xlab("Principle component (PC)") +
    ylab("Explained variation (%)") +
    labs(caption=paste0("produced on ", Sys.time())) +
    geom_text(aes(label = Eigen), vjust = -5, colour = "black", size = 3)
  
  return(list(pca=p1, scree=p2))
}

plot_volcano <- function(title, res, fdr.thres, lfc.thres) {
  res <- res %>% tibble::rownames_to_column("gene_name") %>% as.data.frame()
  res_up <- subset(res, adj.P.Val < fdr.thres & logFC > lfc.thres)
  res_down <- subset(res, adj.P.Val < fdr.thres & logFC < -lfc.thres)
  # Prepare the base plot
  vplot <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), text = gene_name)) +
    geom_point(size = 2) +
    geom_hline(yintercept = -log10(fdr.thres), linetype = "longdash", colour = "grey", size = 1) +
    geom_vline(xintercept = lfc.thres, linetype = "longdash", colour = "#BE684D", size = 1) +
    geom_vline(xintercept = -lfc.thres, linetype = "longdash", colour = "#2C467A", size = 1) +
    labs(title = title, subtitle = "Volcano plot", caption = paste0("Produced on ", Sys.time())) +
    ylab("-log10(FDR)") +
    xlab("log2(FC)") +
    theme_bw() +
    scale_x_continuous(limits = c(-5, 5), oob = scales::squish)
  
  # Highlight points based on thresholds
  if (nrow(res_up) > 0) { vplot <- vplot + geom_point(data = res_up, color = 'red') }
  if (nrow(res_down) > 0) { vplot <- vplot + geom_point(data = res_down, color = 'blue') }
  
  interactive_plot <- plotly::ggplotly(vplot, tooltip = "text")
  return(interactive_plot)
}

table_degs <- function(table_number, comp, sampleinfo, sample_df, resSig, fdr.thres, fc.thres, outDir) {
  comp_grp <- str_split(comp, "-")[[1]]
  
  # get samples in each group
  sample_grp <- lapply(1:2, function(x) sampleinfo$sample[sampleinfo$group %in% comp_grp[x]])
  
  # get normalized log2 cpm of each sample for all DEGs
  diffGenes <- sample_df[rownames(resSig), unlist(sample_grp)]
  diffGenes.df <- as_tibble(diffGenes, rownames = "gene_name")
  
  # add average log2 cpm of each group to the table
  for (x in c(1,2)) { diffGenes.df <- dplyr::mutate(diffGenes.df, !!paste0(comp_grp[[x]],".avg") := rowMeans(across(sample_grp[[x]]))) }
  diffGenes.df <- diffGenes.df %>% dplyr::mutate_if(is.numeric, round, 2)
  
  # add log2FC and p-value of DEGs to the table
  deg.df <- resSig %>% tibble::rownames_to_column("gene_name") %>% dplyr::select(gene_name, logFC, adj.P.Val)
  diffGenes.df <- dplyr::left_join(diffGenes.df, deg.df, by = "gene_name") %>% dplyr::mutate_if(is.numeric, round, 2) %>% dplyr::arrange(desc(abs(logFC)))
  
  # save the table to a csv file
  write_csv(diffGenes.df, file = paste0(outDir, "/", comp, "_deg_all_fdr_", fdr.thres, "_fc_", fc.thres, ".csv"))
  
  diffGenes_up.df <- diffGenes.df %>% dplyr::filter(logFC > log2(fc.thres)) %>% dplyr::arrange(desc(logFC))
  diffGenes_down.df <- diffGenes.df %>% dplyr::filter(logFC < -log2(fc.thres)) %>% dplyr::arrange(logFC)
  
  write_csv(diffGenes_up.df, file = paste0(outDir, "/", comp, "_deg_up_fdr_", fdr.thres, "_fc_", fc.thres, ".csv"))
  write_csv(diffGenes_down.df, file = paste0(outDir, "/", comp, "_deg_down_fdr_", fdr.thres, "_fc_", fc.thres, ".csv"))
  
  # generate interactive and searchable tables
  interactive_tb <- DT::datatable(
    diffGenes.df,
    extensions = c('KeyTable', "FixedHeader"),
    caption = paste0('Table ', table_number, '. DEGs for ', comp),
    options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))
  ) %>% formatRound(columns = c(2:ncol(diffGenes.df)), digits = 2)
  return(interactive_tb)
}

plot_heatmap <- function(comp, sampleinfo, data, resSig, outDir){
  comp_grp <- str_split(comp, "-")[[1]]
  sample_grp <- lapply(1:2, function(x) sampleinfo$sample[sampleinfo$group %in% comp_grp[x]])
  
  diffGenes <- data[rownames(resSig), unlist(sample_grp)]
  if (dim(diffGenes)[1]>0) {
    myheatcolors <- rev(RColorBrewer::brewer.pal(name="RdBu", n=8))
    clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
    clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
    module.assign <- cutree(clustRows, k=2)
    module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9)
    module.color <- module.color[as.vector(module.assign)]
    grp.color <- c(rep("red3",length(sample_grp[[1]])),rep("blue3",length(sample_grp[[2]])))
    comp_str <- gsub("-", " vs ", gsub("\\.", " ", comp))
    
    png(paste0(outDir, "/heatmap_", comp, ".png"), width = 2000, height = 2000, res = 300, units = "px")
    gplots::heatmap.2(
      diffGenes,
      Rowv=as.dendrogram(clustRows),
      Colv=as.dendrogram(clustColumns),
      RowSideColors=module.color,
      ColSideColors=grp.color,
      col=myheatcolors, scale='row', labRow=rownames(diffGenes),
      density.info="histogram", trace="none",
      main=comp_str, 
      srtCol=45
    )
    dev.off()
  }
}

# function for dot plot of enriched gene sets
gs_dotplot <- function(clusterProfResult, fdr_thres, title) {
  # get the top 10 enriched gene sets in each cluster to plot
  res2plot <- clusterProfResult@compareClusterResult %>% dplyr::filter(p.adjust < fdr_thres) %>% dplyr::group_by(Cluster) %>% dplyr::slice_head(n = 10) %>% dplyr::ungroup()
  if (nrow(res2plot)==0) { return(NULL) }
  # calculate gene ratio
  res2plot$GeneRatio <- sapply(strsplit(res2plot$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  # set the levels of x and y axis
  res2plot$Cluster <- factor(res2plot$Cluster, levels = names(clusterProfResult@geneClusters))
  res2plot$Description <- factor(res2plot$Description, levels = unique(res2plot$Description))
  
  p <- ggplot(res2plot, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = GeneRatio, colour = p.adjust), show.legend = T) + 
    theme_bw() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")
    ) +
    labs(title = title, x = "", y = "", caption = "") +
    scale_colour_gradient(name = "FDR adjusted\n p-value\n", limits = c(0, fdr_thres)) +
    scale_size_continuous(name = "Gene Ratio")
  print(p)
}

run_GSEA <- function(comp, res, species, msigdbr_category, msigdbr_subcategory, fdr_thres, outDir) {
  if (species == 'human') {
    msigdbr_species <- "Homo sapiens"
  } else if (species == "mouse") {
    msigdbr_species <- "Mus musculus"
  } else {
    stop("Species not supported")
  }
  # view all available gene sets: print(msigdbr_collections(), n=24)
  # see the full explanation of gene sets at https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
  # fetch gene sets and genes
  gs <- msigdbr::msigdbr(species = msigdbr_species, category = msigdbr_category, subcategory = msigdbr_subcategory) %>% dplyr::select(gs_name, gene_symbol)
  
  comp_grp <- str_split(comp, "-")[[1]]
  
  #make a gene table with decreasing logFC for this contrast group
  mydata.gsea <- res$logFC
  names(mydata.gsea) <- rownames(res)
  mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
  
  #run GSEA
  myGSEA.res <- clusterProfiler::GSEA(mydata.gsea, pvalueCutoff = fdr_thres, TERM2GENE=gs, verbose=FALSE, eps = 0)
  
  if (nrow(myGSEA.res)>0) {
    #add a variable to this result that matches enrichment direction with phenotype
    myGSEA.df <- myGSEA.res@result %>% mutate(phenotype = case_when(NES > 0 ~ comp_grp[1], NES < 0 ~ comp_grp[2]))
    myGSEA.df <- myGSEA.df %>% arrange(desc(abs(NES)))
    
    write_csv(myGSEA.df, paste0(outDir, '/', comp, '_GSEA_table_', paste(c(msigdbr_category, msigdbr_subcategory), collapse = '_'), '.csv'))
    
    gsea_plot = enrichplot::gseaplot2(
      myGSEA.res,
      geneSetID = myGSEA.df$ID[1:min(5, nrow(myGSEA.df))], #plot the first five pathways only
      base_size = 8,
      rel_heights = c(1.5, 0.5, 0.5),
      title = paste0(comp, ' GSEA ', paste(c(msigdbr_category, msigdbr_subcategory), collapse = ' '), ' top5')
    )
    ggsave(paste0(outDir, '/', comp, '_GSEA_plot_', paste(c(msigdbr_category, msigdbr_subcategory), collapse = '_'), '.png'), gsea_plot, width = 7, height = 5)
    
    #create 'bubble plot' to summarize y signatures across x phenotypes
    bubble_plot <- ggplot(myGSEA.df[1:min(dim(myGSEA.df)[1],10),], aes(x=phenotype, y=ID)) + #plot the first ten pathways
      geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
      scale_color_gradient(low="blue", high="red") +
      theme_bw() +
      labs(title = paste0(comp, ' GSEA'), subtitle = paste(c(msigdbr_category, msigdbr_subcategory, "top10"), collapse = ' '), x = "", y = "", caption = "NES: Normalized Enrichment Score")
    ggsave(paste0(outDir, '/', comp, '_GSEA_bubble_plot_', paste(c(msigdbr_category, msigdbr_subcategory), collapse = '_'), '.png'), bubble_plot, width = 7, height = 5)
    
  }else{
    cat(paste0("No significant enrichemnt for ", comp, "."))
  }
}


