if(!dir.exists(file.path(params$report_out_dir, "3.Differential_expression_analysis"))) dir.create(file.path(params$report_out_dir, "3.Differential_expression_analysis"))

if(params$de_method=="DESeq2"){
  cat("We used the standard differential expression analysis steps that were wrapped up in the DESeq2::DESeq() function. ")
  my_contrasts <- lapply(params$comp_pair, function(x) paste0(x[1],"-", x[2])) %>% unlist()
  
  dds <- DESeq2::DESeq(dds)
  res <- lapply(params$comp_pair, function(x) DESeq2::results(dds, contrast=c("group", x[1], x[2]))) %>% setNames(my_contrasts)
  resSig <- lapply(res, function(x) subset(x, padj < params$fdr_thres & abs(log2FoldChange) > log2(params$fc_thres))) %>% setNames(my_contrasts)
  resSig_up <- lapply(resSig, function(x) subset(x, log2FoldChange > log2(params$fc_thres))) %>% setNames(my_contrasts)
  resSig_dn <- lapply(resSig, function(x) subset(x, log2FoldChange < -log2(params$fc_thres))) %>% setNames(my_contrasts)
  
  summary_df <- lapply(resSig, function(x) table(factor(x$log2FoldChange > 0, levels = c(FALSE, TRUE))) )
  summary_df <- do.call(rbind, summary_df) 
  colnames(summary_df) <- c("LFC < 0 (down)","LFC > 0 (up)")
  rownames(summary_df) <- lapply(params$comp_pair, function(x) paste0(x[[1]], "-", x[[2]])) %>% unlist()
  
  res <- lapply(res, function(x) x %>% as.data.frame() %>% dplyr::rename("logFC"="log2FoldChange") %>% dplyr::rename("adj.P.Val"="padj"))
  vp_list <- lapply(my_contrasts, function(x) plot_volcano(x, res[[x]], params$fdr_thres, log2(params$fc_thres)))
  
  resSig_tb <- lapply(resSig, function(x) x %>% as.data.frame() %>% dplyr::rename("logFC"="log2FoldChange") %>% dplyr::rename("adj.P.Val"="padj"))
  tb_list <- lapply(seq(my_contrasts), function(x) table_degs(x, my_contrasts[[x]], sampleinfo, log2.cpm.norm, resSig_tb[[x]], params$fdr_thres, params$fc_thres, file.path(params$report_out_dir, "3.Differential_expression_analysis")))
  lapply(my_contrasts, function(x) plot_heatmap(x, sampleinfo, assay(vsd), resSig_tb[[x]], file.path(params$report_out_dir, "3.Differential_expression_analysis")))
  
}else if(params$de_method=="limma-voom"){
  cat("First, precision weights were applied to TMM-normalized gene counts based on within-group sample-level variance and gene-level mean-variance trends using [VOOM](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29). <br>")
  cat("Second, the count data was fitted to a gene-wise linear model with group status as a coefficient, and contrasts were specified for comparisons of interest. <br>")
  cat("Third, the [Limma](https://academic.oup.com/nar/article/43/7/e47/2414268) empirical Bayes method was used to estimate the posterior odds of differential expression after adjusting for gene-level posterior residual standard deviations. <br>")
  
  # Set up contrasts for pairwise comparisons
  my_contrasts <- lapply(params$comp_pair, function(x) paste0("group", x[1],"-group", x[2])) %>% unlist()

  contrast.matrix <- limma::makeContrasts(contrasts = my_contrasts, levels = colnames(design))
  colnames(contrast.matrix) <- gsub("^group","",colnames(contrast.matrix))
  colnames(contrast.matrix) <- gsub("-group","-",colnames(contrast.matrix))
  
  # Estimate the fold changes and standard errors by fitting a linear model for each gene.
  fit <- limma::lmFit(v, design)
  
  # Compute estimated coefficients and standard errors for a given set of contrasts.
  fit <- limma::contrasts.fit(fit, contrasts = contrast.matrix)
  
  # Apply empirical Bayes smoothing to the standard errors.
  fit <- limma::eBayes(fit, robust = TRUE)
  
  # limma::plotSA(fit, main="Final model: Mean-variance trend")
  
  dt <- limma::decideTests(fit, p.value=params$fdr_thres, lfc = log2(params$fc_thres))
  summary_df <- summary(dt) %>% unclass() %>% as.data.frame() %>% t() %>% as.data.frame()
  
  my_contrasts <- colnames(dt)
  res <- lapply(my_contrasts, function(x) limma::topTable(fit, adjust = "BH", coef = x, number = Inf, sort.by = "P")) %>% setNames(my_contrasts)
  
  resSig <- lapply(my_contrasts, function(x) res[[x]] %>% dplyr::filter(adj.P.Val <= params$fdr_thres, abs(logFC) > log2(params$fc_thres))) %>% setNames(my_contrasts)
  resSig_up <- lapply(my_contrasts, function(x) subset(resSig[[x]], logFC > log2(params$fc_thres))) %>% setNames(my_contrasts)
  resSig_dn <- lapply(my_contrasts, function(x) subset(resSig[[x]], logFC < -log2(params$fc_thres))) %>% setNames(my_contrasts)
  
  vp_list <- lapply(my_contrasts, function(x) plot_volcano(x, res[[x]], params$fdr_thres, log2(params$fc_thres)))
  tb_list <- lapply(1:ncol(dt), function(x) table_degs(x, my_contrasts[[x]], sampleinfo, log2.cpm.norm, resSig[[x]], params$fdr_thres, params$fc_thres, file.path(params$report_out_dir, "3.Differential_expression_analysis")))
  lapply(my_contrasts, function(x) plot_heatmap(x, sampleinfo, log2.cpm.norm, resSig[[x]], file.path(params$report_out_dir, "3.Differential_expression_analysis")))
}
# list available heatmaps in folder
hm_files <- list.files(file.path(params$report_out_dir, "3.Differential_expression_analysis"), pattern = "^heatmap_", full.names = TRUE)
hm_files <- gsub(".*/heatmap_", "", gsub(".png", "", hm_files))

# prepare genes for over-representation analyses
deGenes <- lapply(resSig, rownames)
upGenes <- lapply(resSig_up, rownames)
downGenes <- lapply(resSig_dn, rownames)

# Remove list elements if the row number is 0
deGenes <- deGenes[sapply(deGenes, length) >= 20]
upGenes <- upGenes[sapply(upGenes, length) >= 20]
downGenes <- downGenes[sapply(downGenes, length) >= 20]
