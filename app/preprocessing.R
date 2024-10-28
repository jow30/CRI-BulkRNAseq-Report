# Load in gene count data
files <- file.path(params$nextflow_out_dir,"star_salmon",sampleinfo$sample,"quant.genes.sf") %>% setNames(sampleinfo$sample)
if (!all(file.exists(files))) stop("Samples are not found. Please check whether the samples in the metadata file match the samples provided to the nf-core/rnaseq pipeline")

# Load in gene names
tx2gene<-read.table(file.path(params$nextflow_out_dir,"star_salmon/tx2gene.tsv"), sep="\t")
tx2gene<-tx2gene[,-1] %>% setNames(c("TXNAME","GENENAME"))

# Load in gene annotations
gtf <- rtracklayer::import(params$gtf_file)
geneinfo <- rtracklayer::as.data.frame(gtf)
geneinfo <- geneinfo[geneinfo$type == "gene", c("gene_id", "gene_name", "gene_type")] %>% unique()

if(params$de_method=="DESeq2"){
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene)
  
  if(params$batch_correction) {
    dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = sampleinfo, design = ~ batch + group)
  }else{
    dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = sampleinfo, design = ~ group)
  }
  gene_num_raw <- nrow(dds)

  # if params$log2_cpm_cutoff is not a number, calculate the log2 cpm cutoff; otherwise, use the user-provided value
  if(!is.numeric(params$log2_cpm_cutoff)) {
    log2.cpm.cutoff <- log2(10e6/mean(colSums(counts(dds))))
  }else{
    log2.cpm.cutoff <- params$log2_cpm_cutoff
  }
  p1 <- plot_cpm_density(dds, "Raw", log2.cpm.cutoff)
  
  # pre-filter low count genes
  smallestGroupSize <- min(table(sampleinfo$group))
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  gene_num_filtered <- nrow(dds)
  
  write_csv(as_tibble(counts(dds), rownames='geneName'), file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/gene_counts.csv"))
  cat("Low count genes were filtered out, reducing the number of genes from ", gene_num_raw, " to ", gene_num_filtered, ". ", sep = "")
  cat("Gene counts table with low count genes removed is at `2.Pre-processing_of_raw_reads/gene_counts.csv`. <br><br>", sep = "")
  
  if(params$protein_coding) {
    dds <- dds[rownames(dds) %in% geneinfo$gene_name[geneinfo$gene_type=="protein_coding"],]
    write_csv(as_tibble(counts(dds), rownames='geneName'), file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/protein_coding_gene_counts.csv"))
    cat("We retain only protein-coding genes, reducing the number of genes from ", gene_num_filtered, " to ", nrow(dds), ". ", sep = "")
    cat("Counts of protein-coding genes were saved at `2.Pre-processing_of_raw_reads/protein_coding_gene_counts.csv`. <br><br>", sep = "")
  }
  p2 <- plot_cpm_density(dds, "Filtered", log2.cpm.cutoff)
  
  dds <- DESeq2::estimateSizeFactors(dds)
  dds_norm <- DESeq2::counts(dds, normalized = TRUE)
  write_csv(as_tibble(dds_norm, rownames='geneName'), file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/protein_coding_gene_counts_norm.csv"))
  cat("DESeq2 automatically normalizes count data when it runs differential expression. However, we need to normalize counts aside for certain plots like PCA and heatmaps. We applied a variance stabilizing transformation (VST) to the count data and also normalized with respect to library size when making the plots. ", sep = "")
  cat("The library-size normalized protein-coding gene counts were saved at `2.Pre-processing_of_raw_reads/protein_coding_gene_counts_norm.csv`. <br>", sep = "")
  p3 <- plot_cpm_density(dds_norm, "Filtered, normalized", log2.cpm.cutoff)
  
  log2.cpm.norm <- edgeR::cpm(dds_norm, log=TRUE)
  
  cat("<h3>CPM density plot/violin plot</h3>")
  print(p1)
  print(p2)
  print(p3)

  vsd <- DESeq2::vst(dds)
  pca_attr <- colData(vsd)
  pca_df <- assay(vsd)
  
}else if(params$de_method=="limma-voom"){
  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
  y <- edgeR::DGEList(txi$counts, samples = sampleinfo)
  gene_num_raw <- nrow(y)
  
  if(!is.numeric(params$log2_cpm_cutoff)) {
    # Calculate the log2 cpm cutoff by filterByExpr
    L <- mean(y$samples$lib.size) * 1e-6
    M <- median(y$samples$lib.size) * 1e-6
    log2.cpm.cutoff <- log2(10/M + 2/L)
  }else{
    log2.cpm.cutoff <- params$log2_cpm_cutoff
  }
  
  p1 <- plot_cpm_density(y, "Transcript length bias corrected", log2.cpm.cutoff)
  
  # pre-filter low count genes
  keep <- edgeR::filterByExpr(y, group = y$samples$group)
  y <- y[keep, ]
  gene_num_filtered <- nrow(y)
  
  write_csv(as_tibble(y$counts, rownames='geneName'), file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/gene_counts.csv"))
  cat("Transcript length bias was corrected across samples. Then, low count genes were filtered out, reducing the number of genes from ", gene_num_raw, " to ", gene_num_filtered, ". ", sep = "")
  cat("Table of estimated counts with low-count genes removed is at `2.Pre-processing_of_raw_reads/gene_counts.csv`. <br><br>", sep = "")
  
  if(params$protein_coding) {
    y <- y[rownames(y) %in% geneinfo$gene_name[geneinfo$gene_type=="protein_coding"], ]
    write_csv(as_tibble(y$counts, rownames='geneName'), file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/protein_coding_gene_counts.csv"))
    cat("We retained only protein-coding genes, reducing the number of genes from ", gene_num_filtered, " to ", nrow(y), ".\n", sep = "")
    cat("Counts of protein-coding genes were saved at `2.Pre-processing_of_raw_reads/protein_coding_gene_counts.csv`. <br><br>", sep = "")
  }
  p2 <- plot_cpm_density(y, "Transcript length bias corrected, filtered", log2.cpm.cutoff)

  # Create design matrix
  # refer to https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
  if(params$batch_correction) {
    design <- model.matrix(~0+group+batch, y$samples)
  }else{
    design <- model.matrix(~0+group, y$samples)
  }
    
  # normalize and run voom transformation
  y <- edgeR::calcNormFactors(y, method = "TMM")
  v <- limma::voom(y, design)
  p3 <- plot_cpm_density(y, "Transcript length bias corrected, filtered, and normalized", log2.cpm.cutoff)
  
  # Save normalized counts
  write_csv(as_tibble(y$counts, rownames='geneName'), file = file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/protein_coding_gene_counts_norm.csv"))
  cat("Gene counts were further normalized using the TMM method in edgeR, which were saved at `2.Pre-processing_of_raw_reads/protein_coding_gene_counts_norm.csv`. ", sep = "")

  # Save normalized cpm
  cpm.norm <- edgeR::cpm(y, log=FALSE)
  write_csv(as_tibble(cpm.norm, rownames = "geneName"), file = file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/protein_coding_gene_cpm_norm.csv"))
  cat("Normalized CPM was saved in `2.Pre-processing_of_raw_reads/protein_coding_gene_cpm_norm.csv`. ", sep = "")
  
  # Save log2 normalized cpm
  log2.cpm.norm <- edgeR::cpm(y, log=TRUE)
  write_csv(as_tibble(log2.cpm.norm, rownames = "geneName"), file = file.path(params$report_out_dir, "2.Pre-processing_of_raw_reads/protein_coding_gene_log2_cpm_norm.csv"))
  cat("Normalized log2 CPM was saved in `2.Pre-processing_of_raw_reads/protein_coding_gene_log2_cpm_norm.csv`. <br>", sep = "")
  
  cat("<h3>CPM density plot/violin plot</h3>")
  print(p1)
  print(p2)
  print(p3)
  
  pca_attr <- y$samples
  pca_df <- log2.cpm.norm
  
}else{
  stop("DE test method not supported. Please use DESeq2 or limma-voom. ")
}

