#' Perform differential expression analysis on pseudobulk data
#'
#' @param pseudobulk_data A PseudobulkData object
#' @param design Formula or design matrix for the analysis
#' @param method Analysis method ("DESeq2" or "edgeR")
#' @param contrast Character vector specifying the contrast of interest
#' @param min_counts Minimum count threshold for filtering genes (default: 10)
#' @param min_samples Minimum number of samples where counts > min_counts (default: 3)
#' @param lfc_threshold Log2 fold change threshold (default: 0)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param verbose Logical indicating whether to print progress messages
#'
#' @return A data frame containing differential expression results
#' @export
#'
#' @examples
#' \dontrun{
#' # Using DESeq2
#' de_results <- differential_expression(
#'   pseudobulk_data,
#'   design = ~ condition,
#'   contrast = c("condition", "treated", "control")
#' )
#' }
differential_expression <- function(pseudobulk_data,
  design,
  method = "DESeq2",
  contrast = NULL,
  min_counts = 10,
  min_samples = 3,
  lfc_threshold = 0,
  padj_threshold = 0.05,
  verbose = TRUE) {

# Input validation
if (!inherits(pseudobulk_data, "PseudobulkData")) {
stop("Input must be a PseudobulkData object")
}

if (!method %in% c("DESeq2", "edgeR")) {
stop("Method must be either 'DESeq2' or 'edgeR'")
}

# Extract count matrix and metadata
counts <- pseudobulk_data@counts
metadata <- pseudobulk_data@metadata

# Filter low-expressed genes
if (verbose) message("Filtering low-expressed genes...")
keep_genes <- rowSums(counts >= min_counts) >= min_samples
counts <- counts[keep_genes, ]

if (verbose) {
message(sprintf("Retained %d out of %d genes after filtering",
sum(keep_genes),
length(keep_genes)))
}

# Perform differential expression analysis based on chosen method
if (method == "DESeq2") {
results <- run_deseq2_analysis(
counts = counts,
metadata = metadata,
design = design,
contrast = contrast,
lfc_threshold = lfc_threshold,
verbose = verbose
)
} else {
results <- run_edger_analysis(
counts = counts,
metadata = metadata,
design = design,
contrast = contrast,
lfc_threshold = lfc_threshold,
verbose = verbose
)
}

# Add gene symbols if available
results <- add_gene_annotations(results, pseudobulk_data)

# Filter results based on significance
significant <- which(results$padj < padj_threshold)

if (verbose) {
message(sprintf("Found %d significantly differential genes (padj < %g)",
length(significant),
padj_threshold))
}

# Add volcano plot status
results$diffexpressed <- "NS"
results$diffexpressed[results$padj < padj_threshold & results$log2FoldChange > lfc_threshold] <- "UP"
results$diffexpressed[results$padj < padj_threshold & results$log2FoldChange < -lfc_threshold] <- "DOWN"

return(results)
}

#' Run DESeq2 analysis
#' @keywords internal
run_deseq2_analysis <- function(counts, metadata, design, contrast, lfc_threshold, verbose) {
if (!requireNamespace("DESeq2", quietly = TRUE)) {
stop("DESeq2 package is required for this analysis")
}

if (verbose) message("Running DESeq2 analysis...")

# Create DESeqDataSet object
dds <- DESeq2::DESeqDataSetFromMatrix(
countData = counts,
colData = metadata,
design = design
)

# Run DESeq2
dds <- DESeq2::DESeq(dds, quiet = !verbose)

# Get results with contrast
if (!is.null(contrast)) {
res <- DESeq2::results(dds, 
contrast = contrast,
lfcThreshold = lfc_threshold,
alpha = 0.05)
} else {
res <- DESeq2::results(dds,
lfcThreshold = lfc_threshold,
alpha = 0.05)
}

# Convert to data frame and clean up
results <- as.data.frame(res)
results$gene_id <- rownames(results)

# Add mean counts
results$baseMean <- round(results$baseMean, 2)

return(results)
}

#' Run edgeR analysis
#' @keywords internal
run_edger_analysis <- function(counts, metadata, design, contrast, lfc_threshold, verbose) {
if (!requireNamespace("edgeR", quietly = TRUE)) {
stop("edgeR package is required for this analysis")
}

if (verbose) message("Running edgeR analysis...")

# Create DGEList object
dge <- edgeR::DGEList(counts = counts)

# Calculate normalization factors
dge <- edgeR::calcNormFactors(dge)

# Create design matrix if formula provided
if (inherits(design, "formula")) {
design_matrix <- model.matrix(design, data = metadata)
} else {
design_matrix <- design
}

# Estimate dispersions
dge <- edgeR::estimateDisp(dge, design_matrix)

# Fit model
fit <- edgeR::glmQLFit(dge, design_matrix)

# Define contrast if provided
if (!is.null(contrast)) {
if (is.character(contrast)) {
contrast_matrix <- makeContrasts(contrasts = contrast, 
      levels = design_matrix)
} else {
contrast_matrix <- contrast
}
qlf <- edgeR::glmQLFTest(fit, contrast = contrast_matrix)
} else {
qlf <- edgeR::glmQLFTest(fit)
}

# Get results
results <- as.data.frame(edgeR::topTags(qlf, n = Inf))
results$gene_id <- rownames(results)

# Rename columns to match DESeq2 output
colnames(results)[colnames(results) == "logFC"] <- "log2FoldChange"
colnames(results)[colnames(results) == "PValue"] <- "pvalue"
colnames(results)[colnames(results) == "FDR"] <- "padj"

return(results)
}

#' Add gene annotations to results
#' @keywords internal
add_gene_annotations <- function(results, pseudobulk_data) {
# Check if gene symbols are available in the pseudobulk data object
if ("gene_symbols" %in% names(pseudobulk_data@parameters)) {
gene_symbols <- pseudobulk_data@parameters$gene_symbols
results$gene_symbol <- gene_symbols[results$gene_id]
}

return(results)
}

#' Volcano plot helper function
#' @export
plot_volcano <- function(de_results, 
label_genes = NULL,
max_overlaps = 10) {
if (!requireNamespace("ggplot2", quietly = TRUE)) {
stop("ggplot2 is required for plotting")
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
stop("ggrepel is required for gene labels")
}

# Create basic volcano plot
p <- ggplot2::ggplot(de_results, 
ggplot2::aes(x = log2FoldChange, 
    y = -log10(padj),
    color = diffexpressed)) +
ggplot2::geom_point(size = 1) +
ggplot2::scale_color_manual(values = c("DOWN" = "blue", 
           "UP" = "red", 
           "NS" = "grey")) +
ggplot2::theme_minimal()

# Add labels if specified
if (!is.null(label_genes)) {
de_results$gene_label <- ifelse(de_results$gene_id %in% label_genes |
      de_results$gene_symbol %in% label_genes,
    de_results$gene_symbol, "")

p <- p + ggrepel::geom_text_repel(
ggplot2::aes(label = gene_label),
max.overlaps = max_overlaps
)
}

return(p)
}