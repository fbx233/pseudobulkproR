#' Create pseudobulk data from single-cell data
#'
#' @param data A Seurat object (v4 or v5)
#' @param group_by Character vector specifying the column names in meta.data to group by
#' @param method Character string specifying the aggregation method ("sum" or "mean")
#' @param assay Character string specifying which assay to use (default: "RNA")
#' @param sol Character string specifying which data to use (default: "counts")
#' @param min_cells Integer specifying minimum number of cells required per group (default: 10)
#' @param verbose Logical indicating whether to print progress messages
#'
#' @return A PseudobulkData object containing aggregated data and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # For Seurat object
#' pbulk <- create_pseudobulk(seurat_obj, group_by = c("sample", "celltype"))
#' }
create_pseudobulk <- function(data, 
  group_by, 
  method = "sum",
  assay = "RNA",
  sol = "counts",
  min_cells = 10,
  verbose = TRUE) {

# Input validation
if (!inherits(data, "Seurat")) {
stop("Input must be a Seurat object")
}

if (!all(group_by %in% colnames(data@meta.data))) {
stop("All group_by variables must be present in meta.data")
}

if (!method %in% c("sum", "mean")) {
stop("Method must be either 'sum' or 'mean'")
}

# Detect Seurat version
seurat_version <- package_version(data@version)
is_v5 <- seurat_version >= package_version("5.0.0")

# Helper function to get data based on Seurat version
get_data_matrix <- function(seurat_obj, assay, sol) {
if (is_v5) {
# Seurat V5: Use LayerData
return(LayerData(object = seurat_obj, 
assay = assay, 
layer = sol))
} else {
# Seurat V4: Use GetAssayData
return(GetAssayData(object = seurat_obj, 
assay = assay, 
slot = sol))
}
}

# Get the expression matrix
if (verbose) message("Extracting expression data...")
expr_matrix <- get_data_matrix(data, assay, sol)

# Create grouping factor
if (verbose) message("Creating grouping variable...")
group_cols <- data@meta.data[, group_by, drop = FALSE]
group_factor <- do.call(paste, c(group_cols, list(sep = "_")))

# Get cell counts per group
cell_counts <- table(group_factor)
valid_groups <- names(cell_counts)[cell_counts >= min_cells]

if (length(valid_groups) == 0) {
stop("No groups meet the minimum cell count threshold")
}

if (verbose) {
message(sprintf("Found %d valid groups with >= %d cells", 
length(valid_groups), 
min_cells))
}

# Filter for valid groups
keep_cells <- group_factor %in% valid_groups
expr_matrix <- expr_matrix[, keep_cells]
group_factor <- group_factor[keep_cells]

# Perform aggregation
if (verbose) message("Aggregating data...")

# Convert sparse matrix to dense if necessary
if (inherits(expr_matrix, "dgCMatrix")) {
expr_matrix <- as.matrix(expr_matrix)
}

# Split matrix by groups
split_cols <- split(seq_len(ncol(expr_matrix)), group_factor)

# Aggregate using specified method
if (method == "sum") {
pseudobulk_matrix <- vapply(split_cols, function(idx) {
rowSums(expr_matrix[, idx, drop = FALSE])
}, FUN.VALUE = numeric(nrow(expr_matrix)))
} else { # method == "mean"
pseudobulk_matrix <- vapply(split_cols, function(idx) {
rowMeans(expr_matrix[, idx, drop = FALSE])
}, FUN.VALUE = numeric(nrow(expr_matrix)))
}

# Create metadata for pseudobulk samples
if (verbose) message("Creating metadata...")

metadata <- do.call(rbind, lapply(colnames(pseudobulk_matrix), function(group) {
parts <- strsplit(group, "_")[[1]]
data.frame(
sample_id = group,
setNames(as.list(parts), group_by),
n_cells = cell_counts[group],
stringsAsFactors = FALSE
)
}))

# Create parameters list
parameters <- list(
seurat_version = as.character(seurat_version),
aggregation_method = method,
assay = assay,
sol = sol,
min_cells = min_cells,
group_by = group_by
)

# Create and return PseudobulkData object
if (verbose) message("Creating PseudobulkData object...")

new("PseudobulkData",
counts = pseudobulk_matrix,
metadata = metadata,
parameters = parameters)
}

#' Helper function to validate PseudobulkData object
#'
#' @param object PseudobulkData object
#' @return logical TRUE if valid
#' @keywords internal
validate_pseudobulk_data <- function(object) {
# Check dimensions match
if (ncol(object@counts) != nrow(object@metadata)) {
return("Number of samples in counts matrix doesn't match metadata rows")
}

# Check sample IDs match
if (!all(colnames(object@counts) == object@metadata$sample_id)) {
return("Sample IDs in counts matrix don't match metadata")
}

# Check required parameters
required_params <- c("seurat_version", "aggregation_method", "group_by")
if (!all(required_params %in% names(object@parameters))) {
return("Missing required parameters")
}

TRUE
}

#' Show method for PseudobulkData objects
#'
#' @param object PseudobulkData object
#' @export
setMethod("show", "PseudobulkData", function(object) {
cat("PseudobulkData object\n")
cat("Dimensions:", nrow(object@counts), "genes x", ncol(object@counts), "samples\n")
cat("Aggregation method:", object@parameters$aggregation_method, "\n")
cat("Grouped by:", paste(object@parameters$group_by, collapse = ", "), "\n")
cat("Seurat version:", object@parameters$seurat_version, "\n")
})