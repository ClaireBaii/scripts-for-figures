# fig5_threshold_sensitivity_table.R
# Purpose: Generate Supplementary Table Sx (threshold sensitivity) for Figure 5

suppressPackageStartupMessages({
  library(igraph)
  # library(Hmisc) # Removed
  # library(dplyr) # Removed
  # library(tidyr) # Removed
  # library(readr) # Removed
})

# Custom rcorr function to replace Hmisc::rcorr (Vectorized)
my_rcorr <- function(x, type = "spearman") {
  x <- as.matrix(x)
  n_samples <- nrow(x)
  
  # Calculate correlation matrix
  # Assuming pairwise complete obs is handled by cor, but for p-values we need sample size.
  # For speed and vectorization, we assume n is constant (no NAs) or accept slight inaccuracy if few NAs.
  # rcorr usually handles pairwise NAs.
  
  # For this specific dataset, we assume no NAs after load.
  # If NAs exist, strict pairwise N is needed.
  # Let's use cor(use="pairwise") and approximate N = n_samples for P-value speed.
  
  r_mat <- cor(x, method = type, use = "pairwise.complete.obs")
  
  # Vectorized P-value calculation
  # t = r * sqrt((n-2)/(1-r^2))
  # Avoid r=1 division by zero
  r2 <- r_mat^2
  r2[r2 > 0.99999999] <- 0.99999999
  
  t_stat <- r_mat * sqrt((n_samples - 2) / (1 - r2))
  p_mat <- 2 * pt(-abs(t_stat), df = n_samples - 2)
  
  # Diagonal cleaning
  diag(r_mat) <- 1
  diag(p_mat) <- 0
  
  list(r = r_mat, P = p_mat)
}

compute_edges <- function(mat, min_threshold = 0) {
  
  cor_res <- my_rcorr(as.matrix(mat), type = "spearman")
  r_mat <- cor_res$r
  p_mat <- cor_res$P
  
  # Filter significantly to reduce size
  # Use min_threshold to only keep relevant edges
  # Only upper triangle
  
  mask <- upper.tri(r_mat) & abs(r_mat) >= min_threshold
  idx <- which(mask, arr.ind = TRUE)
  
  if (nrow(idx) == 0) {
     return(data.frame(from=character(), to=character(), r=numeric(), p=numeric(), 
                       p_adj=numeric(), abs_r=numeric(), sign=character(), stringsAsFactors=FALSE))
  }
  
  # Extract values
  r_vals <- r_mat[mask]
  p_vals <- p_mat[mask]
  rows <- row.names(r_mat)[idx[, 1]]
  cols <- colnames(r_mat)[idx[, 2]]
  
  edges <- data.frame(
    from = rows,
    to   = cols,
    r    = r_vals,
    p    = p_vals,
    stringsAsFactors = FALSE
  )
  
  edges$p_adj <- p.adjust(edges$p, method = "BH")
  edges$abs_r <- abs(edges$r)
  edges$sign  <- ifelse(edges$r >= 0, "pos", "neg")
  edges
}


set.seed(123)

# ---- Config ----
cutoffs <- c(0.5, 0.6, 0.7, 0.8)
p_cutoff <- 0.05

f_exudate <- "data/完整数据-分泌物.csv"
f_microbe <- "data/完整数据-微生物.csv"

f_out <- "figure5/FigS5_threshold_sensitivity_table.csv"
dir.create(dirname(f_out), recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----
load_data_matrix <- function(file_path, type = c("exudate","microbe"), top_n = NULL) {
  type <- match.arg(type)
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))

  # readr 对编码/数字列更稳 (Replaced back to read.csv)
  df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)

  if (type == "exudate") {
    if (!"Name" %in% colnames(df)) stop("Exudate file missing 'Name' column")
    sample_cols <- grep("^YS", colnames(df), value = TRUE)
    if (length(sample_cols) == 0) stop("No sample columns found in exudate data (expect colnames start with 'YS')")
    feat <- df$Name
    feat <- make.unique(as.character(feat))
    mat <- t(as.matrix(df[, sample_cols, drop = FALSE]))
    colnames(mat) <- feat
  } else {
    # microbe：尽量兼容 RDYS 或 YS
    sample_cols <- grep("^(RDYS|YS)", colnames(df), value = TRUE)
    if (length(sample_cols) == 0) stop("No sample columns found in microbe data (expect 'RDYS' or 'YS' prefix)")
    id_col <- colnames(df)[1]
    feat <- df[[id_col]]
    feat <- make.unique(as.character(feat))
    mat <- t(as.matrix(df[, sample_cols, drop = FALSE]))
    colnames(mat) <- feat
  }

  # 确保数值
  storage.mode(mat) <- "numeric"

  # remove zero-variance features
  vars <- apply(mat, 2, var, na.rm = TRUE)
  keep_mask <- vars > 0
  mat <- mat[, keep_mask, drop = FALSE]
  
  # Filter top N by variance if requested
  if (!is.null(top_n) && ncol(mat) > top_n) {
    vars <- vars[keep_mask]
    # sort decreasing
    ord <- order(vars, decreasing = TRUE)
    top_idx <- ord[1:top_n]
    mat <- mat[, top_idx, drop = FALSE]
    cat(sprintf("  Filtered to top %d features by variance.\n", top_n))
  }
  
  mat
}

apply_topk_backbone <- function(edges_keep, k = 2) {
  # edges_keep: from,to,r,abs_r,...
  # 对每个节点保留 abs_r 最大的 k 条边（取并集）
  if (nrow(edges_keep) == 0) return(edges_keep)

  # Create long format manually
  e1 <- data.frame(node = edges_keep$from, idx = seq_len(nrow(edges_keep)), abs_r = edges_keep$abs_r, stringsAsFactors = FALSE)
  e2 <- data.frame(node = edges_keep$to,   idx = seq_len(nrow(edges_keep)), abs_r = edges_keep$abs_r, stringsAsFactors = FALSE)
  long <- rbind(e1, e2)
  
  # Group by node and select top k
  long <- long[order(long$node, -long$abs_r), ]
  # Using split-lapply logic
  # Split by node
  spl <- split(long, long$node)
  
  # Get indices of top k per node (no ties handling matching slice_max default exactly, just head)
  keep_vars <- lapply(spl, function(d) {
     head(d$idx, k)
  })
  
  keep_idx <- unique(unlist(keep_vars))
  
  edges_keep[sort(keep_idx), , drop = FALSE]
}

calc_metrics_one <- function(edges_all, cutoff, p_cutoff = 0.05,
                             use_topk = FALSE, topk = 2) {

  edges_keep <- edges_all[abs(edges_all$r) >= cutoff & edges_all$p_adj < p_cutoff, ]

  if (use_topk) {
    edges_keep <- apply_topk_backbone(edges_keep, k = topk)
  }

  if (nrow(edges_keep) == 0) {
    return(data.frame(
      cutoff = cutoff,
      nodes = 0, edges = 0,
      density = 0,
      avg_degree = 0,
      modularity_Q = 0,
      pos_edges = 0, neg_edges = 0, pos_neg_ratio = NA_real_,
      largest_component = 0,
      avg_clustering = NA_real_,
      n_modules = 0,
      stringsAsFactors = FALSE
    ))
  }

  g <- graph_from_data_frame(edges_keep, directed = FALSE)

  # weights: abs correlation
  E(g)$weight <- abs(edges_keep$r)

  n <- vcount(g); m <- ecount(g)
  dens <- edge_density(g, loops = FALSE)
  avg_deg <- if (n > 0) (2*m / n) else 0

  # Louvain modularity (weighted)
  Q <- 0; n_mod <- 0
  if (n > 1 && m > 0) {
    cl <- cluster_louvain(g, weights = E(g)$weight)
    Q <- modularity(cl, weights = E(g)$weight)
    n_mod <- length(cl)
  }

  pos_n <- sum(edges_keep$r > 0)
  neg_n <- sum(edges_keep$r < 0)
  ratio <- if (neg_n > 0) pos_n / neg_n else NA_real_

  comp <- components(g)
  gcc <- max(comp$csize)

  # average clustering; for tiny graphs may return NaN
  clust <- suppressWarnings(transitivity(g, type = "average"))
  if (is.nan(clust)) clust <- NA_real_

  data.frame(
  cutoff = cutoff,
  nodes = n, edges = m,
  density = dens,
  avg_degree = avg_deg,
  modularity_Q = Q,
    pos_edges = pos_n,
    neg_edges = neg_n,
    pos_neg_ratio = ratio,
    largest_component = gcc,
    avg_clustering = clust,
    n_modules = n_mod,
    stringsAsFactors = FALSE
  )
}

# ---- Main ----
run_one_dataset <- function(name, mat, cutoffs, p_cutoff, use_topk = FALSE, topk = 2) {
  min_thr <- min(cutoffs)
  edges_all <- compute_edges(mat, min_threshold = min_thr)
  res_list <- lapply(cutoffs, function(cut) {
    d <- calc_metrics_one(edges_all, cut, p_cutoff, use_topk, topk)
    d$dataset <- name
    d$method <- ifelse(use_topk, paste0("topk=", topk), "full")
    d
  })
  do.call(rbind, res_list)
}

cat("Loading matrices...\n")
mat_ex <- load_data_matrix(f_exudate, "exudate")
mat_mb <- load_data_matrix(f_microbe, "microbe", top_n = 3000)

cat(sprintf("Exudate: %d samples x %d features\n", nrow(mat_ex), ncol(mat_ex)))
cat(sprintf("Microbe : %d samples x %d features\n", nrow(mat_mb), ncol(mat_mb)))

cat("Computing sensitivity tables...\n")
df_ex <- run_one_dataset("exudate", mat_ex, cutoffs, p_cutoff, use_topk = FALSE)
  df_mb <- run_one_dataset("microbe_genus", mat_mb, cutoffs, p_cutoff, use_topk = TRUE, topk = 2)
  
  final_df <- rbind(df_ex, df_mb)
  # Select columns
  cols_keep <- c("dataset", "method", "cutoff",
                 "nodes", "edges", "density", "avg_degree", "modularity_Q",
                 "pos_edges", "neg_edges", "pos_neg_ratio",
                 "largest_component", "avg_clustering", "n_modules")
  final_df <- final_df[, cols_keep]
  
  write.csv(final_df, f_out, row.names = FALSE)
cat(sprintf("\nDone. Saved: %s\n", f_out))
