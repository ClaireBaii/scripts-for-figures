#!/usr/bin/env Rscript

# FigS｜rarefied vs non-rarefied 的 PCoA / PERMANOVA 稳健性对照
# R 4.x (你当前 paper_fig_env 的 R/vegan 版本也能跑)
#
# 输入：
# - data/完整数据-微生物.csv
#
# 输出（写入 figure3/）：
# - FigS_rarefaction_robustness.pdf + FigS_rarefaction_robustness.tiff
# - permanova_robustness.csv

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

# ---------- robust path helpers ----------
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    p <- sub("^--file=", "", file_arg[1])
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  # source() case
  for (i in rev(seq_along(sys.frames()))) {
    of <- attr(sys.frames()[[i]], "ofile")
    if (!is.null(of)) return(normalizePath(of, winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

find_project_root <- function(start_dir, max_up = 8) {
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 0:max_up) {
    if (file.exists(file.path(d, "data", "完整数据-微生物.csv"))) return(d)
    if (dir.exists(file.path(d, "data")) && dir.exists(file.path(d, "figure3"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

read_csv_safe <- function(path) {
  for (enc in c("UTF-8-BOM", "UTF-8", "GB18030")) {
    x <- try(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = enc), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
  }
  stop("Failed to read CSV with UTF-8-BOM/UTF-8/GB18030: ", path, call. = FALSE)
}

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop(
      paste0(
        "Missing R packages: ", paste(miss, collapse = ", "), "\n",
        "建议用 conda-forge 在 paper_fig_env 里安装，例如：\n",
        "  mamba install -n paper_fig_env -c conda-forge r-vegan r-ggplot2 r-dplyr r-patchwork\n"
      ),
      call. = FALSE
    )
  }
}

# ---------- grouping ----------
infer_group <- function(sample_id) {
  m <- regexec("^RDYS_(\\d+)-", sample_id)
  mm <- regmatches(sample_id, m)[[1]]
  if (length(mm) >= 2) {
    grp_num <- as.integer(mm[2])
    return(switch(
      as.character(grp_num),
      "1"="Control","2"="Treat1","3"="Treat2","4"="Treat3","5"="Treat4",
      paste0("Group", grp_num)
    ))
  }
  if (grepl("control", sample_id, ignore.case = TRUE)) return("Control")
  m2 <- regexec("treat\\s*([0-9]+)", sample_id, ignore.case = TRUE)
  mm2 <- regmatches(sample_id, m2)[[1]]
  if (length(mm2) >= 2) return(paste0("Treat", mm2[2]))
  "Unknown"
}

# ---------- distance alignment ----------
align_dist_to_meta <- function(dist_obj, meta) {
  lab <- attr(dist_obj, "Labels")
  if (is.null(lab)) lab <- rownames(as.matrix(dist_obj))

  if (!all(meta$Sample %in% lab)) {
    miss <- meta$Sample[!meta$Sample %in% lab]
    stop("Some meta samples not in distance labels: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  m <- as.matrix(dist_obj)
  m2 <- m[meta$Sample, meta$Sample, drop = FALSE]
  as.dist(m2)
}

# ---------- robust extraction of adonis2 results ----------
extract_adonis2_r2_p <- function(adon, term = "Group") {
  df <- as.data.frame(adon)
  cn <- colnames(df)
  rn <- rownames(df)

  # pick row
  term_i <- which(rn == term)
  if (length(term_i) == 0) term_i <- which(!grepl("Residual", rn, ignore.case = TRUE))
  if (length(term_i) == 0) term_i <- 1
  term_i <- term_i[1]

  # pick R2 column
  r2_col <- which(cn == "R2")
  if (length(r2_col) == 0) r2_col <- grep("R2", cn, fixed = TRUE)
  if (length(r2_col) == 0) stop("Cannot find R2 column in adonis2 output.", call. = FALSE)
  r2_col <- r2_col[1]

  # pick p column (vegan may rename to Pr..F.)
  p_col <- which(cn %in% c("Pr(>F)", "Pr..F.", "Pr(>F).", "Pr(>F )"))
  if (length(p_col) == 0) p_col <- grep("^Pr", cn)
  if (length(p_col) == 0) p_col <- grep("p", tolower(cn))
  if (length(p_col) == 0) stop("Cannot find p-value column in adonis2 output.", call. = FALSE)
  p_col <- p_col[1]

  r2 <- suppressWarnings(as.numeric(df[term_i, r2_col]))
  p  <- suppressWarnings(as.numeric(df[term_i, p_col]))

  list(r2 = r2, p = p, df = df, picked = c(row = rn[term_i], r2_col = cn[r2_col], p_col = cn[p_col]))
}

# ---------- compute PCoA + PERMANOVA ----------
compute_pcoa_permanova <- function(dist_obj, meta, permutations, method_label) {
  dist_obj <- align_dist_to_meta(dist_obj, meta)

  pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE)
  coords <- as.data.frame(pcoa$points)
  colnames(coords) <- c("PC1", "PC2")
  coords$Sample <- rownames(coords)
  coords$Group  <- meta$Group[match(coords$Sample, meta$Sample)]

  eig <- pcoa$eig
  eig_pos <- eig[eig > 0]
  var_expl <- eig_pos / sum(eig_pos) * 100
  pc1_pct <- ifelse(length(var_expl) >= 1, var_expl[1], NA_real_)
  pc2_pct <- ifelse(length(var_expl) >= 2, var_expl[2], NA_real_)

  rownames(meta) <- meta$Sample
  adon <- vegan::adonis2(dist_obj ~ Group, data = meta, permutations = permutations)

  ex <- extract_adonis2_r2_p(adon, term = "Group")
  r2 <- ex$r2
  p  <- ex$p

  # if still NA, print debug table for immediate diagnosis
  if (is.na(r2) || is.na(p)) {
    cat("\n[DEBUG] adonis2 table:\n")
    print(ex$df)
    cat("\n[DEBUG] picked row/cols: ",
        ex$picked["row"], " | ",
        ex$picked["r2_col"], " | ",
        ex$picked["p_col"], "\n", sep = "")
  }

  list(
    coords = coords,
    pc1_pct = pc1_pct,
    pc2_pct = pc2_pct,
    permanova_r2 = r2,
    permanova_p = p,
    method_label = method_label
  )
}

# ====================== MAIN ======================
need_pkgs(c("vegan", "ggplot2", "dplyr", "patchwork"))
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

script_path <- get_script_path()
script_dir  <- if (!is.na(script_path)) dirname(script_path) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
repo_root   <- find_project_root(script_dir)

in_csv <- file.path(repo_root, "data", "完整数据-微生物.csv")
out_dir <- file.path(repo_root, "figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# debug prints
cat("Script path : ", script_path, "\n", sep = "")
cat("Script dir  : ", script_dir,  "\n", sep = "")
cat("Repo root   : ", repo_root,   "\n", sep = "")
cat("Input CSV   : ", in_csv,      "\n", sep = "")

if (!file.exists(in_csv)) {
  cat("\n[DEBUG] data dir listing:\n")
  data_dir <- file.path(repo_root, "data")
  if (dir.exists(data_dir)) print(list.files(data_dir))
  stop("Input not found: ", in_csv, call. = FALSE)
}

permutations <- 9999
rarefy_seed <- 1

# 固定稀释深度：设为 NULL 表示用 min depth；你也可以改成 10000 等
rarefaction_depth_fixed <- NULL

raw <- read_csv_safe(in_csv)

name_lower <- tolower(names(raw))
otu_id_col <- which(name_lower %in% c("otu id", "otu_id", "feature id", "feature_id", "asv id", "asv_id"))
otu_id_col <- if (length(otu_id_col) > 0) otu_id_col[1] else 1
tax_col <- which(name_lower == "taxonomy")
tax_col <- if (length(tax_col) > 0) tax_col[1] else NA_integer_

drop_cols <- c(otu_id_col, tax_col)
drop_cols <- drop_cols[!is.na(drop_cols)]
sample_cols <- setdiff(seq_len(ncol(raw)), drop_cols)
sample_ids <- names(raw)[sample_cols]

cat("Detected sample columns = ", length(sample_ids), "\n", sep = "")
cat("First few sample IDs: ", paste(head(sample_ids, 8), collapse = ", "), "\n", sep = "")

otu_mat <- as.matrix(raw[, sample_cols, drop = FALSE])
mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0
otu_mat <- round(otu_mat)
rownames(otu_mat) <- make.unique(as.character(raw[[otu_id_col]]))
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

meta <- data.frame(
  Sample = sample_ids,
  Group  = vapply(sample_ids, infer_group, character(1)),
  stringsAsFactors = FALSE
)

group_order <- c("Control", "Treat1", "Treat2", "Treat3", "Treat4")
meta$Group <- factor(meta$Group, levels = unique(c(group_order, sort(unique(meta$Group)))))

cat("Group counts:\n\n")
print(table(meta$Group))

if (min(table(meta$Group)) < 2) {
  stop("至少有一个组样本数 < 2，PERMANOVA 可能返回 NA。请检查 sample_id 命名与 infer_group 规则。", call. = FALSE)
}

group_colors <- c(
  Control = "#F8766D",
  Treat1  = "#7CAE00",
  Treat2  = "#00BFC4",
  Treat3  = "#00B0F6",
  Treat4  = "#C77CFF"
)

# -------- non-rarefied：相对丰度 Bray --------
col_sums <- colSums(otu_mat)
if (any(col_sums == 0)) stop("Some samples have zero total counts.", call. = FALSE)

otu_rel <- sweep(otu_mat, 2, col_sums, FUN = "/")
otu_rel[is.na(otu_rel)] <- 0
dist_nonref <- vegan::vegdist(t(otu_rel), method = "bray")
res_nonref <- compute_pcoa_permanova(dist_nonref, meta, permutations, "Non-rarefied (relative abundance)")

# -------- rarefied：等深度稀释后 Bray --------
set.seed(rarefy_seed)
sample_counts <- rowSums(t(otu_mat))
min_depth <- min(sample_counts)

depth_use <- if (is.null(rarefaction_depth_fixed)) {
  min_depth
} else {
  if (min_depth < rarefaction_depth_fixed) {
    stop(sprintf("min depth (%d) < fixed rarefaction depth (%d).", min_depth, rarefaction_depth_fixed), call. = FALSE)
  }
  rarefaction_depth_fixed
}

otu_raref <- vegan::rrarefy(t(otu_mat), sample = depth_use)
dist_raref <- vegan::vegdist(otu_raref, method = "bray")
res_raref <- compute_pcoa_permanova(dist_raref, meta, permutations, sprintf("Rarefied (depth=%d)", depth_use))

# -------- write table --------
robust_tbl <- data.frame(
  method = c("non_rarefied", "rarefied"),
  label  = c(res_nonref$method_label, res_raref$method_label),
  R2     = c(res_nonref$permanova_r2, res_raref$permanova_r2),
  p_value= c(res_nonref$permanova_p,  res_raref$permanova_p),
  permutations = permutations,
  rarefaction_depth = c(NA_integer_, depth_use),
  distance = "Bray-Curtis",
  stringsAsFactors = FALSE
)
write.csv(robust_tbl, file.path(out_dir, "permanova_robustness.csv"), row.names = FALSE)

# -------- plotting --------
permanova_text <- function(r2, p, perm) sprintf("PERMANOVA: R² = %.3f, P = %.3g (perm=%d)", r2, p, perm)

make_plot <- function(res, xlim, ylim, show_legend = TRUE) {
  ggplot(res$coords, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3.0, alpha = 0.95) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(
      title = res$method_label,
      subtitle = permanova_text(res$permanova_r2, res$permanova_p, permutations),
      x = sprintf("PC1 (%.1f%%)", res$pc1_pct),
      y = sprintf("PC2 (%.1f%%)", res$pc2_pct),
      color = "Group"
    ) +
    coord_equal(xlim = xlim, ylim = ylim) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      panel.grid.minor = element_blank()
    )
}

all_coords <- rbind(
  res_nonref$coords[, c("PC1", "PC2")],
  res_raref$coords[, c("PC1", "PC2")]
)
x_pad <- diff(range(all_coords$PC1)) * 0.06
y_pad <- diff(range(all_coords$PC2)) * 0.06
xlim <- range(all_coords$PC1) + c(-x_pad, x_pad)
ylim <- range(all_coords$PC2) + c(-y_pad, y_pad)

p1 <- make_plot(res_nonref, xlim, ylim, show_legend = FALSE)
p2 <- make_plot(res_raref,  xlim, ylim, show_legend = TRUE)
fig <- p1 + p2 + patchwork::plot_layout(widths = c(1, 1))

pdf_path  <- file.path(out_dir, "FigS_rarefaction_robustness.pdf")
tiff_path <- file.path(out_dir, "FigS_rarefaction_robustness.tiff")

ggsave(pdf_path, fig, width = 12, height = 5.2, units = "in")
ggsave(tiff_path, fig, width = 12, height = 5.2, units = "in",
       dpi = 600, device = "tiff", compression = "lzw")

message("Done. Outputs written to: ", out_dir)
