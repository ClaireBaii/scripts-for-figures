#!/usr/bin/env Rscript

# FigS: Bray–Curtis distance-to-control
# 输入: data/完整数据-微生物.csv
# 输出: figure3/FigS_distance_to_control.pdf + .tiff, figure3/distance_to_control_stats.csv

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

# ---------------- utils ----------------
get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop(
      paste0(
        "Missing R packages: ", paste(miss, collapse = ", "), "\n",
        "请不要在脚本里 install.packages。\n",
        "建议用 conda-forge 安装，例如：\n",
        "  mamba install -n paper_fig_env -c conda-forge r-vegan r-ggplot2\n"
      ),
      call. = FALSE
    )
  }
}

infer_group <- function(sample_id) {
  # 期望样本名形如：RDYS_1-xxx / RDYS_2-xxx ...
  m <- regexec("^RDYS_(\\d+)-", sample_id)
  mm <- regmatches(sample_id, m)[[1]]
  if (length(mm) >= 2) {
    grp_num <- as.integer(mm[2])
    return(switch(
      as.character(grp_num),
      "1" = "Control",
      "2" = "Treat1",
      "3" = "Treat2",
      "4" = "Treat3",
      "5" = "Treat4",
      paste0("Group", grp_num)
    ))
  }
  if (grepl("control", sample_id, ignore.case = TRUE)) return("Control")
  m2 <- regexec("treat\\s*([0-9]+)", sample_id, ignore.case = TRUE)
  mm2 <- regmatches(sample_id, m2)[[1]]
  if (length(mm2) >= 2) return(paste0("Treat", mm2[2]))
  "Unknown"
}

sig_star <- function(p) {
  if (is.na(p)) return("NA")
  if (p <= 0.0001) return("****")
  if (p <= 0.001)  return("***")
  if (p <= 0.01)   return("**")
  if (p <= 0.05)   return("*")
  "ns"
}

read_csv_safe <- function(path) {
  # 优先 UTF-8-BOM，失败再 UTF-8，再 GB18030（防止中文路径/内容编码）
  for (enc in c("UTF-8-BOM", "UTF-8", "GB18030")) {
    x <- try(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = enc), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
  }
  stop("Failed to read CSV with UTF-8-BOM/UTF-8/GB18030: ", path, call. = FALSE)
}

# ---------------- main ----------------
need_pkgs(c("vegan", "ggplot2"))
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
})

repo_root <- get_repo_root()
in_csv <- file.path(repo_root, "data", "完整数据-微生物.csv")
out_dir <- file.path(repo_root, "figure3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

control_group <- "Control"

if (!file.exists(in_csv)) stop("Input not found: ", in_csv, call. = FALSE)

raw <- read_csv_safe(in_csv)

name_lower <- tolower(names(raw))
otu_id_col <- which(name_lower %in% c("otu id","otu_id","feature id","feature_id","asv id","asv_id"))
otu_id_col <- if (length(otu_id_col) > 0) otu_id_col[1] else 1
tax_col <- which(name_lower == "taxonomy")
tax_col <- if (length(tax_col) > 0) tax_col[1] else NA_integer_

raw <- raw[!is.na(raw[[otu_id_col]]) & raw[[otu_id_col]] != "", , drop = FALSE]
raw <- raw[!duplicated(raw[[otu_id_col]]), , drop = FALSE]

# 样本列：除了 ID 和 taxonomy 之外的全部列
drop_cols <- c(otu_id_col, tax_col)
drop_cols <- drop_cols[!is.na(drop_cols)]
sample_cols <- setdiff(seq_len(ncol(raw)), drop_cols)
sample_ids <- names(raw)[sample_cols]

if (length(sample_ids) < 5) stop("Too few sample columns detected in microbe CSV.", call. = FALSE)

otu_mat <- as.matrix(raw[, sample_cols, drop = FALSE])
mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0
otu_mat <- round(otu_mat)

rownames(otu_mat) <- make.unique(as.character(raw[[otu_id_col]]))
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]

# meta
meta <- data.frame(
  Sample = sample_ids,
  Group  = vapply(sample_ids, infer_group, character(1)),
  stringsAsFactors = FALSE
)
group_order <- c("Control","Treat1","Treat2","Treat3","Treat4")
meta$Group <- factor(meta$Group, levels = unique(c(group_order, sort(unique(meta$Group)))))

# sanity
tab <- table(meta$Group)
if (!control_group %in% names(tab)) stop("No Control samples inferred. Check sample IDs.", call. = FALSE)
if (tab[[control_group]] < 2) stop("Need >=2 Control samples for leave-one-out PairwiseMean.", call. = FALSE)

# relative abundance (features x samples)
col_sums <- colSums(otu_mat)
if (any(col_sums == 0)) stop("Some samples have zero total counts; cannot compute relative abundance.", call. = FALSE)

otu_rel <- sweep(otu_mat, 2, col_sums, FUN = "/")
otu_rel[is.na(otu_rel)] <- 0

otu_rel_t <- t(otu_rel)  # samples x features
stopifnot(nrow(otu_rel_t) == nrow(meta))
rownames(otu_rel_t) <- sample_ids

control_samples <- meta$Sample[meta$Group == control_group]

# Bray–Curtis distances among all samples
bray <- vegan::vegdist(otu_rel_t, method = "bray")
bray_mat <- as.matrix(bray)

# definition 1: PairwiseMean (leave-one-out within Control)
dist_pairwise <- vapply(sample_ids, function(s) {
  cs <- control_samples
  if (s %in% cs) cs <- setdiff(cs, s)
  if (length(cs) == 0) return(NA_real_)
  mean(bray_mat[s, cs], na.rm = TRUE)
}, numeric(1))

# definition 2: Centroid distance (on relative abundances)
control_centroid <- colMeans(otu_rel_t[control_samples, , drop = FALSE])
dist_centroid <- vapply(sample_ids, function(s) {
  x <- otu_rel_t[s, ]
  0.5 * sum(abs(x - control_centroid))  # denom=2 for compositional vectors summing to 1
}, numeric(1))

# long format (base)
dctl_df <- data.frame(
  Sample = meta$Sample,
  Group  = meta$Group,
  Centroid = as.numeric(dist_centroid[meta$Sample]),
  PairwiseMean = as.numeric(dist_pairwise[meta$Sample]),
  stringsAsFactors = FALSE
)

dctl_long <- rbind(
  data.frame(Sample=dctl_df$Sample, Group=dctl_df$Group, Definition="Centroid", Distance=dctl_df$Centroid),
  data.frame(Sample=dctl_df$Sample, Group=dctl_df$Group, Definition="PairwiseMean", Distance=dctl_df$PairwiseMean)
)
dctl_long$Definition <- factor(dctl_long$Definition, levels = c("Centroid","PairwiseMean"))
dctl_long$Group <- factor(dctl_long$Group, levels = levels(meta$Group))

# ---------------- stats ----------------
# KW per Definition
kw_res <- do.call(rbind, lapply(levels(dctl_long$Definition), function(defn) {
  sub <- dctl_long[dctl_long$Definition == defn & !is.na(dctl_long$Distance), ]
  p <- kruskal.test(Distance ~ Group, data = sub)$p.value
  data.frame(Definition=defn, kw_p=p, stringsAsFactors=FALSE)
}))

# pairwise vs Control using Wilcoxon, BH-adjust within each Definition
pair_res <- do.call(rbind, lapply(levels(dctl_long$Definition), function(defn) {
  sub <- dctl_long[dctl_long$Definition == defn & !is.na(dctl_long$Distance), ]
  groups <- setdiff(levels(sub$Group), control_group)
  groups <- groups[!is.na(groups) & groups != "Unknown"]
  out <- lapply(groups, function(g) {
    x <- sub$Distance[sub$Group == control_group]
    y <- sub$Distance[sub$Group == g]
    if (length(x) < 1 || length(y) < 1) {
      return(data.frame(Definition=defn, Group=g, p=NA_real_, stringsAsFactors=FALSE))
    }
    p <- wilcox.test(x, y, exact = FALSE)$p.value
    data.frame(Definition=defn, Group=g, p=p, stringsAsFactors=FALSE)
  })
  out <- do.call(rbind, out)
  out$q <- p.adjust(out$p, method = "BH")
  out$q_signif <- vapply(out$q, sig_star, character(1))
  out
}))

# summary mean±sd
sum_res <- do.call(rbind, lapply(levels(dctl_long$Definition), function(defn) {
  sub <- dctl_long[dctl_long$Definition == defn & !is.na(dctl_long$Distance), ]
  glev <- levels(sub$Group)
  glev <- glev[!is.na(glev)]
  do.call(rbind, lapply(glev, function(g) {
    v <- sub$Distance[sub$Group == g]
    data.frame(
      Definition=defn,
      Group=g,
      n=sum(!is.na(v)),
      mean=mean(v, na.rm = TRUE),
      sd=sd(v, na.rm = TRUE),
      stringsAsFactors=FALSE
    )
  }))
}))

# merge stats table
stats_tbl <- merge(sum_res, kw_res, by="Definition", all.x=TRUE)
stats_tbl <- merge(stats_tbl, pair_res[, c("Definition","Group","q","q_signif")],
                   by=c("Definition","Group"), all.x=TRUE)
stats_path <- file.path(out_dir, "distance_to_control_stats.csv")
write.csv(stats_tbl, stats_path, row.names = FALSE)

# ---------------- plotting ----------------
# annotate per facet: brackets Control vs each group
make_ann <- function(defn) {
  sub <- dctl_long[dctl_long$Definition == defn & !is.na(dctl_long$Distance), ]
  pr  <- pair_res[pair_res$Definition == defn & !is.na(pair_res$q), ]
  if (nrow(pr) == 0) return(NULL)

  y_max <- max(sub$Distance, na.rm = TRUE)
  y_min <- min(sub$Distance, na.rm = TRUE)
  y_rng <- y_max - y_min
  step <- max(y_rng * 0.12, 0.02)  # adaptive
  base <- y_max + step * 0.35

  # only non-ns if you want hide.ns; here we keep all but you can filter
  pr <- pr[order(pr$Group), ]
  pr$y <- base + (seq_len(nrow(pr)) - 1) * step * 0.35

  data.frame(
    Definition = defn,
    x1 = as.numeric(factor(control_group, levels=levels(meta$Group))),
    x2 = as.numeric(factor(pr$Group, levels=levels(meta$Group))),
    y  = pr$y,
    label = pr$q_signif,
    stringsAsFactors = FALSE
  )
}

ann <- do.call(rbind, lapply(levels(dctl_long$Definition), make_ann))
if (!is.null(ann) && nrow(ann) > 0) {
  ann$Definition <- factor(ann$Definition, levels = levels(dctl_long$Definition))
}

# KW subtitle per facet
kw_map <- setNames(kw_res$kw_p, kw_res$Definition)

p <- ggplot(dctl_long, aes(x = Group, y = Distance)) +
  geom_jitter(width = 0.12, height = 0, size = 2.0, alpha = 0.90, color = "grey35") +
  stat_summary(fun = mean, geom = "point", size = 2.8, color = "black") +
  stat_summary(
    fun.data = function(x) data.frame(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x)),
    geom = "errorbar",
    linewidth = 0.7,
    width = 0.20,
    color = "black"
  ) +
  facet_wrap(~Definition, nrow = 1, scales = "free_y") +
  labs(
    x = NULL,
    y = "Bray–Curtis distance-to-control"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# add KW p in each facet as subtitle-like text
p <- p + geom_text(
  data = data.frame(
    Definition = factor(levels(dctl_long$Definition), levels = levels(dctl_long$Definition)),
    x = 1,
    y = tapply(dctl_long$Distance, dctl_long$Definition, function(v) max(v, na.rm=TRUE)) * 1.08,
    lab = sprintf("Kruskal–Wallis: P = %.3g", kw_map[levels(dctl_long$Definition)])
  ),
  aes(x = x, y = y, label = lab),
  inherit.aes = FALSE,
  hjust = 0,
  size = 3.2
)

# add brackets + stars
if (!is.null(ann) && nrow(ann) > 0) {
  # horizontal line
  p <- p + geom_segment(
    data = ann,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.4
  )
  # vertical ticks
  p <- p + geom_segment(
    data = ann,
    aes(x = x1, xend = x1, y = y, yend = y - 0.01),
    inherit.aes = FALSE,
    linewidth = 0.4
  )
  p <- p + geom_segment(
    data = ann,
    aes(x = x2, xend = x2, y = y, yend = y - 0.01),
    inherit.aes = FALSE,
    linewidth = 0.4
  )
  # labels
  p <- p + geom_text(
    data = ann,
    aes(x = (x1 + x2) / 2, y = y + 0.01, label = label),
    inherit.aes = FALSE,
    size = 3.3
  )
}

pdf_path  <- file.path(out_dir, "FigS_distance_to_control.pdf")
tiff_path <- file.path(out_dir, "FigS_distance_to_control.tiff")

ggsave(pdf_path, p, width = 8.5, height = 4.2, units = "in")
ggsave(tiff_path, p, width = 8.5, height = 4.2, units = "in",
       dpi = 600, device = "tiff", compression = "lzw")

message("Done. Stats: ", stats_path)
message("Done. Figure: ", pdf_path)
