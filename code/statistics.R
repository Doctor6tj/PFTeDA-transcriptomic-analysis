#!/usr/bin/env Rscript

# PFTeDA transcriptomic analysis across three independent study families.
# Family 3 uses five shared doses for the primary summary, with two dose sensitivities.

options(stringsAsFactors = FALSE, warn = 1)

parse_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

config_path <- parse_arg("config")
run_dir <- parse_arg("run-dir")
mode <- parse_arg("mode", "full")
if (is.null(config_path) || is.null(run_dir)) {
  stop("Required arguments: --config=<json> --run-dir=<directory> [--mode=preflight|full]")
}

dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c("metadata", "intermediate", "gene_results", "pathway_results", "qc", "figures", "logs")) {
  dir.create(file.path(run_dir, d), recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(run_dir, "logs", "R_ANALYSIS.log")
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
cfg <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
project_root <- cfg$project_root

required_packages <- c("edgeR", "limma", "singscore", "data.table", "jsonlite", "ggplot2", "matrixStats")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(singscore)
})

set.seed(as.integer(cfg$random_seed))
threshold <- cfg$thresholds
download_root <- file.path(project_root, "data")
geo_root <- file.path(download_root, "geo")
metadata_source <- file.path(project_root, "data", "sample_metadata.csv")
gmt_path <- file.path(download_root, "gene_sets", "MSigDB_Hallmark_2020.gmt")

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
or_else <- function(a, b) if (!is.null(a)) a else b

write_csv <- function(x, path, gzip = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (gzip) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(x, con, row.names = FALSE, na = "")
  } else {
    write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
}

write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    unique(toupper(fields[-c(1, 2)]))
  })
  names(sets) <- vapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]][1], character(1))
  # The retained Enrichr-mirrored Hallmark GMT contains one upstream display-label
  # typo ("Pperoxisome"). Normalize that unique label after parsing while leaving
  # the immutable source snapshot, gene membership and all statistics unchanged.
  names(sets)[names(sets) == "Pperoxisome"] <- "Peroxisome"
  if (anyDuplicated(names(sets))) stop("Pathway-label normalization created duplicate gene-set names")
  sets
}

gene_sets <- read_gmt(gmt_path)
gene_set_sizes <- vapply(gene_sets, length, integer(1))

aggregate_gene_rows <- function(mat) {
  genes <- toupper(sub("_[0-9]+$", "", rownames(mat)))
  rowsum(as.matrix(mat), group = genes, reorder = FALSE)
}

assert_full_rank <- function(design, label) {
  rank <- qr(design)$rank
  if (rank != ncol(design)) stop(label, " design is not full rank: rank=", rank, " columns=", ncol(design))
  invisible(TRUE)
}

contrast_vector <- function(design, weights) {
  out <- setNames(rep(0, ncol(design)), colnames(design))
  missing <- setdiff(names(weights), names(out))
  if (length(missing)) stop("Contrast refers to absent coefficients: ", paste(missing, collapse = ", "))
  out[names(weights)] <- weights
  out
}

add_weight <- function(weights, name, value) {
  if (name %in% names(weights)) weights[[name]] <- weights[[name]] + value else weights[[name]] <- value
  weights
}

dose_key <- function(x) {
  value <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  value <- sub("\\.0+$", "", value)
  gsub("\\.", "p", value)
}

formal_metadata <- function(meta) {
  meta$family_id <- ifelse(
    meta$dataset == "GSE336490", "FAMILY_1",
    ifelse(meta$dataset %in% c("GSE244107", "GSE244109"), "FAMILY_2",
           ifelse(meta$dataset == "GSE145239", "FAMILY_3", "OUT_OF_SCOPE"))
  )
  meta$independent_family_unit <- meta$family_id
  meta$formal_analysis_role <- "OUT_OF_SCOPE"
  meta$formal_analysis_role[meta$dataset == "GSE336490" & meta$model %in% c("HIEC6", "HaCaT")] <- "PRIMARY_DISCOVERY_CONTEXT"
  meta$formal_analysis_role[meta$dataset == "GSE336490" & meta$model == "BEAS2B"] <- "SUPPLEMENTARY_SENSITIVITY_REAL_REPORTED_VALUE_VALIDITY_UNRESOLVED"
  meta$formal_analysis_role[meta$dataset == "GSE336490" & meta$model == "HL7702"] <- "EXPLORATORY_MISIDENTIFIED_CELL_LINE_CONTEXT"
  meta$formal_analysis_role[meta$dataset %in% c("GSE244107", "GSE244109")] <- "PRIMARY_REPLICATION_CONTEXT_WITHIN_ONE_FAMILY"
  meta$formal_analysis_role[meta$dataset == "GSE145239"] <- "PRIMARY_REPLICATION_CONTEXT"
  meta$included_in_family_vote <- meta$formal_analysis_role %in% c(
    "PRIMARY_DISCOVERY_CONTEXT", "PRIMARY_REPLICATION_CONTEXT_WITHIN_ONE_FAMILY", "PRIMARY_REPLICATION_CONTEXT"
  )
  meta
}

build_f2_meta <- function(hash_path, context) {
  h <- read.csv(gzfile(hash_path), check.names = FALSE)
  keep <- h$Chemical_name %in% c("VEH", "Perfluorotetradecanoic acid", "Perfluorooctanoic acid")
  h <- h[keep, , drop = FALSE]
  h$treatment <- c(
    "VEH" = "Vehicle",
    "Perfluorotetradecanoic acid" = "PFTeDA",
    "Perfluorooctanoic acid" = "PFOA"
  )[h$Chemical_name]
  h$dose_uM <- as.numeric(h$Chemical_Concentration_uM)
  h$dose_uM[h$treatment == "Vehicle"] <- 0
  h$context <- context
  h$sample_id <- h$sample_ID
  h
}

build_f2_group <- function(meta) {
  ifelse(meta$treatment == "Vehicle", "VEH", paste0(meta$treatment, "_", vapply(meta$dose_uM, dose_key, character(1))))
}

build_f3_group <- function(meta) {
  ctrl <- sub("^batch[0-9]+:", "", meta$batch)
  ifelse(meta$treatment == "Vehicle", paste0("VEH_", ctrl), paste0("PFTeDA_", vapply(meta$dose_uM, dose_key, character(1))))
}

preflight <- function() {
  log_line("PREFLIGHT_START")
  required <- c(
    metadata_source, gmt_path,
    file.path(geo_root, "GSE336490", "GSE336490_RAW.tar"),
    file.path(geo_root, "GSE244107", "GSE244107_iHep_hash.csv.gz"),
    file.path(geo_root, "GSE244107", "GSE244107_iHep_normcounts.csv.gz"),
    file.path(geo_root, "GSE244109", "GSE244109_iCM_hash.csv.gz"),
    file.path(geo_root, "GSE244109", "GSE244109_iCM_normcounts.csv.gz"),
    file.path(geo_root, "GSE145239", "GSE145239_rawCountsBatch2wID.txt.gz"),
    file.path(geo_root, "GSE145239", "GSE145239_rawCountsBatch4wID.txt.gz")
  )
  if (any(!file.exists(required))) stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = "; "))

  meta <- read.csv(metadata_source, check.names = FALSE)
  # The sample metadata was written with an explicit UTF-8 BOM. Under the
  # batch C locale, R preserves those bytes in the first header name.
  names(meta)[1] <- "dataset"
  selected_metadata <- formal_metadata(meta)
  if (!identical(sort(unique(selected_metadata$family_id[selected_metadata$family_id != "OUT_OF_SCOPE"])), c("FAMILY_1", "FAMILY_2", "FAMILY_3"))) {
    stop("Independent family definition is not exactly FAMILY_1/FAMILY_2/FAMILY_3")
  }
  disc_tab <- table(selected_metadata$model[selected_metadata$dataset == "GSE336490"], selected_metadata$treatment[selected_metadata$dataset == "GSE336490"])
  if (!all(disc_tab == 3)) stop("GSE336490 does not retain exact 3x3 mapping per model")
  if (any(selected_metadata$model == "BEAS2B" & selected_metadata$formal_analysis_role != "SUPPLEMENTARY_SENSITIVITY_REAL_REPORTED_VALUE_VALIDITY_UNRESOLVED")) stop("BEAS role mismatch")
  if (any(selected_metadata$model == "HL7702" & selected_metadata$formal_analysis_role != "EXPLORATORY_MISIDENTIFIED_CELL_LINE_CONTEXT")) stop("HL7702 role mismatch")
  write_csv(selected_metadata, file.path(run_dir, "metadata", "SAMPLE_METADATA.csv"))

  tar_files <- utils::untar(file.path(geo_root, "GSE336490", "GSE336490_RAW.tar"), list = TRUE)
  expected_disc <- paste0(selected_metadata$sample_id[selected_metadata$dataset == "GSE336490"], ".txt.gz")
  if (length(setdiff(expected_disc, basename(tar_files)))) stop("GSE336490 tar is missing selected_metadata sample files")

  f2_design_checks <- list()
  for (spec in list(
    list(accession = "GSE244107", prefix = "iHep", context = "iPSC-Hep"),
    list(accession = "GSE244109", prefix = "iCM", context = "iPSC-CM")
  )) {
    hpath <- file.path(geo_root, spec$accession, paste0(spec$accession, "_", spec$prefix, "_hash.csv.gz"))
    mpath <- file.path(geo_root, spec$accession, paste0(spec$accession, "_", spec$prefix, "_normcounts.csv.gz"))
    h <- build_f2_meta(hpath, spec$context)
    header <- names(read.csv(gzfile(mpath), nrows = 1, check.names = FALSE))
    if (length(setdiff(h$sample_id, header))) stop(spec$accession, " matrix missing selected sample columns")
    h$group <- factor(build_f2_group(h), levels = c("VEH", "PFTeDA_0p1", "PFTeDA_1", "PFTeDA_10", "PFOA_0p1", "PFOA_1", "PFOA_10"))
    h$plate2 <- as.integer(h$Plate == sort(unique(h$Plate))[2])
    design <- model.matrix(~0 + group + plate2, data = h)
    colnames(design) <- sub("^group", "", colnames(design))
    assert_full_rank(design, paste0(spec$accession, " context"))
    f2_design_checks[[spec$accession]] <- list(n = nrow(h), rank = qr(design)$rank, columns = ncol(design))
  }

  f3 <- selected_metadata[selected_metadata$dataset == "GSE145239", , drop = FALSE]
  for (duration in c(24, 240)) {
    m <- f3[as.numeric(f3$duration_h) == duration, , drop = FALSE]
    m$group <- factor(build_f3_group(m))
    design <- model.matrix(~0 + group, data = m)
    colnames(design) <- sub("^group", "", colnames(design))
    assert_full_rank(design, paste0("GSE145239_", duration, "h"))
    batch <- if (duration == 24) 2 else 4
    header <- names(read.delim(gzfile(file.path(geo_root, "GSE145239", paste0("GSE145239_rawCountsBatch", batch, "wID.txt.gz"))), nrows = 1, check.names = FALSE))
    if (length(setdiff(m$sample_id, header))) stop("GSE145239 count matrix missing selected_metadata samples at ", duration, "h")
  }

  packages <- data.frame(
    package = required_packages,
    version = vapply(required_packages, function(p) as.character(utils::packageVersion(p)), character(1)),
    stringsAsFactors = FALSE
  )
  write_csv(packages, file.path(run_dir, "metadata", "R_PACKAGE_VERSIONS.csv"))
  manifest <- data.frame(
    path = substring(gsub("\\\\", "/", required), nchar(gsub("\\\\", "/", project_root)) + 2L),
    bytes = file.info(required)$size,
    modified_time = format(file.info(required)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, file.path(run_dir, "metadata", "FORMAL_INPUT_MANIFEST.csv"))
  report <- list(
    project_id = cfg$project_id, mode = "preflight", verdict = "PASS",
    independent_family_count = 3,
    primary_family1_contexts = c("HIEC6", "HaCaT"),
    beas_role = cfg$beas_role, hl7702_role = cfg$hl7702_role, pfoa_role = cfg$pfoa_role,
    f2_design_checks = f2_design_checks,
    package_versions = as.list(setNames(packages$version, packages$package)),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write_json(report, file.path(run_dir, "PREFLIGHT_REPORT.json"))
  log_line("PREFLIGHT_PASS")
  invisible(list(meta = selected_metadata, report = report))
}

read_discovery_counts <- function(meta) {
  tar_path <- file.path(geo_root, "GSE336490", "GSE336490_RAW.tar")
  extract_dir <- file.path(run_dir, "intermediate", "GSE336490_extracted_counts")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(list.files(extract_dir, pattern = "\\.txt\\.gz$"))) utils::untar(tar_path, exdir = extract_dir)
  vectors <- list()
  genes_ref <- NULL
  for (sample in meta$sample_id) {
    path <- file.path(extract_dir, paste0(sample, ".txt.gz"))
    if (!file.exists(path)) stop("Missing extracted count file: ", path)
    tab <- read.delim(gzfile(path), check.names = FALSE)
    if (ncol(tab) != 2) stop("Unexpected discovery count columns: ", path)
    genes <- as.character(tab[[1]])
    if (is.null(genes_ref)) genes_ref <- genes
    if (!identical(genes_ref, genes)) stop("Discovery count row order differs across samples")
    vectors[[sample]] <- as.numeric(tab[[2]])
  }
  mat <- do.call(cbind, vectors)
  rownames(mat) <- genes_ref
  aggregate_gene_rows(mat)
}

read_ipsc_matrix <- function(accession, prefix) {
  path <- file.path(geo_root, accession, paste0(accession, "_", prefix, "_normcounts.csv.gz"))
  mat <- read.csv(gzfile(path), row.names = 1, check.names = FALSE)
  aggregate_gene_rows(as.matrix(mat))
}

read_liver_counts <- function(batch) {
  path <- file.path(geo_root, "GSE145239", paste0("GSE145239_rawCountsBatch", batch, "wID.txt.gz"))
  mat <- read.delim(gzfile(path), row.names = 1, check.names = FALSE)
  aggregate_gene_rows(as.matrix(mat))
}

prepare_sets <- function(genes) {
  genes <- toupper(genes)
  overlap <- lapply(gene_sets, function(s) intersect(genes, s))
  n_overlap <- vapply(overlap, length, integer(1))
  coverage <- n_overlap / gene_set_sizes[names(n_overlap)]
  eligible <- n_overlap >= as.integer(threshold$min_pathway_tested_genes) &
    coverage >= as.numeric(threshold$min_pathway_coverage_fraction)
  table <- data.frame(
    pathway = names(gene_sets),
    set_size = as.integer(gene_set_sizes[names(gene_sets)]),
    overlap_genes = as.integer(n_overlap[names(gene_sets)]),
    coverage_fraction = as.numeric(coverage[names(gene_sets)]),
    eligible = as.logical(eligible[names(gene_sets)]),
    stringsAsFactors = FALSE
  )
  list(index = overlap[eligible], table = table)
}

score_gene_sets <- function(expr, set_index) {
  ranked <- singscore::rankGenes(expr)
  scores <- matrix(NA_real_, nrow = length(set_index), ncol = ncol(expr),
                   dimnames = list(names(set_index), colnames(expr)))
  for (i in seq_along(set_index)) {
    scored <- singscore::simpleScore(ranked, upSet = set_index[[i]], centerScore = TRUE, knownDirection = TRUE)
    values <- if ("TotalScore" %in% colnames(scored)) scored$TotalScore else scored[[1]]
    names(values) <- rownames(scored)
    scores[i, names(values)] <- as.numeric(values)
  }
  scores
}

pathway_contrast <- function(camera_y, expr, design, contrast, context_info, contrast_name, contrast_type, score_matrix = NULL) {
  prep <- prepare_sets(rownames(expr))
  eligible_index <- prep$index
  base <- prep$table
  if (!length(eligible_index)) stop("No eligible Hallmark pathways for ", context_info$context)
  cam <- limma::camera(
    camera_y, index = eligible_index, design = design, contrast = contrast,
    inter.gene.cor = as.numeric(threshold$camera_inter_gene_correlation), sort = FALSE
  )
  cam$pathway <- rownames(cam)
  rownames(cam) <- NULL
  names(cam)[names(cam) == "NGenes"] <- "camera_n_genes"
  names(cam)[names(cam) == "Direction"] <- "camera_direction"
  names(cam)[names(cam) == "PValue"] <- "camera_p_value"
  names(cam)[names(cam) == "FDR"] <- "camera_fdr"

  if (is.null(score_matrix)) score_matrix <- score_gene_sets(expr, eligible_index)
  score_fit <- limma::lmFit(score_matrix, design)
  score_contrast <- limma::contrasts.fit(score_fit, contrasts = contrast)
  score_contrast <- limma::eBayes(score_contrast, robust = TRUE)
  score_tab <- limma::topTable(score_contrast, coef = 1, number = Inf, sort.by = "none")
  score_tab$pathway <- rownames(score_tab)
  sigma <- score_contrast$sigma[match(score_tab$pathway, rownames(score_contrast$coefficients))]
  score_out <- data.frame(
    pathway = score_tab$pathway,
    score_effect = score_tab$logFC,
    score_standardized_effect = score_tab$logFC / sigma,
    score_t = score_tab$t,
    score_p_value = score_tab$P.Value,
    score_fdr = score_tab$adj.P.Val,
    score_direction = ifelse(score_tab$logFC >= 0, "UP", "DOWN"),
    stringsAsFactors = FALSE
  )
  out <- merge(base, cam[, c("pathway", "camera_n_genes", "camera_direction", "camera_p_value", "camera_fdr")],
               by = "pathway", all.x = TRUE, sort = FALSE)
  out <- merge(out, score_out, by = "pathway", all.x = TRUE, sort = FALSE)
  out$family_id <- context_info$family_id
  out$dataset <- context_info$dataset
  out$context <- context_info$context
  out$role <- context_info$role
  out$contrast <- contrast_name
  out$contrast_type <- contrast_type
  out$camera_direction <- toupper(out$camera_direction)
  out[, c("family_id", "dataset", "context", "role", "contrast", "contrast_type",
          "pathway", "set_size", "overlap_genes", "coverage_fraction", "eligible",
          "camera_n_genes", "camera_direction", "camera_p_value", "camera_fdr",
          "score_effect", "score_standardized_effect", "score_t", "score_p_value",
          "score_fdr", "score_direction")]
}

plot_pca <- function(expr, meta, file, title) {
  if (ncol(expr) < 3) return(invisible(NULL))
  v <- apply(expr, 1, stats::var)
  use <- order(v, decreasing = TRUE)[seq_len(min(1000, length(v)))]
  p <- stats::prcomp(t(expr[use, , drop = FALSE]), center = TRUE, scale. = FALSE)
  label_group <- or_else(meta$treatment, meta$group)
  group_factor <- factor(label_group)
  grDevices::png(file, width = 1400, height = 1100, res = 160)
  plot(p$x[, 1], p$x[, 2], pch = 19, col = as.integer(group_factor),
       xlab = "PC1", ylab = "PC2", main = title)
  text(p$x[, 1], p$x[, 2], labels = meta$sample_id, pos = 3, cex = 0.55)
  legend("topright", legend = levels(group_factor), col = seq_along(levels(group_factor)), pch = 19, cex = 0.75)
  grDevices::dev.off()
}

run_raw_model <- function(counts, meta, design, contrast_specs, context_info, model_id) {
  stopifnot(identical(colnames(counts), meta$sample_id))
  rownames(design) <- meta$sample_id
  assert_full_rank(design, model_id)
  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design)
  if (sum(keep) < 100) stop(model_id, " has too few genes after filterByExpr")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  expr <- edgeR::cpm(dge, log = TRUE, prior.count = 0.5)
  prep <- prepare_sets(rownames(expr))
  score_matrix <- score_gene_sets(expr, prep$index)
  path_rows <- list()
  gene_files <- character()
  for (name in names(contrast_specs)) {
    spec <- contrast_specs[[name]]
    test <- edgeR::glmTreat(fit, contrast = spec$vector, lfc = as.numeric(threshold$gene_abs_log2fc))
    tab <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
    tab$gene <- rownames(tab)
    rownames(tab) <- NULL
    tab$family_id <- context_info$family_id
    tab$dataset <- context_info$dataset
    tab$context <- context_info$context
    tab$role <- context_info$role
    tab$contrast <- name
    tab$contrast_type <- spec$type
    tab$test_statistic <- if ("F" %in% names(tab)) {
      tab$F
    } else if ("LR" %in% names(tab)) {
      tab$LR
    } else {
      rep(NA_real_, nrow(tab))
    }
    tab$gene_call <- tab$FDR <= as.numeric(threshold$gene_fdr) &
      abs(tab$logFC) >= as.numeric(threshold$gene_abs_log2fc)
    tab <- tab[, c("family_id", "dataset", "context", "role", "contrast", "contrast_type",
                   "gene", "logFC", "logCPM", "test_statistic", "PValue", "FDR", "gene_call")]
    gpath <- file.path(run_dir, "gene_results", paste0(safe_name(model_id), "__", safe_name(name), ".csv.gz"))
    write_csv(tab, gpath, gzip = TRUE)
    gene_files <- c(gene_files, gpath)
    path_rows[[name]] <- pathway_contrast(dge, expr, design, spec$vector, context_info, name, spec$type, score_matrix)
  }
  path_tab <- do.call(rbind, path_rows)
  write_csv(path_tab, file.path(run_dir, "pathway_results", paste0(safe_name(model_id), "__PATHWAYS.csv")))
  saveRDS(
    list(dge = dge, design = design, fit = fit, contrasts = contrast_specs,
         score_matrix = score_matrix, metadata = meta, context_info = context_info),
    file.path(run_dir, "intermediate", paste0(safe_name(model_id), "__MODEL_OBJECTS.rds")), compress = "xz"
  )
  plot_pca(expr, meta, file.path(run_dir, "figures", paste0(safe_name(model_id), "__PCA.png")), model_id)
  qc <- data.frame(
    model_id = model_id, family_id = context_info$family_id, dataset = context_info$dataset,
    context = context_info$context, samples = ncol(counts), genes_input = nrow(counts),
    genes_tested = nrow(dge), design_columns = ncol(design), design_rank = qr(design)$rank,
    residual_df_min = min(fit$df.residual), normalization = "TMM_edgeR_QL", stringsAsFactors = FALSE
  )
  list(pathways = path_tab, qc = qc,
       model = list(fit = fit, design = design, expr = expr, score_matrix = score_matrix,
                    context_info = context_info, model_id = model_id),
       gene_files = gene_files)
}

run_log_model <- function(expr, meta, design, contrast_specs, context_info, model_id) {
  stopifnot(identical(colnames(expr), meta$sample_id))
  rownames(design) <- meta$sample_id
  assert_full_rank(design, model_id)
  keep <- apply(expr, 1, stats::var) > 0 & rowMeans(2^expr - 1) >= 1
  expr <- expr[keep, , drop = FALSE]
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  prep <- prepare_sets(rownames(expr))
  score_matrix <- score_gene_sets(expr, prep$index)
  path_rows <- list()
  gene_files <- character()
  for (name in names(contrast_specs)) {
    spec <- contrast_specs[[name]]
    cfit <- limma::contrasts.fit(limma::lmFit(expr, design), contrasts = spec$vector)
    cfit <- limma::treat(cfit, lfc = as.numeric(threshold$gene_abs_log2fc), trend = TRUE, robust = TRUE)
    tab <- limma::topTreat(cfit, coef = 1, number = Inf, sort.by = "none")
    tab$gene <- rownames(tab)
    rownames(tab) <- NULL
    tab$family_id <- context_info$family_id
    tab$dataset <- context_info$dataset
    tab$context <- context_info$context
    tab$role <- context_info$role
    tab$contrast <- name
    tab$contrast_type <- spec$type
    tab$gene_call <- tab$adj.P.Val <= as.numeric(threshold$gene_fdr) &
      abs(tab$logFC) >= as.numeric(threshold$gene_abs_log2fc)
    tab <- tab[, c("family_id", "dataset", "context", "role", "contrast", "contrast_type",
                   "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "gene_call")]
    names(tab)[names(tab) == "P.Value"] <- "PValue"
    names(tab)[names(tab) == "adj.P.Val"] <- "FDR"
    gpath <- file.path(run_dir, "gene_results", paste0(safe_name(model_id), "__", safe_name(name), ".csv.gz"))
    write_csv(tab, gpath, gzip = TRUE)
    gene_files <- c(gene_files, gpath)
    path_rows[[name]] <- pathway_contrast(expr, expr, design, spec$vector, context_info, name, spec$type, score_matrix)
  }
  path_tab <- do.call(rbind, path_rows)
  write_csv(path_tab, file.path(run_dir, "pathway_results", paste0(safe_name(model_id), "__PATHWAYS.csv")))
  saveRDS(
    list(expr = expr, design = design, fit = fit, contrasts = contrast_specs,
         score_matrix = score_matrix, metadata = meta, context_info = context_info),
    file.path(run_dir, "intermediate", paste0(safe_name(model_id), "__MODEL_OBJECTS.rds")), compress = "xz"
  )
  plot_pca(expr, meta, file.path(run_dir, "figures", paste0(safe_name(model_id), "__PCA.png")), model_id)
  qc <- data.frame(
    model_id = model_id, family_id = context_info$family_id, dataset = context_info$dataset,
    context = context_info$context, samples = ncol(expr), genes_input = nrow(expr),
    genes_tested = nrow(expr), design_columns = ncol(design), design_rank = qr(design)$rank,
    residual_df_min = min(fit$df.residual),
    normalization = "repository_normalized_counts_log2_plus1_limma_trend", stringsAsFactors = FALSE
  )
  list(pathways = path_tab, qc = qc,
       model = list(fit = fit, design = design, expr = expr, score_matrix = score_matrix,
                    context_info = context_info, model_id = model_id),
       gene_files = gene_files)
}

write_multidf_omnibus <- function(model, contrast_matrix, name, context_info) {
  fit <- limma::lmFit(model$expr, model$design)
  cfit <- limma::contrasts.fit(fit, contrasts = contrast_matrix)
  cfit <- limma::eBayes(cfit, trend = TRUE, robust = TRUE)
  gene <- limma::topTable(cfit, coef = seq_len(ncol(contrast_matrix)), number = Inf, sort.by = "none")
  gene$gene <- rownames(gene)
  gene$family_id <- context_info$family_id
  gene$dataset <- context_info$dataset
  gene$context <- context_info$context
  gene$contrast <- name
  write_csv(gene, file.path(run_dir, "gene_results", paste0(safe_name(model$model_id), "__", safe_name(name), "__OMNIBUS.csv.gz")), gzip = TRUE)
  sfit <- limma::lmFit(model$score_matrix, model$design)
  sc <- limma::contrasts.fit(sfit, contrasts = contrast_matrix)
  sc <- limma::eBayes(sc, robust = TRUE)
  score <- limma::topTable(sc, coef = seq_len(ncol(contrast_matrix)), number = Inf, sort.by = "none")
  score$pathway <- rownames(score)
  score$family_id <- context_info$family_id
  score$dataset <- context_info$dataset
  score$context <- context_info$context
  score$contrast <- name
  write_csv(score, file.path(run_dir, "pathway_results", paste0(safe_name(model$model_id), "__", safe_name(name), "__SCORE_OMNIBUS.csv")))
}

f2_contrast_specs <- function(design) {
  specs <- list()
  for (chem in c("PFTeDA", "PFOA")) {
    for (dose in c("0p1", "1", "10")) {
      w <- c(1, -1); names(w) <- c(paste0(chem, "_", dose), "VEH")
      specs[[paste0(chem, "_", dose, "_vs_VEH")]] <- list(vector = contrast_vector(design, w), type = paste0(chem, "_DOSE_RESPONSE"))
    }
    w <- c(1/3, 1/3, 1/3, -1)
    names(w) <- c(paste0(chem, "_0p1"), paste0(chem, "_1"), paste0(chem, "_10"), "VEH")
    ctype <- if (chem == "PFTeDA") "PFTeDA_CONTEXT_SUMMARY" else "PFOA_CONTEXT_RESPONSE"
    specs[[paste0(chem, "_summary_vs_VEH")]] <- list(vector = contrast_vector(design, w), type = ctype)
  }
  for (dose in c("0p1", "1", "10")) {
    w <- c(1, -1); names(w) <- c(paste0("PFTeDA_", dose), paste0("PFOA_", dose))
    specs[[paste0("PFTeDA_minus_PFOA_", dose)]] <- list(vector = contrast_vector(design, w), type = "DIRECT_COMPARATOR")
  }
  w <- c(1/3, 1/3, 1/3, -1/3, -1/3, -1/3)
  names(w) <- c("PFTeDA_0p1", "PFTeDA_1", "PFTeDA_10", "PFOA_0p1", "PFOA_1", "PFOA_10")
  specs[["PFTeDA_minus_PFOA_mean_shared_doses"]] <- list(vector = contrast_vector(design, w), type = "DIRECT_COMPARATOR")
  w1 <- c(1, -1, -1, 1); names(w1) <- c("PFTeDA_1", "PFOA_1", "PFTeDA_0p1", "PFOA_0p1")
  w2 <- c(1, -1, -1, 1); names(w2) <- c("PFTeDA_10", "PFOA_10", "PFTeDA_0p1", "PFOA_0p1")
  specs[["chemical_by_dose_1_vs_0p1"]] <- list(vector = contrast_vector(design, w1), type = "CHEMICAL_DOSE_INTERACTION_BASIS")
  specs[["chemical_by_dose_10_vs_0p1"]] <- list(vector = contrast_vector(design, w2), type = "CHEMICAL_DOSE_INTERACTION_BASIS")
  specs
}

f3_contrast_specs <- function(meta, design, true_common_doses, lowest_four_shared_doses) {
  doses <- sort(unique(meta$dose_uM[meta$treatment == "PFTeDA"]))
  specs <- list()
  all_w <- numeric()
  common_w <- numeric()
  lowest_four_w <- numeric()
  for (dose in doses) {
    key <- dose_key(dose)
    ctrl <- unique(sub("^batch[0-9]+:", "", meta$batch[meta$treatment == "PFTeDA" & meta$dose_uM == dose]))
    if (length(ctrl) != 1) stop("PFTeDA dose does not map to one control type: ", dose)
    w <- c(1, -1); names(w) <- c(paste0("PFTeDA_", key), paste0("VEH_", ctrl))
    specs[[paste0("PFTeDA_", key, "_vs_matched_VEH")]] <- list(vector = contrast_vector(design, w), type = "PFTeDA_DOSE_RESPONSE")
    all_w <- add_weight(all_w, paste0("PFTeDA_", key), 1/length(doses))
    all_w <- add_weight(all_w, paste0("VEH_", ctrl), -1/length(doses))
    if (dose %in% true_common_doses) {
      common_w <- add_weight(common_w, paste0("PFTeDA_", key), 1/length(true_common_doses))
      common_w <- add_weight(common_w, paste0("VEH_", ctrl), -1/length(true_common_doses))
    }
    if (dose %in% lowest_four_shared_doses) {
      lowest_four_w <- add_weight(lowest_four_w, paste0("PFTeDA_", key), 1/length(lowest_four_shared_doses))
      lowest_four_w <- add_weight(lowest_four_w, paste0("VEH_", ctrl), -1/length(lowest_four_shared_doses))
    }
  }
  specs[["PFTeDA_common_five_dose_summary_vs_matched_VEH"]] <- list(
    vector = contrast_vector(design, common_w), type = "PFTeDA_CONTEXT_SUMMARY"
  )
  specs[["PFTeDA_all_available_dose_summary_vs_matched_VEH"]] <- list(
    vector = contrast_vector(design, all_w), type = "PFTeDA_F3_ALL_AVAILABLE_DOSE_SUMMARY"
  )
  specs[["PFTeDA_lowest_four_shared_dose_summary_vs_matched_VEH"]] <- list(
    vector = contrast_vector(design, lowest_four_w), type = "PFTeDA_F3_LOWEST_FOUR_SHARED_DOSE_SUMMARY"
  )
  specs
}

run_full <- function(selected_metadata) {
  log_line("FORMAL_STATISTICS_START")
  all_pathways <- list()
  all_qc <- list()
  context_manifest <- list()

  disc_meta <- selected_metadata[selected_metadata$dataset == "GSE336490", , drop = FALSE]
  disc_counts <- read_discovery_counts(disc_meta)
  for (context in c("HIEC6", "HaCaT", "BEAS2B", "HL7702")) {
    m <- disc_meta[disc_meta$model == context, , drop = FALSE]
    m$group <- factor(m$treatment, levels = c("Vehicle", "PFOA", "PFTeDA"))
    design <- model.matrix(~0 + group, data = m)
    colnames(design) <- sub("^group", "", colnames(design))
    counts <- disc_counts[, m$sample_id, drop = FALSE]
    primary <- context %in% c("HIEC6", "HaCaT")
    role <- if (primary) "PRIMARY" else if (context == "BEAS2B") "SUPPLEMENTARY_SENSITIVITY" else "EXPLORATORY_MISIDENTIFIED_CONTEXT"
    info <- list(family_id = "FAMILY_1", dataset = "GSE336490", context = context, role = role)
    specs <- list(
      PFTeDA_vs_Vehicle = list(vector = contrast_vector(design, c(PFTeDA = 1, Vehicle = -1)),
                               type = if (primary) "PFTeDA_CONTEXT_SUMMARY" else "SUPPLEMENTARY_PFTeDA_RESPONSE"),
      PFOA_vs_Vehicle = list(vector = contrast_vector(design, c(PFOA = 1, Vehicle = -1)),
                             type = if (primary) "PFOA_CONTEXT_RESPONSE" else "SUPPLEMENTARY_PFOA_RESPONSE"),
      PFTeDA_minus_PFOA = list(vector = contrast_vector(design, c(PFTeDA = 1, PFOA = -1)),
                               type = if (primary) "DIRECT_COMPARATOR_EFFECT_PATTERN_ONLY" else "SUPPLEMENTARY_DIRECT_COMPARATOR")
    )
    result <- run_raw_model(counts, m, design, specs, info, paste0("F1_", context))
    all_pathways[[paste0("F1_", context)]] <- result$pathways
    all_qc[[paste0("F1_", context)]] <- result$qc
    context_manifest[[length(context_manifest) + 1]] <- data.frame(
      family_id = "FAMILY_1", dataset = "GSE336490", context = context, role = role,
      independent_family_unit = "FAMILY_1", stringsAsFactors = FALSE
    )
  }

  m <- disc_meta[disc_meta$model %in% c("HIEC6", "HaCaT"), , drop = FALSE]
  m$group <- factor(paste(m$model, m$treatment, sep = "_"))
  design <- model.matrix(~0 + group, data = m)
  colnames(design) <- sub("^group", "", colnames(design))
  w <- c(1, -1, -1, 1)
  names(w) <- c("HIEC6_PFTeDA", "HIEC6_Vehicle", "HaCaT_PFTeDA", "HaCaT_Vehicle")
  specs <- list(F1_HIEC6_minus_HaCaT_PFTeDA_response = list(vector = contrast_vector(design, w), type = "HETEROGENEITY"))
  info <- list(family_id = "FAMILY_1", dataset = "GSE336490",
               context = "HIEC6_vs_HaCaT_study_context", role = "FORMAL_HETEROGENEITY")
  result <- run_raw_model(disc_counts[, m$sample_id, drop = FALSE], m, design, specs, info, "F1_PRIMARY_CONTEXT_INTERACTION")
  all_pathways[["F1_INTERACTION"]] <- result$pathways
  all_qc[["F1_INTERACTION"]] <- result$qc

  f2_data <- list()
  for (spec0 in list(
    list(accession = "GSE244107", prefix = "iHep", context = "iPSC-Hep"),
    list(accession = "GSE244109", prefix = "iCM", context = "iPSC-CM")
  )) {
    mat <- read_ipsc_matrix(spec0$accession, spec0$prefix)
    h <- build_f2_meta(file.path(geo_root, spec0$accession, paste0(spec0$accession, "_", spec0$prefix, "_hash.csv.gz")), spec0$context)
    if (length(setdiff(h$sample_id, colnames(mat)))) stop(spec0$accession, " missing matrix samples")
    expr <- log2(mat[, h$sample_id, drop = FALSE] + 1)
    h$group <- factor(build_f2_group(h), levels = c("VEH", "PFTeDA_0p1", "PFTeDA_1", "PFTeDA_10", "PFOA_0p1", "PFOA_1", "PFOA_10"))
    h$plate2 <- as.integer(h$Plate == sort(unique(h$Plate))[2])
    design <- model.matrix(~0 + group + plate2, data = h)
    colnames(design) <- sub("^group", "", colnames(design))
    specs <- f2_contrast_specs(design)
    info <- list(family_id = "FAMILY_2", dataset = spec0$accession, context = spec0$context, role = "PRIMARY")
    result <- run_log_model(expr, h, design, specs, info, paste0("F2_", spec0$context))
    all_pathways[[paste0("F2_", spec0$context)]] <- result$pathways
    all_qc[[paste0("F2_", spec0$context)]] <- result$qc
    f2_data[[spec0$context]] <- list(expr = expr, meta = h, accession = spec0$accession)
    intmat <- cbind(specs[["chemical_by_dose_1_vs_0p1"]]$vector, specs[["chemical_by_dose_10_vs_0p1"]]$vector)
    colnames(intmat) <- c("interaction_1_vs_0p1", "interaction_10_vs_0p1")
    write_multidf_omnibus(result$model, intmat, "chemical_by_dose_omnibus", info)
    context_manifest[[length(context_manifest) + 1]] <- data.frame(
      family_id = "FAMILY_2", dataset = spec0$accession, context = spec0$context,
      role = "PRIMARY", independent_family_unit = "FAMILY_2", stringsAsFactors = FALSE
    )
  }

  common_genes <- Reduce(intersect, lapply(f2_data, function(x) rownames(x$expr)))
  expr_combined <- cbind(f2_data[["iPSC-Hep"]]$expr[common_genes, , drop = FALSE],
                         f2_data[["iPSC-CM"]]$expr[common_genes, , drop = FALSE])
  mh <- f2_data[["iPSC-Hep"]]$meta
  mc <- f2_data[["iPSC-CM"]]$meta
  mh$context_prefix <- "Hep"
  mc$context_prefix <- "CM"
  m <- rbind(mh, mc)
  m$group <- factor(paste0(m$context_prefix, "_", build_f2_group(m)))
  m$cm_plate2 <- as.integer(m$context_prefix == "CM" & m$Plate == sort(unique(mc$Plate))[2])
  m$hep_plate2 <- as.integer(m$context_prefix == "Hep" & m$Plate == sort(unique(mh$Plate))[2])
  design <- model.matrix(~0 + group + cm_plate2 + hep_plate2, data = m)
  colnames(design) <- sub("^group", "", colnames(design))
  w <- numeric()
  for (dose in c("0p1", "1", "10")) {
    w <- add_weight(w, paste0("Hep_PFTeDA_", dose), 1/3)
    w <- add_weight(w, paste0("CM_PFTeDA_", dose), -1/3)
  }
  w <- add_weight(w, "Hep_VEH", -1)
  w <- add_weight(w, "CM_VEH", 1)
  specs <- list(F2_Hep_minus_CM_PFTeDA_summary_response = list(vector = contrast_vector(design, w), type = "HETEROGENEITY"))
  info <- list(family_id = "FAMILY_2", dataset = "GSE244107+GSE244109",
               context = "iPSC_Hep_vs_CM_protocol_context", role = "FORMAL_HETEROGENEITY")
  result <- run_log_model(expr_combined[, m$sample_id, drop = FALSE], m, design, specs, info, "F2_CONTEXT_INTERACTION")
  all_pathways[["F2_INTERACTION"]] <- result$pathways
  all_qc[["F2_INTERACTION"]] <- result$qc

  f3_meta <- selected_metadata[selected_metadata$dataset == "GSE145239", , drop = FALSE]
  f3_doses_by_duration <- split(
    f3_meta$dose_uM[f3_meta$treatment == "PFTeDA"],
    f3_meta$duration_h[f3_meta$treatment == "PFTeDA"]
  )
  f3_doses_by_duration <- lapply(f3_doses_by_duration, function(x) sort(unique(as.numeric(x))))
  f3_true_common_doses <- sort(Reduce(intersect, f3_doses_by_duration))
  f3_expected_common_doses <- c(0.06, 0.67, 3.35, 6.7, 17)
  if (!identical(f3_true_common_doses, f3_expected_common_doses)) {
    stop(
      "Unexpected Family 3 dose intersection: ",
      paste(f3_true_common_doses, collapse = ", "),
      "; expected ", paste(f3_expected_common_doses, collapse = ", ")
    )
  }
  f3_lowest_four_shared_doses <- f3_true_common_doses[seq_len(4)]
  f3_dose_sets <- do.call(rbind, lapply(names(f3_doses_by_duration), function(duration) {
    data.frame(
      duration_h = as.numeric(duration),
      dose_uM = f3_doses_by_duration[[duration]],
      in_true_common_five = f3_doses_by_duration[[duration]] %in% f3_true_common_doses,
      in_lowest_four_sensitivity = f3_doses_by_duration[[duration]] %in% f3_lowest_four_shared_doses,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(f3_dose_sets, file.path(run_dir, "metadata", "FAMILY3_DOSE_SETS.csv"))
  write_json(
    list(
      primary_summary = "TRUE_COMMON_FIVE_DOSES",
      primary_doses_uM = f3_true_common_doses,
      sensitivity_all_available = lapply(f3_doses_by_duration, as.numeric),
      sensitivity_lowest_four_shared_uM = f3_lowest_four_shared_doses,
      dose_unit = "uM",
      duration_contexts_h = c(24, 240)
    ),
    file.path(run_dir, "metadata", "FAMILY3_DOSE_DEFINITION.json")
  )
  f3_data <- list()
  for (duration in c(24, 240)) {
    batch <- if (duration == 24) 2 else 4
    counts_all <- read_liver_counts(batch)
    m <- f3_meta[as.numeric(f3_meta$duration_h) == duration, , drop = FALSE]
    m$group <- factor(build_f3_group(m))
    design <- model.matrix(~0 + group, data = m)
    colnames(design) <- sub("^group", "", colnames(design))
    counts <- counts_all[, m$sample_id, drop = FALSE]
    specs <- f3_contrast_specs(m, design, f3_true_common_doses, f3_lowest_four_shared_doses)
    info <- list(family_id = "FAMILY_3", dataset = "GSE145239",
                 context = paste0("liver_spheroid_", duration, "h"), role = "PRIMARY")
    result <- run_raw_model(counts, m, design, specs, info, paste0("F3_liver_", duration, "h"))
    all_pathways[[paste0("F3_", duration)]] <- result$pathways
    all_qc[[paste0("F3_", duration)]] <- result$qc
    f3_data[[as.character(duration)]] <- list(counts = counts, meta = m, specs = specs)
    context_manifest[[length(context_manifest) + 1]] <- data.frame(
      family_id = "FAMILY_3", dataset = "GSE145239", context = paste0("liver_spheroid_", duration, "h"),
      role = "PRIMARY", independent_family_unit = "FAMILY_3", stringsAsFactors = FALSE
    )
  }

  common_genes <- intersect(rownames(f3_data[["24"]]$counts), rownames(f3_data[["240"]]$counts))
  counts_combined <- cbind(f3_data[["24"]]$counts[common_genes, , drop = FALSE],
                           f3_data[["240"]]$counts[common_genes, , drop = FALSE])
  m24 <- f3_data[["24"]]$meta
  m240 <- f3_data[["240"]]$meta
  m24$duration_prefix <- "D24"
  m240$duration_prefix <- "D240"
  m <- rbind(m24, m240)
  m$group <- factor(paste0(m$duration_prefix, "_", build_f3_group(m)))
  design <- model.matrix(~0 + group, data = m)
  colnames(design) <- sub("^group", "", colnames(design))
  f3_interaction_vector <- function(dose_set) {
    w <- numeric()
    for (duration_prefix in c("D24", "D240")) {
      sign <- if (duration_prefix == "D24") 1 else -1
      source_meta <- if (duration_prefix == "D24") m24 else m240
      for (dose in dose_set) {
        ctrl <- unique(sub("^batch[0-9]+:", "", source_meta$batch[source_meta$treatment == "PFTeDA" & source_meta$dose_uM == dose]))
        key <- dose_key(dose)
        w <- add_weight(w, paste0(duration_prefix, "_PFTeDA_", key), sign/length(dose_set))
        w <- add_weight(w, paste0(duration_prefix, "_VEH_", ctrl), -sign/length(dose_set))
      }
    }
    contrast_vector(design, w)
  }
  specs <- list(
    F3_24h_minus_240h_true_common_five_dose_response = list(
      vector = f3_interaction_vector(f3_true_common_doses), type = "HETEROGENEITY"
    ),
    F3_24h_minus_240h_lowest_four_shared_dose_response = list(
      vector = f3_interaction_vector(f3_lowest_four_shared_doses),
      type = "HETEROGENEITY_SENSITIVITY_FOUR_DOSES"
    )
  )
  info <- list(family_id = "FAMILY_3", dataset = "GSE145239",
               context = "liver_24h_vs_240h_duration_batch", role = "FORMAL_HETEROGENEITY")
  result <- run_raw_model(counts_combined[, m$sample_id, drop = FALSE], m, design, specs, info, "F3_DURATION_BATCH_INTERACTION")
  all_pathways[["F3_INTERACTION"]] <- result$pathways
  all_qc[["F3_INTERACTION"]] <- result$qc

  pathway_all <- do.call(rbind, all_pathways)
  rownames(pathway_all) <- NULL
  write_csv(pathway_all, file.path(run_dir, "pathway_results", "ALL_PATHWAY_RESULTS.csv"))
  saveRDS(pathway_all, file.path(run_dir, "intermediate", "ALL_PATHWAY_RESULTS.rds"), compress = "xz")
  qc_all <- do.call(rbind, all_qc)
  rownames(qc_all) <- NULL
  write_csv(qc_all, file.path(run_dir, "qc", "MODEL_QC_SUMMARY.csv"))
  contexts <- do.call(rbind, context_manifest)
  write_csv(contexts, file.path(run_dir, "metadata", "CONTEXT_FAMILY_MANIFEST.csv"))
  capture.output(sessionInfo(), file = file.path(run_dir, "metadata", "R_SESSION_INFO.txt"))
  summary <- list(
    project_id = cfg$project_id, status = "STATISTICS_COMPLETE",
    independent_family_count = 3, model_count = nrow(qc_all), pathway_result_rows = nrow(pathway_all),
    gene_result_file_count = length(list.files(file.path(run_dir, "gene_results"), pattern = "\\.csv\\.gz$")),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write_json(summary, file.path(run_dir, "STATISTICAL_RUN_SUMMARY.json"))
  log_line("FORMAL_STATISTICS_COMPLETE", "models=", nrow(qc_all), "pathway_rows=", nrow(pathway_all))
}

main <- function() {
  pf <- preflight()
  if (mode == "preflight") return(invisible(NULL))
  if (mode != "full") stop("Unknown mode: ", mode)
  run_full(pf$meta)
}

tryCatch(
  {
    main()
    write_json(list(project_id = cfg$project_id, mode = mode, status = "PASS",
                    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
               file.path(run_dir, "RUN_STATUS.json"))
  },
  error = function(e) {
    log_line("RUN_ERROR", conditionMessage(e))
    write_json(list(project_id = cfg$project_id, mode = mode, status = "FAIL",
                    error = conditionMessage(e),
                    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
               file.path(run_dir, "RUN_STATUS.json"))
    stop(e)
  }
)
