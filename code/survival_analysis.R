#!/usr/bin/env Rscript
# Gene expression survival analysis: Kaplan-Meier + log-rank + Cox.
# Base R + survival + ggplot2 only. Paths honor DATA_DIR (default /data) and RESULTS_DIR (default /results).

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

# ---------------------------------------------------------------- parameters
defaults <- list(
  "gene"       = "TP53",
  "split"      = "median",
  "covariates" = "age,stage",
  "time-col"   = "os_months",
  "event-col"  = "os_event"
)

parse_args <- function(argv, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("-h", "--help")) {
      cat("Usage: code/run [--gene TP53[,EGFR,...]] [--split median|quartile]\n",
          "                [--covariates age,stage|none] [--time-col os_months] [--event-col os_event]\n", sep = "")
      quit(status = 0)
    }
    if (!startsWith(a, "--")) stop("Unexpected argument '", a, "'. Flags look like --gene TP53.", call. = FALSE)
    key <- sub("^--", "", a)
    val <- NULL
    if (grepl("=", key, fixed = TRUE)) {
      val <- sub("^[^=]*=", "", key)
      key <- sub("=.*$", "", key)
    } else {
      if (i == length(argv)) stop("Flag --", key, " needs a value.", call. = FALSE)
      i <- i + 1
      val <- argv[i]
    }
    if (!key %in% names(defaults)) {
      stop("Unknown flag --", key, ". Valid flags: ", paste0("--", names(defaults), collapse = ", "), call. = FALSE)
    }
    opts[[key]] <- val
    i <- i + 1
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)

DATA_DIR    <- Sys.getenv("DATA_DIR", "/data")
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", "/results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

split_mode <- tolower(trimws(opts[["split"]]))
if (!split_mode %in% c("median", "quartile")) {
  stop("--split must be 'median' or 'quartile' (got '", opts[["split"]], "').", call. = FALSE)
}
time_col  <- trimws(opts[["time-col"]])
event_col <- trimws(opts[["event-col"]])
cov_raw   <- trimws(opts[["covariates"]])
covariates <- if (tolower(cov_raw) %in% c("none", "")) character(0) else trimws(strsplit(cov_raw, ",", fixed = TRUE)[[1]])
covariates <- covariates[nzchar(covariates)]
gene_queries <- unique(trimws(strsplit(opts[["gene"]], ",", fixed = TRUE)[[1]]))
gene_queries <- gene_queries[nzchar(gene_queries)]
if (length(gene_queries) == 0) stop("--gene is empty. Example: --gene TP53 or --gene TP53,EGFR", call. = FALSE)

cat("== Gene survival analysis (Kaplan-Meier + Cox) ==\n")
cat("Parameters:\n")
cat("  --gene        ", paste(gene_queries, collapse = ","), "\n", sep = "")
cat("  --split       ", split_mode, "\n", sep = "")
cat("  --covariates  ", if (length(covariates)) paste(covariates, collapse = ",") else "none", "\n", sep = "")
cat("  --time-col    ", time_col, "\n", sep = "")
cat("  --event-col   ", event_col, "\n", sep = "")
cat("  DATA_DIR      ", DATA_DIR, "\n", sep = "")
cat("  RESULTS_DIR   ", RESULTS_DIR, "\n\n", sep = "")

# ---------------------------------------------------------------- inputs
cohort_dir <- file.path(DATA_DIR, "synthetic-tumor-cohort")
annot_dir  <- file.path(DATA_DIR, "synthetic-gene-annotation")
paths <- c(expression = file.path(cohort_dir, "expression.tsv"),
           clinical   = file.path(cohort_dir, "clinical.tsv"),
           annotation = file.path(annot_dir, "genes.tsv"))
missing_files <- paths[!file.exists(paths)]
if (length(missing_files)) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "),
       ". Attach the 'synthetic-tumor-cohort' and 'synthetic-gene-annotation' data assets.", call. = FALSE)
}

read_tsv <- function(p) read.delim(p, sep = "\t", header = TRUE, check.names = FALSE,
                                   stringsAsFactors = FALSE, quote = "")
annot    <- read_tsv(paths[["annotation"]])
clinical <- read_tsv(paths[["clinical"]])
expr_df  <- read_tsv(paths[["expression"]])

expr <- as.matrix(expr_df[, setdiff(names(expr_df), "gene_id"), drop = FALSE])
rownames(expr) <- expr_df$gene_id

# clinical columns needed
need <- unique(c("patient_id", time_col, event_col, covariates))
absent <- setdiff(need, names(clinical))
if (length(absent)) {
  stop("Column(s) not found in clinical.tsv: ", paste(absent, collapse = ", "),
       ". Available columns: ", paste(names(clinical), collapse = ", "), call. = FALSE)
}
if (!all(clinical[[event_col]] %in% c(0, 1, NA))) {
  stop("Event column '", event_col, "' must be coded 1 = event, 0 = censored.", call. = FALSE)
}

# align patients between expression and clinical
common <- intersect(clinical$patient_id, colnames(expr))
if (length(common) == 0) stop("No patient_id values match between clinical.tsv and expression.tsv.", call. = FALSE)
clinical <- clinical[match(common, clinical$patient_id), , drop = FALSE]
expr <- expr[, common, drop = FALSE]

# complete cases for everything used in the models
used_cols <- unique(c(time_col, event_col, covariates))
keep <- complete.cases(clinical[, used_cols, drop = FALSE])
n_dropped <- sum(!keep)
clinical <- clinical[keep, , drop = FALSE]
expr <- expr[, clinical$patient_id, drop = FALSE]
n_patients <- nrow(clinical)
n_events <- sum(clinical[[event_col]] == 1)
cat(sprintf("Cohort: %d patients, %d events (%d patients dropped for missing values), %d genes measured\n\n",
            n_patients, n_events, n_dropped, nrow(expr)))

# covariates: character columns become factors (stage keeps I < II < III < IV order for the reference level)
for (cv in covariates) {
  if (is.character(clinical[[cv]])) {
    lv <- unique(clinical[[cv]])
    if (cv == "stage") lv <- intersect(c("I", "II", "III", "IV"), lv)
    clinical[[cv]] <- factor(clinical[[cv]], levels = if (cv == "stage") lv else sort(lv))
  }
}

# ---------------------------------------------------------------- resolve genes
cohort_ids <- rownames(expr)
cohort_annot <- annot[match(cohort_ids, annot$gene_id), c("gene_id", "symbol"), drop = FALSE]

resolve_gene <- function(query) {
  q <- toupper(query)
  hit <- which(toupper(annot$symbol) == q | toupper(annot$gene_id) == q)
  if (length(hit) == 0) {
    return(list(ok = FALSE, msg = sprintf("'%s' is not in the annotation (synthetic-gene-annotation).", query)))
  }
  gid <- annot$gene_id[hit[1]]
  if (!gid %in% cohort_ids) {
    return(list(ok = FALSE, msg = sprintf(
      "'%s' (%s) is annotated but not measured in this %d-gene cohort.", query, gid, length(cohort_ids))))
  }
  list(ok = TRUE, gene_id = gid, symbol = annot$symbol[hit[1]])
}

close_matches <- function(query, n = 5) {
  q <- toupper(query)
  syms <- cohort_annot$symbol
  d <- drop(adist(q, toupper(syms)))
  sub_hit <- grepl(q, toupper(syms), fixed = TRUE) | startsWith(toupper(syms), substr(q, 1, 2))
  ord <- order(d - 2 * sub_hit, syms)
  head(syms[ord], n)
}

resolved <- lapply(gene_queries, resolve_gene)
bad <- which(!vapply(resolved, function(r) r$ok, logical(1)))
if (length(bad)) {
  for (b in bad) {
    message("ERROR: ", resolved[[b]]$msg)
    message("  Closest genes in the cohort: ", paste(close_matches(gene_queries[b]), collapse = ", "))
  }
  known_syms <- cohort_annot$symbol[!grepl("^SYNG", cohort_annot$symbol)]
  message("  Named genes available in this cohort: ", paste(known_syms, collapse = " "))
  stop("Could not resolve ", length(bad), " gene(s); nothing was run. Use a gene symbol (case-insensitive) or a gene_id.",
       call. = FALSE)
}

# ---------------------------------------------------------------- analysis helpers
fmt_p <- function(p) if (is.na(p)) "p = NA" else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

km_curve_df <- function(fit, labels) {
  strata_n <- fit$strata
  idx <- rep(seq_along(strata_n), strata_n)
  df <- data.frame(time = fit$time, surv = fit$surv, lower = fit$lower, upper = fit$upper,
                   n.censor = fit$n.censor, g = idx)
  out <- lapply(seq_along(strata_n), function(k) {
    d <- df[df$g == k, ]
    rbind(data.frame(time = 0, surv = 1, lower = 1, upper = 1, n.censor = 0, g = k), d)
  })
  out <- do.call(rbind, out)
  out$group <- factor(labels[out$g], levels = labels)
  out
}

analyse_gene <- function(res) {
  gid <- res$gene_id; sym <- res$symbol
  x <- as.numeric(expr[gid, ])
  dat <- clinical
  dat$expr <- x
  dat$time <- dat[[time_col]]
  dat$event <- dat[[event_col]]

  # high / low groups
  if (split_mode == "median") {
    cut_hi <- cut_lo <- median(x)
    dat$group <- ifelse(x > cut_hi, "High", "Low")
  } else {
    cut_lo <- unname(quantile(x, 0.25)); cut_hi <- unname(quantile(x, 0.75))
    dat$group <- ifelse(x >= cut_hi, "High", ifelse(x <= cut_lo, "Low", NA))
  }
  sub <- dat[!is.na(dat$group), ]
  n_high <- sum(sub$group == "High"); n_low <- sum(sub$group == "Low")
  sub$group <- factor(sub$group, levels = c("High", "Low"))

  # Kaplan-Meier + log-rank
  fit <- survfit(Surv(time, event) ~ group, data = sub)
  lr <- tryCatch(survdiff(Surv(time, event) ~ group, data = sub), error = function(e) NULL)
  logrank_p <- if (is.null(lr)) NA_real_ else pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

  # Cox on continuous expression (all patients), adjusted for covariates
  rhs <- paste(c("expr", covariates), collapse = " + ")
  cox <- coxph(as.formula(paste("Surv(time, event) ~", rhs)), data = dat)
  sm <- summary(cox)
  hr <- sm$conf.int["expr", ]
  cox_p <- sm$coefficients["expr", "Pr(>|z|)"]

  # KM plot
  labels <- c(sprintf("High (n=%d)", n_high), sprintf("Low (n=%d)", n_low))
  cdf <- km_curve_df(fit, labels)
  cens <- cdf[cdf$n.censor > 0, ]
  cols <- setNames(c("#D55E00", "#0072B2"), labels)
  ann <- sprintf("Log-rank %s\nCox HR = %.2f (%.2f-%.2f), %s\n(per 1 log2 unit; adjusted for %s)",
                 fmt_p(logrank_p), hr[["exp(coef)"]], hr[["lower .95"]], hr[["upper .95"]], fmt_p(cox_p),
                 if (length(covariates)) paste(covariates, collapse = ", ") else "nothing")
  split_txt <- if (split_mode == "median") "median split" else "top vs bottom quartile"
  p <- ggplot(cdf, aes(time, surv, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), stat = "identity", alpha = 0.12, colour = NA) +
    geom_step(linewidth = 0.9) +
    geom_point(data = cens, aes(time, surv), shape = 3, size = 1.8, show.legend = FALSE) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    annotate("text", x = max(cdf$time) * 0.98, y = 0.98, label = ann, hjust = 1, vjust = 1, size = 3.3) +
    labs(title = sprintf("%s expression and overall survival", sym),
         subtitle = sprintf("%s (%s); %d patients, %d events in plotted groups", gid, split_txt,
                            n_high + n_low, sum(sub$event)),
         x = sprintf("Time (%s)", time_col), y = "Survival probability", colour = NULL, fill = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  safe_sym <- gsub("[^A-Za-z0-9._-]", "_", sym)
  png_name <- paste0("km_", safe_sym, ".png")
  ggsave(file.path(RESULTS_DIR, png_name), p, width = 7, height = 5, dpi = 150, device = "png")

  row <- data.frame(gene = gid, symbol = sym,
                    HR = unname(hr[["exp(coef)"]]), CI_low = unname(hr[["lower .95"]]), CI_high = unname(hr[["upper .95"]]),
                    p = unname(cox_p), logrank_p = logrank_p, n_high = n_high, n_low = n_low,
                    stringsAsFactors = FALSE)
  list(row = row, png = png_name, cut_lo = cut_lo, cut_hi = cut_hi)
}

# ---------------------------------------------------------------- run
results <- list()
for (r in resolved) {
  cat(sprintf("Analysing %s (%s)...\n", r$symbol, r$gene_id))
  results[[r$symbol]] <- withCallingHandlers(
    analyse_gene(r),
    warning = function(w) { message("  warning: ", conditionMessage(w)); invokeRestart("muffleWarning") })
}

tab <- do.call(rbind, lapply(results, function(x) x$row))
rownames(tab) <- NULL
out_tab <- tab
num_cols <- c("HR", "CI_low", "CI_high", "p", "logrank_p")
out_tab[num_cols] <- lapply(out_tab[num_cols], signif, digits = 4)
write.csv(out_tab, file.path(RESULTS_DIR, "cox_results.csv"), row.names = FALSE, quote = FALSE)

# summary.txt
lines <- c(
  "Gene expression survival analysis (Kaplan-Meier + Cox)",
  strrep("=", 54),
  sprintf("Genes:        %s", paste(tab$symbol, collapse = ", ")),
  sprintf("Split:        %s%s", split_mode,
          if (split_mode == "quartile") " (top vs bottom quartile, middle 50% dropped)" else " (above vs at/below median)"),
  sprintf("Covariates:   %s", if (length(covariates)) paste(covariates, collapse = ", ") else "none"),
  sprintf("Time / event: %s / %s", time_col, event_col),
  sprintf("Cohort:       %d patients, %d events", n_patients, n_events),
  "",
  "Cox HR is per 1 log2 unit of continuous expression (all patients). HR > 1 = higher expression, worse survival.",
  "Log-rank p compares the high vs low groups shown in the KM plot.",
  ""
)
for (i in seq_len(nrow(tab))) {
  t <- tab[i, ]
  dir_txt <- if (t$p < 0.05) {
    if (t$HR > 1) "higher expression associated with WORSE survival" else "higher expression associated with BETTER survival"
  } else "no significant association at p < 0.05"
  lines <- c(lines,
    sprintf("%s (%s)", t$symbol, t$gene),
    sprintf("  Cox HR = %.3f (95%% CI %.3f-%.3f), %s", t$HR, t$CI_low, t$CI_high, fmt_p(t$p)),
    sprintf("  Log-rank %s; n_high = %d, n_low = %d", fmt_p(t$logrank_p), t$n_high, t$n_low),
    sprintf("  -> %s", dir_txt),
    "")
}
writeLines(lines, file.path(RESULTS_DIR, "summary.txt"))

# run_summary.json (hand-written to avoid an extra dependency)
jstr <- function(x) paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x)), '"')
jarr <- function(x) paste0("[", paste(jstr(x), collapse = ", "), "]")
png_files <- vapply(results, function(x) x$png, character(1), USE.NAMES = FALSE)
json <- c(
  "{",
  '  "capsule": "gene-survival-km",',
  '  "parameters": {',
  sprintf('    "gene": %s,', jstr(paste(gene_queries, collapse = ","))),
  sprintf('    "split": %s,', jstr(split_mode)),
  sprintf('    "covariates": %s,', jstr(if (length(covariates)) paste(covariates, collapse = ",") else "none")),
  sprintf('    "time_col": %s,', jstr(time_col)),
  sprintf('    "event_col": %s', jstr(event_col)),
  "  },",
  sprintf('  "genes_analysed": %s,', jarr(tab$symbol)),
  sprintf('  "n_patients": %d,', n_patients),
  sprintf('  "n_events": %d,', n_events),
  sprintf('  "output_files": %s', jarr(c(png_files, "cox_results.csv", "summary.txt", "run_summary.json"))),
  "}"
)
writeLines(json, file.path(RESULTS_DIR, "run_summary.json"))

cat("\n")
cat(paste(readLines(file.path(RESULTS_DIR, "cox_results.csv")), collapse = "\n"), "\n", sep = "")
cat(sprintf("\nDone. Wrote %d KM plot(s), cox_results.csv, summary.txt, run_summary.json to %s\n",
            length(png_files), RESULTS_DIR))
