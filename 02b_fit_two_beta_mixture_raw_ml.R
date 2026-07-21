#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(dplyr)
    library(readr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 02b_fit_two_beta_mixture_raw_ml.R <config_path.R>")
}
source(args[1])

ml_prob_base <- "/dcs04/hicks/data/sparthib/drna_retinalrs/processed_data/04_pileup/01_ml_probability"
out_dir <- file.path(
    ml_prob_base,
    paste0(config$basecaller, "_", config$align),
    "02b_two_beta_mixture_raw_ml"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Output dir: ", out_dir)

samples <- tibble(
    sample_id = config$sample_ids,
    condition = config$conditions
)
mod_names <- config$mod_names

min_prob <- 1e-6
max_prob <- 1 - 1e-6
min_calls_per_mod <- 500L
max_calls_per_sample_mod <- 50000L
set.seed(1)

get_mod_block_sizes <- function(mm_str) {
    if (length(mm_str) == 0L || is.na(mm_str) || !nzchar(mm_str)) return(NULL)
    blocks <- strsplit(mm_str, ";", fixed = TRUE)[[1]]
    blocks <- blocks[nzchar(blocks)]
    result <- lapply(blocks, function(block) {
        m <- regmatches(block, regexec("^[ACGT]\\+([^.,]+)\\.?,?(.*)", block))[[1]]
        if (length(m) < 3L) return(NULL)
        deltas <- m[3]
        list(
            mod_code = gsub("[.?]", "", m[2]),
            n_pos = if (nzchar(deltas)) length(strsplit(deltas, ",", fixed = TRUE)[[1]]) else 0L
        )
    })
    Filter(Negate(is.null), result)
}

reservoir_update <- function(existing, new_values, max_n) {
    new_values <- as.numeric(new_values)
    new_values <- new_values[is.finite(new_values)]
    if (!length(new_values)) return(existing)
    combined <- c(existing, new_values)
    if (length(combined) > max_n) combined <- sample(combined, max_n)
    combined
}

weighted_beta_mle <- function(x, w, start, mean_bounds, max_concentration = 200) {
    if (!all(is.finite(w)) || sum(w) <= 1e-8) return(start)

    objective <- function(log_par) {
        a <- exp(log_par[1])
        b <- exp(log_par[2])
        value <- -sum(w * dbeta(x, a, b, log = TRUE))
        if (is.finite(value)) value else .Machine$double.xmax / 100
    }

    opt <- tryCatch(
        optim(log(start), objective, method = "BFGS", control = list(maxit = 1000)),
        error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$value)) return(start)

    a <- exp(opt$par[1]); b <- exp(opt$par[2])
    mu <- a / (a + b)
    concentration <- min(a + b, max_concentration)
    mu <- min(max(mu, mean_bounds[1]), mean_bounds[2])
    c(mu * concentration, (1 - mu) * concentration)
}

fit_two_beta_mixture <- function(x, n_init = 20L, max_iter = 500L, tol = 1e-7) {
    x <- pmin(pmax(as.numeric(x), min_prob), max_prob)
    if (!length(x)) stop("No finite probabilities supplied to mixture model")

    best <- NULL
    for (init in seq_len(n_init)) {
        pi_low <- runif(1, 0.55, 0.90)
        pi_high <- 1 - pi_low
        low_par <- c(runif(1, 0.7, 2), runif(1, 4, 12))
        high_par <- c(runif(1, 4, 12), runif(1, 0.7, 2))
        old_loglik <- -Inf

        for (iter in seq_len(max_iter)) {
            d_low <- pi_low * dbeta(x, low_par[1], low_par[2])
            d_high <- pi_high * dbeta(x, high_par[1], high_par[2])
            denom <- pmax(d_low + d_high, 1e-300)
            r_low <- d_low / denom
            r_high <- d_high / denom

            pi_low <- min(max(mean(r_low), 0.01), 0.99)
            pi_high <- 1 - pi_low
            low_par <- weighted_beta_mle(x, r_low, low_par, c(0.001, 0.45), 200)
            high_par <- weighted_beta_mle(x, r_high, high_par, c(0.55, 0.999), 200)

            # Recompute likelihood using the updated parameters.
            denom_updated <- pmax(
                pi_low * dbeta(x, low_par[1], low_par[2]) +
                    pi_high * dbeta(x, high_par[1], high_par[2]),
                1e-300
            )
            loglik <- sum(log(denom_updated))
            if (is.finite(old_loglik) && abs(loglik - old_loglik) < tol) break
            old_loglik <- loglik
        }

        d_low <- pi_low * dbeta(x, low_par[1], low_par[2])
        d_high <- pi_high * dbeta(x, high_par[1], high_par[2])
        denom <- pmax(d_low + d_high, 1e-300)
        posterior <- cbind(
            prob_unmodified_component = d_low / denom,
            prob_modified_component = d_high / denom
        )
        loglik <- sum(log(denom))

        low_mean <- low_par[1] / sum(low_par)
        high_mean <- high_par[1] / sum(high_par)
        valid_fit <- is.finite(loglik) && low_mean < high_mean &&
            low_mean <= 0.45 && high_mean >= 0.55

        if (valid_fit && (is.null(best) || loglik > best$logLik)) {
            best <- list(
                x = x,
                posterior = posterior,
                logLik = loglik,
                iterations = iter,
                params = tibble(
                    component = c("unmodified_low", "modified_high"),
                    distribution = "beta",
                    weight = c(pi_low, pi_high),
                    alpha = c(low_par[1], high_par[1]),
                    beta = c(low_par[2], high_par[2]),
                    mean = c(low_mean, high_mean)
                )
            )
        }
    }
    if (is.null(best)) stop("All two-Beta fits failed component-separation checks")
    best
}

all_probs <- setNames(vector("list", length(mod_names)), names(mod_names))
for (mc in names(mod_names)) {
    all_probs[[mc]] <- tibble(sample_id = character(), ml_probability_raw = numeric())
}

for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    bam_path <- file.path(config$bam_dir, config$bam_filename(sid))
    if (!file.exists(bam_path)) stop("BAM not found: ", bam_path)
    message("Processing: ", bam_path)

    sample_probs <- setNames(lapply(mod_names, function(x) numeric()), names(mod_names))
    bf <- BamFile(bam_path, yieldSize = 100000L)
    param <- ScanBamParam(
        tag = c("MM", "ML"),
        flag = scanBamFlag(
            isSecondaryAlignment = FALSE,
            isSupplementaryAlignment = FALSE,
            isUnmappedQuery = FALSE,
            isDuplicate = FALSE
        )
    )

    open(bf)
    on.exit(try(close(bf), silent = TRUE), add = TRUE)
    repeat {
        chunk <- scanBam(bf, param = param)[[1]]
        n_reads <- length(chunk$tag$MM)
        if (n_reads == 0L) break

        for (j in seq_len(n_reads)) {
            mod_info <- get_mod_block_sizes(chunk$tag$MM[j])
            ml_vec <- chunk$tag$ML[[j]]
            if (is.null(mod_info) || !length(mod_info) || is.null(ml_vec) || !length(ml_vec)) next

            ml_idx <- 1L
            for (mi in mod_info) {
                np <- mi$n_pos
                if (np == 0L) next
                end_idx <- ml_idx + np - 1L
                if (end_idx > length(ml_vec)) break

                # Preserve the raw modkit probability exactly as ML / 255.
                probs_raw <- as.numeric(ml_vec[ml_idx:end_idx]) / 255
                if (mi$mod_code %in% names(sample_probs)) {
                    sample_probs[[mi$mod_code]] <- reservoir_update(
                        sample_probs[[mi$mod_code]], probs_raw, max_calls_per_sample_mod
                    )
                }
                ml_idx <- end_idx + 1L
            }
        }

        counts <- vapply(sample_probs, length, integer(1))
        if (all(counts >= max_calls_per_sample_mod)) break
    }
    close(bf)

    for (mc in names(mod_names)) {
        all_probs[[mc]] <- bind_rows(
            all_probs[[mc]],
            tibble(sample_id = sid, ml_probability_raw = sample_probs[[mc]])
        )
        message(sid, " ", mod_names[[mc]], ": retained ", length(sample_probs[[mc]]), " calls")
    }
}

summary_rows <- list()
for (mc in names(mod_names)) {
    mod_label <- mod_names[[mc]]
    prob_df <- all_probs[[mc]] %>% filter(is.finite(ml_probability_raw))
    if (nrow(prob_df) < min_calls_per_mod) {
        warning("Skipping ", mod_label, ": fewer than ", min_calls_per_mod, " calls")
        next
    }

    message("Fitting constrained two-Beta mixture for ", mod_label)
    model_x <- pmin(pmax(prob_df$ml_probability_raw, min_prob), max_prob)
    fit <- fit_two_beta_mixture(model_x)

    # The posterior is calculated from model_x, but the output retains the untouched raw value.
    posterior_out <- prob_df %>%
        mutate(
            modification = mod_label,
            mod_code = mc,
            ml_probability = ml_probability_raw,
            ml_probability_model = model_x,
            prob_unmodified_component = fit$posterior[, "prob_unmodified_component"],
            prob_modified_component = fit$posterior[, "prob_modified_component"],
            assigned_component = if_else(
                prob_modified_component >= prob_unmodified_component,
                "modified_high", "unmodified_low"
            ),
            high_conf_modified = prob_modified_component >= 0.90,
            high_conf_unmodified = prob_unmodified_component >= 0.90
        ) %>%
        select(-ml_probability_raw)

    params_out <- fit$params %>%
        mutate(
            basecaller = config$basecaller,
            align = config$align,
            modification = mod_label,
            mod_code = mc,
            n_calls_used = nrow(prob_df),
            logLik = fit$logLik,
            iterations = fit$iterations
        )

    raw_summary <- posterior_out %>%
        group_by(sample_id) %>%
        summarise(
            n = n(),
            min_raw = min(ml_probability),
            q25_raw = quantile(ml_probability, 0.25),
            median_raw = median(ml_probability),
            mean_raw = mean(ml_probability),
            q75_raw = quantile(ml_probability, 0.75),
            max_raw = max(ml_probability),
            n_unique_raw = n_distinct(ml_probability),
            fraction_raw_equal_1 = mean(ml_probability == 1),
            .groups = "drop"
        )

    write_tsv(params_out, file.path(out_dir, paste0("two_beta_mixture_params_", mod_label, ".tsv")))
    write_tsv(posterior_out, file.path(out_dir, paste0("two_beta_mixture_posteriors_", mod_label, ".tsv")))
    write_tsv(raw_summary, file.path(out_dir, paste0("raw_probability_diagnostics_", mod_label, ".tsv")))
    summary_rows[[length(summary_rows) + 1L]] <- params_out
}

write_tsv(bind_rows(summary_rows), file.path(out_dir, "two_beta_mixture_params_all_mods.tsv"))
message("Done.")
