#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Rsamtools)
    library(data.table)
    library(readr)
    library(dplyr)
    library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Usage: Rscript 02_fit_beta_mixture_raw_ml.R <config_path.R>")
}

source(args[1])

ml_prob_base <- "/dcs04/hicks/data/sparthib/drna_retinalrs/processed_data/04_pileup/01_ml_probability"

out_dir <- file.path(
    ml_prob_base,
    paste0(config$basecaller, "_", config$align),
    "02_beta_mixture_raw_ml"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output dir: ", out_dir)

samples <- data.frame(
    sample_id = config$sample_ids,
    condition = config$conditions,
    stringsAsFactors = FALSE
)

mod_names <- config$mod_names

min_prob <- 1e-6
max_prob <- 1 - 1e-6
min_calls_per_mod <- 500
max_calls_per_sample_mod <- 50000

get_mod_block_sizes <- function(mm_str) {
    if (is.na(mm_str)) {
        return(NULL)
    }

    blocks <- strsplit(mm_str, ";", fixed = TRUE)[[1]]
    blocks <- blocks[nzchar(blocks)]

    result <- vector("list", length(blocks))

    for (k in seq_along(blocks)) {
        m <- regmatches(
            blocks[k],
            regexec("^[ACGT]\\+([^.,]+)\\.?,?(.*)", blocks[k])
        )[[1]]

        if (length(m) < 3) next

        mod_code <- gsub("[.?]", "", m[2])
        deltas_str <- m[3]

        n_pos <- if (nzchar(deltas_str)) {
            length(strsplit(deltas_str, ",", fixed = TRUE)[[1]])
        } else {
            0L
        }

        result[[k]] <- list(
            mod_code = mod_code,
            n_pos = n_pos
        )
    }

    Filter(Negate(is.null), result)
}

reservoir_update <- function(existing, new_values, max_n) {
    new_values <- new_values[is.finite(new_values)]
    new_values <- new_values[new_values > 0 & new_values < 1]

    if (length(new_values) == 0) {
        return(existing)
    }

    combined <- c(existing, new_values)

    if (length(combined) > max_n) {
        combined <- sample(combined, max_n)
    }

    combined
}

weighted_beta_mle <- function(
    x,
    w,
    alpha_start,
    beta_start,
    max_concentration = 200,
    min_mean = 0.001,
    max_mean = 0.98
) {
    if (sum(w) <= 1e-8) {
        return(c(alpha_start, beta_start))
    }

    neg_loglik <- function(par) {
        a <- exp(par[1])
        b <- exp(par[2])
        -sum(w * dbeta(x, a, b, log = TRUE))
    }

    opt <- optim(
        par = log(c(alpha_start, beta_start)),
        fn = neg_loglik,
        method = "BFGS",
        control = list(maxit = 1000)
    )

    a <- exp(opt$par[1])
    b <- exp(opt$par[2])

    mu <- a / (a + b)
    conc <- a + b

    mu <- min(max(mu, min_mean), max_mean)
    conc <- min(conc, max_concentration)

    a <- mu * conc
    b <- (1 - mu) * conc

    c(a, b)
}

fit_beta_uniform_mixture <- function(x, n_init = 10, max_iter = 500, tol = 1e-6) {
    x <- pmin(pmax(x, min_prob), max_prob)
    x <- x[is.finite(x)]

    best <- NULL

    for (init in seq_len(n_init)) {
        pi_low <- runif(1, 0.3, 0.8)
        pi_high <- runif(1, 0.01, min(0.4, 1 - pi_low - 0.05))
        pi_uniform <- max(0.05, 1 - pi_low - pi_high)

        alpha_low <- runif(1, 0.5, 2)
        beta_low <- runif(1, 3, 10)

        alpha_high <- runif(1, 3, 10)
        beta_high <- runif(1, 0.5, 2)

        loglik_old <- -Inf

        for (iter in seq_len(max_iter)) {
            d_low <- pi_low * dbeta(x, alpha_low, beta_low)
            d_uniform <- pi_uniform * dunif(x, 0, 1)
            d_high <- pi_high * dbeta(x, alpha_high, beta_high)

            denom <- pmax(d_low + d_uniform + d_high, 1e-300)

            r_low <- d_low / denom
            r_uniform <- d_uniform / denom
            r_high <- d_high / denom

            pi_low <- mean(r_low)
            pi_uniform <- mean(r_uniform)
            pi_high <- mean(r_high)

            low_params <- weighted_beta_mle(
                x,
                r_low,
                alpha_low,
                beta_low,
                max_concentration = 200,
                min_mean = 0.001,
                max_mean = 0.40
            )

            high_params <- weighted_beta_mle(
                x,
                r_high,
                alpha_high,
                beta_high,
                max_concentration = 100,
                min_mean = 0.70,
                max_mean = 0.98
            )

            alpha_low <- low_params[1]
            beta_low <- low_params[2]

            alpha_high <- high_params[1]
            beta_high <- high_params[2]

            loglik <- sum(log(denom))

            if (abs(loglik - loglik_old) < tol) {
                break
            }

            loglik_old <- loglik
        }

        if (is.null(best) || loglik > best$logLik) {
            posterior <- cbind(
                prob_unmodified_component = r_low,
                prob_ambiguous_uniform_component = r_uniform,
                prob_modified_component = r_high
            )

            params <- tibble(
                component = c(
                    "unmodified_low",
                    "ambiguous_uniform",
                    "modified_high"
                ),
                distribution = c("beta", "uniform", "beta"),
                weight = c(pi_low, pi_uniform, pi_high),
                alpha = c(alpha_low, NA_real_, alpha_high),
                beta = c(beta_low, NA_real_, beta_high),
                mean = c(
                    alpha_low / (alpha_low + beta_low),
                    0.5,
                    alpha_high / (alpha_high + beta_high)
                )
            )

            best <- list(
                params = params,
                posterior = posterior,
                x = x,
                logLik = loglik,
                iterations = iter
            )
        }
    }

    best
}

all_probs <- list()

for (mc in names(mod_names)) {
    all_probs[[mc]] <- tibble(
        sample_id = character(),
        ml_probability = numeric()
    )
}

for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]

    bam_path <- file.path(config$bam_dir, config$bam_filename(sid))

    if (!file.exists(bam_path)) {
        stop("BAM not found: ", bam_path)
    }

    message("Processing: ", bam_path)

    sample_probs <- list()

    for (mc in names(mod_names)) {
        sample_probs[[mc]] <- numeric(0)
    }

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

    repeat {
        chunk <- scanBam(bf, param = param)[[1]]
        n_reads <- length(chunk$tag$MM)

        if (n_reads == 0L) break

        mm_tags <- chunk$tag$MM
        ml_tags <- chunk$tag$ML

        for (j in seq_len(n_reads)) {
            mod_info <- get_mod_block_sizes(mm_tags[j])

            if (is.null(mod_info) || length(mod_info) == 0) next

            ml_vec <- ml_tags[[j]]

            if (is.null(ml_vec) || length(ml_vec) == 0) next

            ml_idx <- 1L

            for (mi in mod_info) {
                mc <- mi$mod_code
                np <- mi$n_pos

                if (np == 0L) next

                end_idx <- ml_idx + np - 1L

                if (end_idx > length(ml_vec)) break

                probs <- ml_vec[ml_idx:end_idx] / 255
                probs <- pmin(pmax(probs, min_prob), max_prob)

                if (mc %in% names(sample_probs)) {
                    sample_probs[[mc]] <- reservoir_update(
                        sample_probs[[mc]],
                        probs,
                        max_calls_per_sample_mod
                    )
                }

                ml_idx <- end_idx + 1L
            }
        }

        current_counts <- vapply(sample_probs, length, integer(1))

        message(
            sid,
            " retained calls so far: ",
            paste(names(current_counts), current_counts, sep = "=", collapse = ", ")
        )

        if (all(current_counts >= max_calls_per_sample_mod)) {
            message(
                "Reached ",
                max_calls_per_sample_mod,
                " calls for all mods in sample ",
                sid,
                "; stopping early."
            )
            break
        }
    }

    close(bf)

    for (mc in names(mod_names)) {
        all_probs[[mc]] <- bind_rows(
            all_probs[[mc]],
            tibble(
                sample_id = sid,
                ml_probability = sample_probs[[mc]]
            )
        )

        message(
            sid,
            " ",
            mod_names[[mc]],
            ": retained ",
            length(sample_probs[[mc]]),
            " raw ML probabilities"
        )
    }
}

summary_rows <- list()

for (mc in names(mod_names)) {
    mod_label <- mod_names[[mc]]
    prob_df <- all_probs[[mc]]
    x <- prob_df$ml_probability

    message("Fitting Beta + Uniform + Beta mixture for ", mod_label)

    if (length(x) < min_calls_per_mod) {
        warning(
            "Skipping ",
            mod_label,
            ": fewer than ",
            min_calls_per_mod,
            " calls"
        )
        next
    }

    fit <- fit_beta_uniform_mixture(x)

    params_out <- fit$params %>%
        mutate(
            basecaller = config$basecaller,
            align = config$align,
            modification = mod_label,
            mod_code = mc,
            n_calls_used = length(fit$x),
            logLik = fit$logLik,
            iterations = fit$iterations
        )

    posterior_out <- tibble(
        sample_id = prob_df$sample_id,
        modification = mod_label,
        mod_code = mc,
        ml_probability = fit$x,
        prob_unmodified_component =
            fit$posterior[, "prob_unmodified_component"],
        prob_ambiguous_uniform_component =
            fit$posterior[, "prob_ambiguous_uniform_component"],
        prob_modified_component =
            fit$posterior[, "prob_modified_component"],
        assigned_component = c(
            "unmodified_low",
            "ambiguous_uniform",
            "modified_high"
        )[max.col(fit$posterior)],
        high_conf_modified = prob_modified_component >= 0.90,
        high_conf_unmodified = prob_unmodified_component >= 0.90
    )

    write_tsv(
        params_out,
        file.path(out_dir, paste0("beta_mixture_params_", mod_label, ".tsv"))
    )

    write_tsv(
        posterior_out,
        file.path(out_dir, paste0("beta_mixture_posteriors_", mod_label, ".tsv"))
    )

    summary_rows[[length(summary_rows) + 1L]] <- params_out
}

summary_dt <- bind_rows(summary_rows)

write_tsv(
    summary_dt,
    file.path(out_dir, "beta_mixture_params_all_mods.tsv")
)

message("Done.")