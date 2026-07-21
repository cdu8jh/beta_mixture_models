#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 03b_plot_two_beta_mixture_raw_ml.R <two_beta_output_dir>")
out_dir <- args[1]
if (!dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir)

posterior_files <- list.files(
    out_dir,
    pattern = "^two_beta_mixture_posteriors_.*\\.tsv$",
    full.names = TRUE
)
if (!length(posterior_files)) stop("No two-Beta posterior files found in: ", out_dir)

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

for (posterior_file in posterior_files) {
    mod_label <- basename(posterior_file) |>
        str_remove("^two_beta_mixture_posteriors_") |>
        str_remove("\\.tsv$")
    posts <- read_tsv(posterior_file, show_col_types = FALSE)

    required <- c("sample_id", "ml_probability", "prob_modified_component")
    missing <- setdiff(required, names(posts))
    if (length(missing)) stop("Missing columns in ", posterior_file, ": ", paste(missing, collapse = ", "))

    diagnostics <- posts %>%
        group_by(sample_id) %>%
        summarise(
            n = n(),
            min_raw = min(ml_probability, na.rm = TRUE),
            median_raw = median(ml_probability, na.rm = TRUE),
            max_raw = max(ml_probability, na.rm = TRUE),
            n_unique_raw = n_distinct(ml_probability),
            min_posterior = min(prob_modified_component, na.rm = TRUE),
            median_posterior = median(prob_modified_component, na.rm = TRUE),
            max_posterior = max(prob_modified_component, na.rm = TRUE),
            n_unique_posterior = n_distinct(prob_modified_component),
            .groups = "drop"
        )
    write_tsv(diagnostics, file.path(plot_dir, paste0("plot_input_diagnostics_", mod_label, ".tsv")))

    if (all(posts$ml_probability == 1, na.rm = TRUE)) {
        stop(
            "Every raw ml_probability is 1 for ", mod_label,
            ". This is an upstream TSV/extraction problem, not a ggplot problem."
        )
    }

    # No jitter: plotted x values are exactly the raw values stored in the TSV.
    p <- ggplot(posts, aes(x = ml_probability, y = prob_modified_component)) +
        geom_point(alpha = 0.20, size = 0.30) +
        facet_wrap(~sample_id, ncol = 2) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        theme_bw() +
        labs(
            title = paste0("Posterior modified probability vs modkit probability: ", mod_label),
            subtitle = "Separated by sample; no jitter or transformation of raw probabilities",
            x = "Original modkit ML probability (ML / 255)",
            y = "Posterior probability of modified component"
        )

    ggsave(
        file.path(plot_dir, paste0("posterior_vs_modkit_by_sample_", mod_label, ".pdf")),
        p, width = 10, height = 8
    )
}
message("Done.")
