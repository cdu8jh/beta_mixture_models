#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Usage: Rscript 03_plot_beta_mixture_raw_ml.R <beta_mixture_output_dir>")
}

out_dir <- args[1]

if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir)
}

posterior_files <- list.files(
    out_dir,
    pattern = "^beta_mixture_posteriors_.*\\.tsv$",
    full.names = TRUE
)

if (length(posterior_files) == 0) {
    stop("No posterior files found in: ", out_dir)
}

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading posterior files from: ", out_dir)
message("Writing plots to: ", plot_dir)

for (posterior_file in posterior_files) {
    mod_label <- basename(posterior_file) |>
        str_remove("^beta_mixture_posteriors_") |>
        str_remove("\\.tsv$")

    message("Plotting ", mod_label)

    posts <- read_tsv(posterior_file, show_col_types = FALSE)

    required_cols <- c(
        "sample_id",
        "ml_probability",
        "prob_modified_component"
    )

    missing_cols <- setdiff(required_cols, colnames(posts))

    if (length(missing_cols) > 0) {
        warning(
            "Skipping ",
            mod_label,
            ": missing columns: ",
            paste(missing_cols, collapse = ", ")
        )
        next
    }

    p_scatter <- ggplot(
        posts,
        aes(
            x = ml_probability,
            y = prob_modified_component,
            color = sample_id
        )
    ) +
        geom_point(
            alpha = 0.25,
            size = 0.25,
            position = position_jitter(width = 0.0015, height = 0.0015)
        ) +
        theme_bw() +
        labs(
            title = paste0(
                "Posterior modified probability vs raw modkit probability: ",
                mod_label
            ),
            subtitle = "All samples overlaid; color indicates sample",
            x = "Original modkit ML probability",
            y = "Posterior probability of modified component",
            color = "Sample"
        ) 


    ggsave(
        file.path(
            plot_dir,
            paste0("posterior_vs_modkit_all_samples_", mod_label, ".pdf")
        ),
        p_scatter,
        width = 8,
        height = 6
    )

    message("Wrote scatterplot for ", mod_label)
}

message("Done.")
