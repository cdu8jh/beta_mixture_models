RNA Modification Probability Mixture Models

This repository fits mixture models to per-read RNA modification probabilities (for example, m6A probabilities reported by modkit).  The models separate a low-probability, likely-unmodified population from a high-probability, likely-modified population and an ambiguous background population.  Their main outputs are posterior modification probabilities, local-FDR-based probability cutoffs, and coverage-aware site-level modification estimates.

The repository is intended for Oxford Nanopore direct RNA-sequencing data, but the statistical workflow is applicable to any per-read probability scores constrained to the interval (0, 1).
