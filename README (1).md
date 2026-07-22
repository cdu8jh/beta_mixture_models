# RNA Modification Probability Mixture Models

This repository fits mixture models to per-read RNA modification probabilities (for example, m6A probabilities reported by `modkit`).  The models separate a low-probability, likely-unmodified population from a high-probability, likely-modified population and an ambiguous background population.  Their main outputs are posterior modification probabilities, local-FDR-based probability cutoffs, and coverage-aware site-level modification estimates.

The repository is intended for Oxford Nanopore direct RNA-sequencing data, but the statistical workflow is applicable to any per-read probability scores constrained to the interval `(0, 1)`.

## Model

The default model has three components:

\[
f(p) = \pi_0\operatorname{Beta}(p;\alpha_0,\beta_0) +
       \pi_1\operatorname{Beta}(p;\alpha_1,\beta_1) +
       \pi_u\operatorname{Uniform}(p;0,1)
\]

where:

| Component | Interpretation |
| --- | --- |
| `Beta(α0, β0)` | Low probabilities; likely unmodified reads |
| `Beta(α1, β1)` | High probabilities; likely modified reads |
| `Uniform(0,1)` | Ambiguous or poorly separated reads |

The model is fitted by expectation-maximization (EM). For each read probability `p`, it estimates a posterior probability of belonging to the modified component. The uniform component prevents uncertain observations from being forced into either biological class.

## Repository layout

```text
.
├── 01_fit_global_beta_mixture.R
├── 02_fit_site_beta_mixture.R
├── 03_plot_beta_mixture.R
├── run_global_beta_mixture.sh
├── run_site_beta_mixture.sh
└── README.md
```

Script names may differ slightly in a downstream project; the expected roles are:

| Script | Purpose |
| --- | --- |
| Global fitting script | Fits one mixture model using all eligible reads for a modification/sample set. Useful for QC and an overall probability calibration. |
| Site-level fitting script | Fits the same model independently at each genomic or transcriptomic site. Use this for site-level calls and rates. |
| Plotting script | Creates probability histograms and raw-versus-posterior diagnostic scatterplots. |
| Shell scripts | Submit the R workflows to a SLURM cluster or run them reproducibly from the command line. |

## Requirements

Install R (version 4.2 or later recommended) and the packages used by the scripts. Typical dependencies are:

```r
install.packages(c("data.table", "readr", "dplyr", "ggplot2", "stringr"))
```

The workflow expects tab-delimited input and does not require access to raw FAST5/POD5 files. Upstream basecalling, alignment, and `modkit` extraction must already be complete.

## Input data

Each input row should represent one read-level probability call at one candidate modification site. The fitting scripts need a numeric probability column with values from 0 to 1, commonly named `ml_prob`, `probability`, or an equivalent configured in the script.

For site-level analysis, input must also include identifiers that define a site, such as:

```text
chrom    position    strand    sample    ml_prob
```

If working on transcript coordinates, use a stable transcript identifier and transcript position instead of chromosome and genomic position. Do not mix coordinate systems in the same analysis.

Before fitting, remove missing, non-finite, and out-of-range probabilities. Values exactly equal to 0 or 1 should be clipped very slightly inward (for example, to `1e-6` and `1 - 1e-6`) because beta densities are evaluated on the open interval `(0, 1)`.

## Quick start

### Global model

Use the global model first to understand the overall distribution of basecaller probabilities and assess whether the components are well separated.

```bash
Rscript 01_fit_global_beta_mixture.R <input.tsv> <output_directory>
Rscript 03_plot_beta_mixture.R <output_directory>
```

For a SLURM workflow, adapt the paths and resources in the supplied submission script, then run:

```bash
sbatch run_global_beta_mixture.sh
```

### Site-level model

Use the site-level model when the goal is to call or quantify modification at individual sites.

```bash
Rscript 02_fit_site_beta_mixture.R <input.tsv> <output_directory>
```

The site-level model should enforce a minimum number of reads per site. A reasonable initial value is `20`, but it should be selected based on the experiment’s coverage distribution and computational budget. Sites below the threshold should be retained in an output status table as `insufficient_coverage`, rather than treated as confidently unmodified.

## Key outputs

Exact filenames depend on the scripts, but the analyses produce these concepts:

| Output | Meaning |
| --- | --- |
| Fitted parameters | Mixture proportions (`π`) and beta parameters (`α`, `β`) for each fit. |
| Per-read posterior table | Original probability plus posterior probability that the read is modified. |
| Local-FDR table / cutoff | A data-driven probability threshold at a chosen local FDR target. |
| Site summary table | Coverage, fitted status, cutoff, posterior-weighted modification rate, and other site-level metrics. |
| Diagnostic plots | Raw probability histograms and raw-versus-posterior scatterplots. |

## Calling reads with local FDR

For a read with probability `p`, the local false discovery rate is calculated as:

\[
\operatorname{lfdr}(p) =
\frac{\pi_0 f_0(p) + \pi_u f_u(p)}
     {\pi_0 f_0(p) + \pi_1 f_1(p) + \pi_u f_u(p)}
\]

The cutoff `τ` is the smallest probability for which `lfdr(p) ≤ α`, where `α` is the desired local-FDR target (commonly 0.10). Reads with `p ≥ τ` are treated as modified for a thresholded analysis.

This is preferable to choosing a fixed raw-probability threshold because it accounts for the observed probability distribution and the estimated ambiguous background. The cutoff is nevertheless model-dependent and should be reported with the fitted parameters and diagnostics.

## Posterior-weighted site modification rate

Rather than simply counting reads above a cutoff, the site-level workflow can estimate a soft modification rate:

\[
\hat{\theta}_{site} = \frac{\sum_i w_i p_i}{\sum_i w_i},
\qquad w_i = 1 - \operatorname{lfdr}(p_i)
\]

Here, high-confidence modified reads receive large weights, while ambiguous and likely-unmodified reads contribute little. This provides a continuous estimate, but it is not a direct measurement of the fraction of molecules modified; it remains calibrated to the basecaller’s probability scale and model fit.

## Quality control and interpretation

Inspect the diagnostic plots for every sample or analysis batch.

- A clear low-probability mode and high-probability mode support separation of unmodified and modified reads.
- In raw-versus-posterior plots, high raw probabilities should generally correspond to high modified posterior probabilities.
- A broad intermediate cloud is expected when many calls are ambiguous; the uniform component should absorb much of it.
- If beta components pile up at 0 or 1, overlap heavily, or vary strongly between EM initializations, the fit may be unstable. Treat resulting thresholds cautiously.
- A failed or insufficient-coverage site is not evidence that the site is unmodified.

The global model is primarily a calibration and QC tool. For biological comparisons of specific genomic/transcriptomic locations, use site-level estimates and account for coverage, replicate structure, and multiple testing in downstream analyses.

## Recommended analysis sequence

1. Extract per-read m6A probabilities from aligned reads.
2. Run the global model by sample/basecaller to inspect probability distributions and choose an initial local-FDR target.
3. Review histograms, scatterplots, fitted parameters, and convergence diagnostics.
4. Run the site-level model using a prespecified minimum read threshold.
5. Filter or label sites by fit status and coverage before comparing conditions.
6. Use posterior-weighted rates or thresholded calls in downstream reproducibility and differential-modification analyses.

## Limitations

- The model classifies basecaller probabilities; it cannot independently validate a chemical modification.
- Results depend on basecaller calibration, sequence context, read quality, and coverage.
- A fixed three-component form may not fit every batch equally well.
- Per-site fits with few reads are inherently unstable, even when they converge.
- Technical replicates and repeated observations from the same biological material are not automatically independent biological replicates.

## Reproducibility

For each run, retain the input manifest, script version or Git commit, R session information, command line, model settings (`alpha`, minimum reads, probability-column name), fitted parameter tables, and diagnostic plots. Set a random seed before multi-start EM fitting so results can be reproduced exactly.

## Citation

If this repository supports a publication, cite the basecaller and `modkit` version used to generate probabilities, along with the repository commit/version and the statistical framework described here.
