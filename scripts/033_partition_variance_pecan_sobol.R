#!/usr/bin/env Rscript
# 033_partition_variance_pecan_sobol.R
# PEcAn-backed variance decomposition: consume sobol.indices.*.Rdata and
# partition forecast variance by PEcAn root input types before applying the
# existing local parameter drilldown within the `param` component.

library(PEcAn.settings, include.only = "read.settings")
library(readr, include.only = c("write_csv"))
library(yaml, include.only = "read_yaml")
library(PEcAn.logger)

source("R/variance_decomposition.R")
source("R/variance_decomposition_pecan_sobol.R")

opts <- list(
  optparse::make_option(c("-c", "--config"),
    default = "000-config.yml",
    help = "Path to project config YAML [default: %default]"
  ),
  optparse::make_option(c("-p", "--pecan-outdir"),
    default = NULL,
    help = "Override PEcAn output directory from config"
  ),
  optparse::make_option(c("-f", "--force"),
    action = "store_true", default = FALSE,
    help = "Overwrite existing output"
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = opts))
cfg <- yaml::read_yaml(args$config)

cfg_paths <- cfg$default$paths
data_dir <- if (!is.null(cfg_paths$data_dir)) cfg_paths$data_dir else "data"
pecan_outdir <- if (!is.null(args$pecan_outdir)) {
  args$pecan_outdir
} else if (!is.null(cfg_paths$pecan_outdir)) {
  cfg_paths$pecan_outdir
} else {
  "output"
}

pecan_data_dir <- file.path(data_dir, "pecan_sobol")
final_output <- file.path(pecan_data_dir, "variance_partition_site_level.csv")

if (!args$force && file.exists(final_output)) {
  PEcAn.logger::logger.info(
    "Output exists: ", final_output, ". Use --force to regenerate. Skipping."
  )
  quit(save = "no", status = 0)
}

if (!dir.exists(pecan_data_dir)) {
  dir.create(pecan_data_dir, recursive = TRUE)
}

sobol_indices <- load_pecan_sobol_indices(pecan_outdir)
readr::write_csv(sobol_indices, file.path(pecan_data_dir, "type_level_indices.csv"))

settings_path <- file.path(pecan_outdir, "pecan.CONFIGS.xml")
if (!file.exists(settings_path)) {
  PEcAn.logger::logger.severe("Settings not found: ", settings_path)
}
settings <- PEcAn.settings::read.settings(settings_path)

ensemble_variance <- calculate_ensemble_variance(
  output_dir = pecan_outdir,
  run_ids = unique(sobol_indices$runid),
  settings = settings
)

if (is.null(ensemble_variance) || nrow(ensemble_variance) == 0) {
  PEcAn.logger::logger.severe("No ensemble variance computed.")
}

readr::write_csv(
  ensemble_variance,
  file.path(pecan_data_dir, "ensemble_variance.csv")
)

variance_partition_site <- partition_pecan_type_variance(
  sobol_indices = sobol_indices,
  ensemble_variance = ensemble_variance
)

readr::write_csv(variance_partition_site, final_output)

local_sa <- safe_read_csv(file.path(data_dir, "aggregated_sensitivity.csv"))

variance_partition_params <- partition_parameter_variance_local(
  local_sa = local_sa,
  variance_partition_site = variance_partition_site,
  ensemble_variance = ensemble_variance,
  parameter_category = "param"
)

readr::write_csv(
  variance_partition_params,
  file.path(pecan_data_dir, "variance_partition_parameters.csv")
)

type_summary <- variance_partition_site |>
  dplyr::filter(.data$category != "interaction") |>
  dplyr::summarize(
    mean_frac = mean(.data$frac_of_total, na.rm = TRUE),
    median_frac = median(.data$frac_of_total, na.rm = TRUE),
    sd_frac = sd(.data$frac_of_total, na.rm = TRUE),
    n_sites = dplyr::n(),
    .by = c("variable", "category")
  ) |>
  dplyr::arrange(.data$variable, dplyr::desc(.data$mean_frac))

readr::write_csv(
  type_summary,
  file.path(pecan_data_dir, "variance_partition_summary.csv")
)

param_summary <- variance_partition_params |>
  dplyr::filter(!is.na(.data$Var_parameter_param)) |>
  dplyr::summarize(
    mean_frac = mean(.data$frac_of_total, na.rm = TRUE),
    median_frac = median(.data$frac_of_total, na.rm = TRUE),
    sd_frac = sd(.data$frac_of_total, na.rm = TRUE),
    n_sites = dplyr::n(),
    .by = c("response_var", "parameter")
  ) |>
  dplyr::arrange(.data$response_var, dplyr::desc(.data$mean_frac))

readr::write_csv(
  param_summary,
  file.path(pecan_data_dir, "variance_partition_params_summary.csv")
)

PEcAn.logger::logger.info(
  "PEcAn-backed variance decomposition complete: ", pecan_data_dir
)
