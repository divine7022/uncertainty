#' Parse metadata from a canonical PEcAn Sobol filename
#'
#' PEcAn stores run-level identifiers in analysis filenames produced from
#' `ensemble.filename()`. Factor/type metadata is read from the saved objects
#' themselves; only file keys are parsed here.
#'
#' Expected pattern:
#'   sobol.indices.<runid>.<variable>.<start_year>.<end_year>.Rdata
#'
#' @param path File path to a Sobol indices artifact.
#' @return Named list with parsed metadata.
parse_pecan_sobol_filename <- function(path) {
  parts <- strsplit(basename(path), "\\.")[[1]]

  if (length(parts) < 7 || parts[[1]] != "sobol" || parts[[2]] != "indices") {
    PEcAn.logger::logger.severe("Unrecognized Sobol filename: ", path)
  }

  list(
    runid = parts[[3]],
    variable = parts[[4]],
    start_year = parts[[5]],
    end_year = parts[[6]],
    sobol_file = basename(path)
  )
}


#' Load factor metadata from a sibling sobol.design file
#'
#' @param indices_path Path to the sobol.indices file.
#' @return Tibble with `parameters`, `source_type`, and `source_tag` when
#'   available, otherwise `NULL`.
load_pecan_factor_metadata <- function(indices_path) {
  design_path <- file.path(
    dirname(indices_path),
    sub("^sobol\\.indices\\.", "sobol.design.", basename(indices_path))
  )

  if (!file.exists(design_path)) {
    return(NULL)
  }

  design_env <- new.env(parent = emptyenv())
  load(design_path, envir = design_env)

  if (!exists("sobol_design", envir = design_env)) {
    return(NULL)
  }

  factor_metadata <- design_env$sobol_design$factor_metadata
  if (is.null(factor_metadata)) {
    return(NULL)
  }

  tibble::as_tibble(factor_metadata) |>
    dplyr::rename(parameters = "factor") |>
    dplyr::select(dplyr::any_of(c("parameters", "source_type", "source_tag")))
}


#' Merge PEcAn factor metadata into saved Sobol results
#'
#' @param sobol_results Tibble loaded from a sobol.indices file.
#' @param factor_metadata Tibble from [load_pecan_factor_metadata()].
#' @return Tibble with `source_type` / `source_tag` filled from the PEcAn
#'   object first and the sibling design file second.
merge_pecan_factor_metadata <- function(sobol_results, factor_metadata) {
  if (!("source_type" %in% names(sobol_results))) {
    sobol_results$source_type <- NA_character_
  }
  if (!("source_tag" %in% names(sobol_results))) {
    sobol_results$source_tag <- NA_character_
  }
  if (is.null(factor_metadata)) {
    return(sobol_results)
  }

  sobol_results |>
    dplyr::left_join(factor_metadata, by = "parameters", suffix = c("", ".design")) |>
    dplyr::mutate(
      source_type = dplyr::coalesce(.data$source_type, .data$source_type.design),
      source_tag = dplyr::coalesce(.data$source_tag, .data$source_tag.design)
    ) |>
    dplyr::select(-dplyr::any_of(c("source_type.design", "source_tag.design")))
}


#' Load PEcAn Sobol indices saved as Rdata artifacts
#'
#' Reads all `sobol.indices.*.Rdata` files in a PEcAn outdir, returning the
#' long-format Sobol results table with canonical filename metadata added.
#' Factor/type metadata comes from `sobol_results` when present and falls back
#' to the sibling `sobol.design.*.Rdata` artifact only when needed.
#'
#' @param output_dir Path to a PEcAn output directory.
#' @return Tibble with Sobol indices across all matching files.
load_pecan_sobol_indices <- function(output_dir) {
  sobol_files <- list.files(
    output_dir,
    pattern = "^sobol\\.indices\\..*\\.Rdata$",
    full.names = TRUE
  )

  if (length(sobol_files) == 0) {
    PEcAn.logger::logger.severe(
      "No sobol.indices.*.Rdata files found in ", output_dir
    )
  }

  purrr::map_dfr(sobol_files, \(path) {
    sobol_env <- new.env(parent = emptyenv())
    load(path, envir = sobol_env)

    sobol_results <- if (exists("sobol_results", envir = sobol_env)) {
      tibble::as_tibble(sobol_env$sobol_results)
    } else if (exists("sobol_indices_result", envir = sobol_env)) {
      tibble::as_tibble(sobol_env$sobol_indices_result$results)
    } else {
      PEcAn.logger::logger.severe(
        "Neither `sobol_results` nor `sobol_indices_result` found in ", path
      )
    }

    if (!("parameters" %in% names(sobol_results))) {
      PEcAn.logger::logger.severe(
        "`sobol_results` is missing the `parameters` column in ", path
      )
    }

    factor_metadata <- NULL
    needs_metadata <- !("source_type" %in% names(sobol_results)) ||
      !("source_tag" %in% names(sobol_results)) ||
      any(is.na(sobol_results$source_type))

    if (needs_metadata) {
      factor_metadata <- load_pecan_factor_metadata(path)
    }
    sobol_results <- merge_pecan_factor_metadata(sobol_results, factor_metadata)

    parsed <- parse_pecan_sobol_filename(path)

    sobol_results |>
      dplyr::mutate(
        runid = parsed$runid,
        variable = parsed$variable,
        start_year = parsed$start_year,
        end_year = parsed$end_year,
        sobol_file = parsed$sobol_file,
        .before = 1
      )
  })
}


#' Partition ensemble variance using PEcAn type-level Sobol indices
#'
#' Consumes PEcAn-saved Sobol first-order indices where each factor already
#' corresponds to a root uncertainty type (`param`, `met`, `poolinitcond`,
#' `events`, etc.). Interactions are represented as the residual not explained
#' by the non-negative first-order contributions.
#'
#' @param sobol_indices Tibble from [load_pecan_sobol_indices()].
#' @param ensemble_variance Tibble from [calculate_ensemble_variance()].
#' @return Tibble with the same schema as the legacy partition output.
partition_pecan_type_variance <- function(sobol_indices, ensemble_variance) {
  required_sobol_cols <- c(
    "runid", "variable", "parameters", "sensitivity", "original"
  )
  missing_sobol <- setdiff(required_sobol_cols, names(sobol_indices))
  if (length(missing_sobol) > 0) {
    PEcAn.logger::logger.severe(
      "sobol_indices missing columns: ", paste(missing_sobol, collapse = ", ")
    )
  }

  required_var_cols <- c("runid", "site_id", "variable", "ensemble_variance")
  missing_var <- setdiff(required_var_cols, names(ensemble_variance))
  if (length(missing_var) > 0) {
    PEcAn.logger::logger.severe(
      "ensemble_variance missing columns: ", paste(missing_var, collapse = ", ")
    )
  }

  first_order <- sobol_indices |>
    dplyr::filter(.data$sensitivity == "Si") |>
    dplyr::mutate(
      category = dplyr::coalesce(.data$source_type, .data$parameters)
    )

  sobol_joined <- first_order |>
    dplyr::left_join(ensemble_variance, by = c("runid", "variable")) |>
    dplyr::mutate(
      Var_total = .data$ensemble_variance,
      first_order_index = pmax(.data$original, 0),
      Var_i = .data$first_order_index * .data$Var_total
    )

  var_cat <- sobol_joined |>
    dplyr::summarize(
      Var_category = sum(.data$Var_i, na.rm = TRUE),
      .by = c("runid", "site_id", "variable", "category")
    )

  var_total_lookup <- sobol_joined |>
    dplyr::distinct(.data$runid, .data$site_id, .data$variable, .data$Var_total)

  var_cat <- var_cat |>
    dplyr::left_join(var_total_lookup, by = c("runid", "site_id", "variable")) |>
    dplyr::mutate(
      Var_first_sum = sum(.data$Var_category, na.rm = TRUE),
      Var_interaction = pmax(.data$Var_total - .data$Var_first_sum, 0),
      .by = c("runid", "site_id", "variable")
    )

  interaction_rows <- var_cat |>
    dplyr::distinct(
      .data$runid,
      .data$site_id,
      .data$variable,
      .data$Var_total,
      .data$Var_interaction
    ) |>
    dplyr::mutate(
      category = "interaction",
      Var_category = .data$Var_interaction
    )

  var_cat |>
    dplyr::select(
      "runid", "site_id", "variable", "category",
      "Var_category", "Var_total", "Var_interaction"
    ) |>
    dplyr::bind_rows(interaction_rows) |>
    dplyr::mutate(
      frac_of_total = dplyr::if_else(
        .data$Var_total > 0,
        .data$Var_category / .data$Var_total,
        NA_real_
      )
    ) |>
    dplyr::arrange(.data$runid, .data$variable, .data$category)
}
