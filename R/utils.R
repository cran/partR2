#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' Captures and suppresses (still to find out why) warnings of an expression
#'
#' This function is used within partR2 to capture lme4 model fitting warnings in the
#' bootstrap and permutation procedures.
#'
#' @param expr An expression, such as the sequence of code used by rptR to calculate
#' bootstrap or permutation estimates
#' @return List of warnings.
#' @keywords internal

with_warnings <- function(expr) {
  myWarnings <- NULL
  myMessages <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(warnings = myWarnings)
}

#' Create list of combination of variables.
#'
#' @inheritParams partR2
#'
#' @keywords internal
#'
#' @return list with all combinations of predictors specified in partvars/partbatch
#'
make_combs <- function(partvars, partbatch, max_level) {
  # create list of all unique combinations except for the full model
  if (!is.null(partvars)) {
    if (length(partvars) > 1) {
      all_comb <- unlist(lapply(
        1:(length(partvars)),
        function(x) utils::combn(partvars, x, simplify = FALSE)
      ),
      recursive = FALSE
      )
    } else if (length(partvars) == 1) {
      all_comb <- as.list(partvars)
    }
  } else {
    all_comb <- NA
  }

  # check batches
  if (!is.null(partbatch)) {
    if (!is.list(partbatch)) stop("partbatch must be a list")
    if (is.null(all_comb)) all_comb <- NA
    # first, make combinations of partbatches
    combs_num <- purrr::map(
      1:length(partbatch),
      function(m) {
        utils::combn(length(partbatch), m,
          simplify = FALSE
        )
      }
    ) %>%
      unlist(recursive = FALSE)
    comb_batches <- purrr::map(combs_num, function(x) unlist(partbatch[x]))
    # now add those to partvar combs
    comb_batches2 <- purrr::map(comb_batches, function(x) {
      purrr::map(all_comb, function(z) c(z, x))
    }) %>%
      unlist(recursive = FALSE)
    # now add those to all_combs
    all_comb <- c(partbatch, all_comb, comb_batches2)
    # check for duplicates or NA and remove in case
    all_comb <- purrr::map(all_comb, function(x) x[!(duplicated(x) | is.na(x))])
    # last step remove any empty list elements
    all_comb[purrr::map(all_comb, length) == 0] <- NULL
    # remove potential duplicates
    all_comb <- all_comb[!(duplicated(purrr::map(all_comb, function(x) as.character(sort(x)))))]
    # change all formats to unnamed character vecotr
    all_comb <- purrr::map(all_comb, function(x) as.character(unname(x)))
  }

  # commonality coefficients up to max_level (e.g. 3 for
  # the cc of 3 predictors)
  if (!is.null(partbatch) & !is.null(max_level)) {
    stop("Argument max_level does currently not work in combination with argument partbatch,
                please use partvars or leave max_level at NULL")
  }
  if (!is.null(max_level)) {
    remove_combs <- purrr::map_lgl(all_comb, function(x) length(x) > max_level)
    all_comb[remove_combs] <- NULL
  }

  all_comb
}


#' List to data.frame with bootstrap samples per row
#'
#'
#' @param lcol A bootstrap list-column
#' @return data.frame with bootstrap estimates as list-column
#' @keywords internal

boot_to_df <- function(lcol) {
  out <- lcol %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$term) %>%
    dplyr::summarise(boot_estimates = list(.data$estimate))
  out
}


#' Modify term names if partbatch is a named list
#'
#'
#' @param partbatch list with batches
#' @param part_names character vector with names of terms (and their combinations)
#' @return Modified names for partvars, where partvar combinations are replaced
#' with partbatch names.
#' @keywords internal
#'
# change partbatch with names, if present
mod_names_partbatch <- function(partbatch, part_names) {
  added_batches <- purrr::map(partbatch, function(x) {
    out <- paste(x, collapse = "+")
    out
  })

  for (i in 1:length(added_batches)) {
    if (is.null(names(added_batches)[i])) next()
    part_names <- gsub(added_batches[[i]], names(added_batches)[i], part_names, fixed = TRUE)
  }

  part_names
}

#' Adds columns for CI
#' @param df data.frame with point estimates
#' @return data.frame with two additional columns for lower and upper CI containing NAs
#' @keywords internal
add_CI_cols <- function(df) {
  df %>% dplyr::mutate(CI_lower = NA, CI_upper = NA)
}
