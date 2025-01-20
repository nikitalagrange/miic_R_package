#' Compute complexity NML non Markov equivalent between two variables with
#' conditioning set
#' @description
#' The complexity term is calculated knowing that X is the parent of Y
#' such that X -> Y. So the conditioning set corresponds to the parents of node
#' Y
#' @param X [vector]
#' A vector that contains the observational data of the first variable.
#' @param Y [vector]
#' A vector that contains the observational data of the second variable.
#' @param df_conditioning [data frame]
#' The data frame of the observations of the conditioning variables. Parents of
#' Y node
#' @return Numeric value of the complexity NML term
#' @export
#' @useDynLib miic
#' @importFrom stats density sd
#'
#' @examples
#' N <- 1000
#' X <- sample(x = c(1, 2, 3, 4, 5), size = N, replace = TRUE)
#' Y <- sample(x = c(1, 2, 3), size = N, replace = TRUE)
#' Z <- sample(x = c(5, 6), size = N, replace = TRUE)
#' res <- computeComplexityOrient(as.factor(X), as.factor(Y), df_conditioning = as.factor(Z))
#' message("k(X;Y|Z) = ", res)
computeComplexityOrient <- function(X, Y, df_conditioning = NULL) {
  input_data <- data.frame(X, Y)
  if (!is.null(df_conditioning)) {
    input_data <- data.frame(input_data, df_conditioning)
  }

  n_samples <- nrow(input_data)
  n_nodes <- ncol(input_data)

  is_continuous <- sapply(input_data, is.numeric)

  # Numeric factor matrix, level starts from 0
  input_factor <- as.matrix(apply(
    input_data, 2,
    function(x) (as.numeric(factor(x, levels = unique(x))) - 1)
  ))
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  # Data list, numeric for continuous columns, -1 for discrete columns
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  # Order list, order(column) for continuous columns (index starting from 0),
  # -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1:n_nodes)) {
    input_double[, i] <- rep_len(-1, n_samples)
    input_order[, i] <- rep_len(-1, n_samples)
  }

  arg_list <- list(
    "is_continuous" = is_continuous,
    "levels" = max_level_list,
    "n_eff" = -1,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples
  )

  cpp_input <- list(
    "factor" = as.vector(input_factor),
    "double" = as.vector(input_double),
    "order" = as.vector(input_order)
  )
  # Call cpp code
  rescpp <- complexOrient(cpp_input, arg_list)
  return(rescpp)
}


#' Compute regret term
#' @description
#' The logarithm of the denominator of Normalized Maximum Likelihood (NML) is
#' often called the regret
#' @references
#' \itemize{
#' \item Silander et al., \emph{AISTATS 2018.}
#'   https://proceedings.mlr.press/v84/silander18a/silander18a.pdf}
#' @param N [an integer]
#' Length of the categotical data
#' @param r [an integer]
#' Number of different categories (variable level) data
#' @param data [a data frame]
#' A data frame of the same length N as the variable studied
#' @return Numeric value of the regret term
#' @export
#' @useDynLib miic
#' @importFrom stats density sd
#'
#' @examples
#' library(miic)
#' N <- 1000
#' r <- 45
#' X <- sample(x = c(1, 2, 3, 4, 5), size = N, replace = TRUE)
#' res <- regret(N, r, data.frame(X))
#' message("reg(N, r) = ", res)
#'
regret <- function(N, r, data) {
  n_samples <- nrow(data)
  n_nodes <- ncol(data)
  input_factor <- apply(data, 2, function(x) {
    (as.numeric(factor(x, levels = unique(x))) - 1)
  })
  input_factor[is.na(input_factor)] <- -1
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  input_factor <- as.vector(as.matrix(input_factor))
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  # Order list, order(column) for continuous columns (index starting from 0, NA
  # mapped to -1), -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  is_continuous <- rep(0, n_nodes)
  for (i in c(1:ncol(data))) {
    input_double[, i] <- rep_len(-1, n_samples)
    input_order[, i] <- rep_len(-1, n_samples)
  }
  input_order <- as.vector(input_order)
  input_double <- as.vector(input_double)
  arg_list <- list(
    "is_continuous" = as.numeric(is_continuous),
    "levels" = max_level_list,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples
  )
  cpp_input <- list(
    "factor" = input_factor, "double" = input_double,
    "order" = input_order
  )
  res <- regret_term(N, r, cpp_input, arg_list)
  return(res)
}
