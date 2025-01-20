#' Compute edge score (mutual information) with conditioning set
#' @description
#' Compute edge score (mutual information) with conditioning set and Markov
#' equivalent (symmetric) complexity NML term
#' @param X [vector]
#' A vector that contains the observational data of the first variable.
#' @param Y [vector]
#' A vector that contains the observational data of the second variable.
#' @param df_conditioning [data frame]
#' The data frame of the observations of the conditioning variables.
#' @return A list that contains :
#' \itemize{
#' \item info: The estimation of (conditional) mutual information without the
#' complexity cost.
#' \item infok: The estimation of (conditional) mutual information with the
#' complexity cost (\eqn{Ik = -I + cplx}).}
#' @export
#' @useDynLib miic
#' @importFrom stats density sd
#' @examples
#' N <- 1000
#' X <- sample(x = c(1, 2, 3, 4, 5), size = N, replace = TRUE)
#' Y <- sample(x = c(1, 2, 3), size = N, replace = TRUE)
#' res <- computeEdgeScore(as.factor(X), as.factor(Y))
#' message("I(X;Y) = ", res$info)
#' message("Ik(X;Y) = ", res$infok)
#' Z <- sample(x = c(5, 6), size = N, replace = TRUE)
#' res <- computeEdgeScore(as.factor(X), as.factor(Y), df_conditioning = as.factor(Z))
#' message("I(X;Y|Z) = ", res$info)
#' message("Ik(X;Y|Z) = ", res$infok)

computeEdgeScore <- function(X, Y, df_conditioning = NULL) {
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
  rescpp <- EdgeScore(cpp_input, arg_list)
  result <- list()
  result$info <- rescpp$info
  result$infok <- rescpp$infok
  return(result)
}
