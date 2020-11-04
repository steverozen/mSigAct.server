# The code below are adapted from SignatureEstimation R package
# Original code source: https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz
# Original paper: https://academic.oup.com/bioinformatics/article-pdf/34/2/330/25114255/btx604.pdf
# Reference manual for SignatureEstimation R package:https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.pdf


#' decomposeQP Function
#'
#' This function allows to get the optimal solution by using dual method to solve the quadratic programming problem.
#' @param m observed tumor profile vector for a single patient/sample, 96 by 1. m is normalized.
#' @param P signature profile matrix, 96 by N(N = # signatures, COSMIC: N=30)
#' @param ... some control parameter that can be passed into the solve.QP function
#' @keywords quadratic programming
#' @export
decomposeQP <- function(m, P, ...){
  # N: how many signatures are selected
  N = ncol(P)
  # G: matrix appearing in the quatric programming objective function
  G = t(P) %*% P
  # C: matrix constraints under which we want to minimize the quatric programming objective function.
  C <- cbind(rep(1,N), diag(N))
  # b: vector containing the values of b_0.
  b <- c(1,rep(0,N))
  # d: vector appearing in the quatric programming objective function
  d <- t(m) %*% P
  
  #Solve quadratic programming problem
  out = quadprog::solve.QP(Dmat = G, dvec = d, Amat = C, bvec = b, meq = 1)
  
  #Some exposure values are negative, but very close to 0
  #Change these neagtive values to zero and renormalized
  exposures = out$solution
  exposures[exposures < 0] = 0
  exposures = exposures/sum(exposures)
  
  # return the exposures
  return(exposures)
}

#' decomposeSA Function
#'
#' This function allows to get the optimal solution by using simulated annealing to solve the optimization problem.
#' @param m observed tumor profile vector for a single patient/sample, 96 by 1. m is normalized.
#' @param P signature profile matrix, 96 by N(N = # signatures, COSMIC: N=30)
#' @param control some control parameter that can be passed into the GenSA function
#' @keywords simulated annealing
#' @export
#Wrapper for GenSA function to run simulated annealing
decomposeSA <- function(m, P, control = list()) {
  #objective function to be minimized
  #local version of Frobenius norm to simplify and speed-up the objective function
  FrobeniusNorm.local <- function(exposures) {
    estimate = P %*% exposures
    return(sqrt(sum((m - (estimate / sum(estimate)))^2)))
  }
  # N: how many signatures are selected
  N = ncol(P)
  #change our suggestion to control GenSA function based on user's requirements
  our.control = list(maxit=1000, temperature=10, nb.stop.improvement=1000, simple.function=TRUE)
  our.control[names(control)] = control
  #Solve the problem using simulated annealing package GenSA
  sa = GenSA::GenSA(lower=rep(0.0,N), upper=rep(1.0,N), fn=FrobeniusNorm.local, control=our.control)
  #Normalize the solution
  exposures = sa$par/sum(sa$par)
  
  return(exposures)
}


#' findSigExposures Function
#' wrapper function
#' This function allows to obtain the optimal solution by specifying quadratic programming or simulated annealing to solve the optimization problem.
#' @param M observed tumor profile matrix for all the patient/sample, 96 by G. G is the number of patients. Each column can be mutation
#' counts, or mutation probabilities. Each column will be normalized to sum up to 1.
#' @param P signature profile matrix, 96 by N(N = # signatures, COSMIC: N=30)
#' @param decompostion.method which method is selected to get the optimal solution: decomposeQP or decomposeSA
#' @keywords optimal method: QP or SA
#' @keywords internal
findSigExposures <- function(M, P, decomposition.method = decomposeQP, ...) {
  ## process and check function parameters
  ## M, P
  M = as.matrix(M)
  P = as.matrix(P)
  if(nrow(M) != nrow(P))
    stop("Matrices 'M' and 'P' must have the same number of rows (mutations types).")
  if(any(rownames(M) != rownames(P)))
    stop("Matrices 'M' and 'P' must have the same row names (mutations types).")
  if(ncol(P) == 1)
    stop("Matrices 'P' must have at least 2 columns (signatures).")
  ## decomposition.method
  if(!is.function(decomposition.method))
    stop("Parameter 'decomposition.method' must be a function.")
  ## normalize M by column (just in case it is not normalized)
  M = apply(M, 2, function(x){x/sum(x)})
  
  ## find solutions
  ## matrix of signature exposures per sample/patient (column)
  exposures = apply(M, 2, decomposition.method, P, ...)
  rownames(exposures) = colnames(P)
  colnames(exposures) = colnames(M)
  
  ## compute estimation error for each sample/patient (Frobenius norm)
  errors = sapply(seq(ncol(M)), function(i) FrobeniusNorm(M[,i],P,exposures[,i]))
  names(errors) = colnames(M)
  
  return(list(exposures=exposures, errors=errors))
}

#' bootstrapSigExposures Function
#' This function allows to obtain the bootstrap distribution of the signature exposures of a certain tumor sample
#' @param m observed tumor profile vector for a patient/sample, 96 by 1. It can be mutation
#' counts, or mutation probabilities.
#' @param P signature profile matrix, 96 by N (N = # signatures, COSMIC: N=30).
#' @param mutation.count if m is a vector of counts, then mutation.count equals the summation of all the counts.
#'  If m is probabilities, then mutation.count has to be specified.
#' @param R The number of bootstrap replicates.
#' @param decompostion.method which method is selected to get the optimal solution: decomposeQP or decomposeSA
#' @keywords bootstrap
#' @keywords internal
bootstrapSigExposures <- function(m, P, R, mutation.count = NULL, decomposition.method = decomposeQP, ...) {
  ## process and check function parameters
  ## m, P
  P = as.matrix(P)
  if(length(m) != nrow(P))
    stop("Length of vector 'm' and number of rows of matrix 'P' must be the same.")
  if(any(names(m) != rownames(P)))
    stop("Elements of vector 'm' and rows of matrix 'P' must have the same names (mutations types).")
  if(ncol(P) == 1)
    stop("Matrices 'P' must have at least 2 columns (signatures).")
  
  ## if 'mutation.count' is not specified the 'm' has to contain counts
  if(is.null(mutation.count)) {
    if(all(is.wholenumber(m))) # alternative test might be (assuming only 2 possibilities): sum(m) > 1
      mutation.count = sum(m)
    else
      stop("Please specify the parameter 'mutation.count' in the function call or provide mutation counts in parameter 'm'.")
  }
  
  ## normalize m to be a vector of probabilities.
  m = m / sum(m)
  
  ## find optimal solutions using using provided decomposition method for each bootstrap replicate
  ## matrix of signature exposures per replicate (column)
  K = length(m) #number of mutation types
  exposures = replicate(R, {
    mutations_sampled = sample(seq(K), mutation.count, replace = TRUE, prob = m)
    m_sampled = as.numeric(table(factor(mutations_sampled, levels=seq(K))))
    m_sampled = m_sampled / sum(m_sampled)
    decomposition.method(m_sampled, P, ...)
  })
  rownames(exposures) = colnames(P)
  colnames(exposures) = paste0('Replicate_',seq(R))
  
  ## compute estimation error for each replicate/trial (Frobenius norm)
  errors = apply(exposures, 2, function(e) FrobeniusNorm(m,P,e))
  names(errors) = colnames(exposures)
  
  return(list(exposures=exposures, errors=errors))
}

#' suboptimalSigExposures Function
#' This function allows to obtain the simulated annealing distribution of the signature exposures of a certain tumor sample
#' @param m observed tumor profile vector for a patient/sample, 96 by 1. It can be mutation
#' counts, or mutation probabilities.
#' @param P signature profile matrix, 96 by N(N = # signatures, COSMIC: N=30)
#' @param mutation.count if m is a vector of counts, then mutation.count equals the summation of all the counts.
#'  If m is probabilities, then mutation.count has to be specified.
#' @param R The number of replicates/trials of simulated annealing.
#' @param optimal.error if it is NULL, then use SA method to obtain. Or you can provide one.
#' @param suboptimal.factor suboptimal error.
#' @param control some control parameter that can be passed into the function
#' @keywords suboptimal 'simulated annealing'
#' @keywords internal
suboptimalSigExposures <- function(m, P, R, optimal.error = NULL, suboptimal.factor = 1.05, control = list()) {
  ## process and check function parameters
  ## m, P
  P = as.matrix(P)
  if(length(m) != nrow(P))
    stop("Length of vector 'm' and number of rows of matrix 'P' must be the same.")
  if(any(names(m) != rownames(P)))
    stop("Elements of vector 'm' and rows of matrix 'P' must have the same names (mutations types).")
  if(ncol(P) == 1)
    stop("Matrices 'P' must have at least 2 columns (signatures).")
  
  ## normalize m to be a vector of probabilities.
  m = m/sum(m)
  
  ## If optimal.error is null, we compute it by default SA method ('decomposeSA').
  if(is.null(optimal.error)) {
    sa_expos = decomposeSA(m, P, control)
    optimal.error = FrobeniusNorm(m, P, sa_expos)
  }
  
  ## find suboptimal solutions using simulated annealing with predefined error (with respect to optimal error)
  control$threshold.stop = suboptimal.factor * optimal.error
  ## matrix of signature exposures per replicate/trial (column)
  exposures = replicate(R, decomposeSA(m, P, control))
  rownames(exposures) = colnames(P)
  colnames(exposures) = paste0('Replicate_',seq(R))
  
  ## compute estimation error for each replicate/trial (Frobenius norm)
  errors = apply(exposures, 2, function(e) FrobeniusNorm(m,P,e))
  names(errors) = colnames(exposures)
  
  return(list(exposures=exposures, errors=errors))
}

FrobeniusNorm <- function(M, P, E) {
  sqrt(sum((M-P%*%E)^2))
}

is.wholenumber <- function(x, tol = .Machine$double.eps) {
  abs(x - round(x)) < tol
}
