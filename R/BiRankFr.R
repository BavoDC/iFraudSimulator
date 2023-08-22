#' BiRank algorithm
#'
#' BiRank algorithm to compute fraud scores.
#'
#' @param sNetwork Data frame containing the variables startNode, endNode and possibly Date.
#' @param fraudMat Data frame containing the variables FraudInd (binary variable indicating fraudulent claims) and possibly Date.
#' @param Today Date of analysis (e.g. 1/1/2020), default is \code{Sys.Date()}. Supply either as object with class \code{Date}
#' or as character string in the format \%d/\%m/\%Y.
#' @param decayR Parameter for exponential decay of recency of relation (in weeks).
#' @param decayF  Parameter for exponential decay of recency of fraud (in weeks).
#' @param alpha Damping factor for propagation algorithm (return to start).
#' @param maxiter Maximum number of iterations for propagation algorithm.
#' @param Epsilon Positive convergence tolerance \eqn{\epsilon}.
#' @param printProgress Logical, indicates whether progress of the algorithm has to be printed.
#' @param pInit Initial values for the party score vector \bold{\eqn{p}}.
#' @param cInit Initial values for the fraud score vector \bold{\eqn{c}}.
#' @param ConvCriterion Which convergence criterion to use. \code{"Sep"} uses
#' \eqn{||p - p_{old}||_{2} / ||p_{old}||_{2} < \epsilon} and \eqn{||c - c_{old}||_{2} / ||c_{old}||_{2} < \epsilon}.
#' \code{"Whole"} uses \eqn{||x - x_{old}||_{2} / ||x_{old}||_{2} < \epsilon} with \code{x = c(c, p)}. \code{"Order"} uses
#' the same convergence criterion as \code{"Sep"} and checks if the order of the elements of \eqn{\bold{p}} and \eqn{\bold{c}}
#' has not changed since the previous iteration.
#'
#' @return A list with the following components:
#' @return \item{ResultsClaims}{A data frame containing the claim IDs, fraud scores, scaled and normalized fraud scores.}
#' @return \item{ResultsParties}{A data frame containing the party IDs, party scores, scaled and normalized party scores.}
#' @return \item{AdjacencyMatrix}{Adjancency/weight matrix indicating which nodes are connected.}
#' @return \item{iter}{Number of iterations that the algorithm needed to converge.}
#' @return \item{Converged}{Logical, indicating whether the algorithm converged.}
#'
#'
#' @examples
#' library(iFraudSimulator)
#' NetwLabel = data.frame(
#' startNode = c('P2', 'P3', 'P3', 'C1', 'C1', 'C5', 'P1', 'P4', 'C2', 'C2'),
#' endNode = c('C3', 'C3', 'C4', 'P2', 'P3', 'P3', 'C1', 'C5', 'P1', 'P4'),
#' stringsAsFactors = FALSE
#' )
#' NetwLabel[grepl("C", NetwLabel$startNode), 1:2] = NetwLabel[grepl("C", NetwLabel$startNode), 2:1]
#' NetwLabel = NetwLabel[order(NetwLabel$startNode), ]
#' NetwLabel$startNode = as.numeric(gsub("P", "", NetwLabel$startNode))
#' NetwLabel$endNode   = as.numeric(gsub("C", "", NetwLabel$endNode))
#' NetwLabel
#'
#' c0      = c(rep(0, 3), 1, 0)
#' Results = BiRankFr(NetwLabel, data.frame(FraudInd = c0))
#' Results
BiRankFr <- function(sNetwork, fraudMat, Today = Sys.Date(), decayR = 0, decayF = 0,
                     alpha = 0.85, maxiter = 1e3, Epsilon = 1e-14,
                     printProgress = FALSE, pInit = NULL, cInit = NULL, ConvCriterion = c("Sep", "Whole", "Order")) {
  if(!is.null(pInit) | !is.null(cInit)) {
    if(any(!sapply(c(pInit, cInit), is.vector)))
      stop("Provide vectors for pInit/cInit.")
  }
  ConvCriterion = match.arg(ConvCriterion)
  ConvCriterion =
    if(ConvCriterion == "Sep") {
      function(c, p, cOld, pOld, Epsilon = Epsilon) {
        as.vector(sqrt(crossprod(c - cOld)) / sqrt(crossprod(cOld))) < Epsilon &
        as.vector(sqrt(crossprod(p - pOld)) / sqrt(crossprod(pOld))) < Epsilon
      }
    } else if(ConvCriterion == "Whole") {
      function(c, p, cOld, pOld, Epsilon = Epsilon) {
        xi    = as.vector(rbind(c, p))
        xiOld = as.vector(rbind(cOld, pOld))
        as.vector(sqrt(crossprod(xi - xiOld)) / sqrt(crossprod(xiOld))) < Epsilon
      }
    } else {
      function(c, p, cOld, pOld, Epsilon = Epsilon) {
        as.vector(sqrt(crossprod(c - cOld)) / sqrt(crossprod(cOld))) < Epsilon &
        as.vector(sqrt(crossprod(p - pOld)) / sqrt(crossprod(pOld))) < Epsilon &
        crossprod(order(c) - order(cOld)) == 0 & crossprod(order(p) - order(pOld)) == 0
      }
    }

  #### 1. Function ####

  #### 1. Preparations ####
  sNetwork$Date =
    if (decayR != 0 & !is.Date(sNetwork$Date)) {
      as.Date(sNetwork$Date, format = "%d/%m/%Y")
    } else if(decayR == 0) {
      1
    }
  if(any(colnames(fraudMat) == "Date")) {
    fraudMat$Date =
      if (decayF != 0 & !is.Date(fraudMat$Date)) {
        as.Date(fraudMat$Date, format = "%d/%m/%Y")
      } else if (decayF == 0) {
        NULL
      }
  }

  if(printProgress)
    cat("\n\nSetting up adjacency matrix.\n\n")
  aMat = .AdjMat(sNetwork, decay = decayR, Today = Today)

  if(printProgress)
    cat("\n\nNormalizing matrix.\n\n")
  S    = .SNMM(aMat)

  if(printProgress)
    cat("\n\nInitiating query vector.\n\n")
  c0   = .QueryVector(fraudMat, decay = decayF, Today = Today)

  if(!is.null(pInit) | !is.null(cInit)) {
    if(length(pInit) != ncol(S) | length(cInit) != nrow(S))
      stop("Wrong length of pInit/cInit vectors")
  }

  pOld  = if(!is.null(pInit)) pInit else as.vector(runif(ncol(S)))
  cOld  = if(!is.null(cInit)) cInit else as.vector(runif(nrow(S)))
  Conv  = F
  iter  = 1

  if(printProgress)
    cat("\n\nRunning algorithm.\n\n")
  while(!Conv) {
    c = alpha * (S %*% pOld) + (1 - alpha) * c0
    p = t(t(c) %*% S)
    if(ConvCriterion(c, p, cOld, pOld, Epsilon))
      break
    cOld = c
    pOld = p
    iter = iter + 1
    if(iter > maxiter) {
      warning("Maximum number of iterations has been reached.", immediate. = T)
      break
    }
    if(iter %% 1e2 == 0 & printProgress)
      cat("\n\nIteration number", iter, "\n\n")
  }

  c = as.vector(c)
  p = as.vector(p)

  ResultsClaims  = cbind.data.frame(ID          = seq_len(ncol(aMat)),
                                    Score       = c,
                                    StdScore    = scale(c),
                                    ScaledScore = (c - min(c)) / diff(range(c)))
  ResultsParties = cbind.data.frame(ID          = seq_len(nrow(aMat)),
                                    Score       = p,
                                    StdScore    = scale(p),
                                    ScaledScore = (p - min(p)) / diff(range(p)))
  ResultsClaims  = ResultsClaims[order(ResultsClaims$Score, decreasing = T), ]
  ResultsParties = ResultsParties[order(ResultsParties$Score, decreasing = T), ]

  Results = list(ResultsClaims   = ResultsClaims,
                 ResultsParties  = ResultsParties,
                 AdjacencyMatrix = aMat,
                 iter = iter,
                 Converged = iter < maxiter)
  class(Results) = "BiRankFr"
  return(Results)
}


