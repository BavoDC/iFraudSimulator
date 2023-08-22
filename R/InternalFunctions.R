#' Internal function to build the adjacency matrix.
#'
#' @param Network Network dataset.
#' @param decay Decay factor.
#' @param Today Date.
#'
#' @return An adjacency/weight matrix
.AdjMat <- function(Network, decay = 0, Today = Sys.Date()) {
  aMat =
    if (decay != 0) {
      w        = exp(-decay * as.numeric(difftime(Today, Network$Date, units = "weeks")))
      w[w > 1] = 1
      sparseMatrix(i = Network$startNode,
                   j = Network$endNode,
                   x = w)
    } else {
      sparseMatrix(i = Network$startNode,
                   j = Network$endNode,
                   x = 1)
    }
  return(aMat)
}

#' Internal function to set up symmetrically normalized matrix
#'
#' @param aMat Adjacency matrix.
#'
#' @return symmetrically normalized matrix
.SNMM   <- function(aMat) {
  Dp = sparseMatrix(1:nrow(aMat), 1:nrow(aMat), x = 1 / sqrt(rowSums(aMat)))
  Dc = sparseMatrix(1:ncol(aMat), 1:ncol(aMat), x = 1 / sqrt(colSums(aMat)))
  S  = t(Dp %*% aMat %*% Dc)
  return(S)
}

#' Internal function to set up query vector
#'
#' @param fraudMat Data frame with binary indicator for fraudulent claims
#' @param decay Decay parameter.
#' @param Today Date
#'
#' @return Query vector
.QueryVector <- function(fraudMat, decay, Today) {
  QueryValueOne = F
  c0 =
    if(decay != 0) {
      Tmp = exp(-decay * as.numeric(difftime(Today, fraudMat$Date, units = "weeks"))) * fraudMat$FraudInd
      Tmp =
        if(sum(Tmp) > 1e-14 & !QueryValueOne) {
          Tmp / sum(Tmp)
        } else {
          Tmp
        }
    } else {
      if(!QueryValueOne) {
        fraudMat$FraudInd / sum(fraudMat$FraudInd)
      } else {
        fraudMat$FraudInd
      }
    }
  return(c0)
}
