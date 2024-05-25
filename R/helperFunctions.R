#' Min-max scaling to normalize a vector
#'
#' @param x an object of type vector.
#' @param Type a character vector indicating whic type of normalization has to be used. \code{Type = "1"} scales
#' the vector into the range [-1, 1] and \code{Type = "2"} scales the vector into the range [0, 1].
#'
#' @return the normalized vector
#' @export
#'
#' @examples
#' x <- rnorm(10)
#' normalize(x)
normalize <- function(x, Type = c("1", "2")) {
  Type = match.arg(Type)
  if(!is.vector(x))
    stop("Has to be of type vector.")
  if(anyNA(x))
    x = x[!is.na(x)]
  if(Type == 2)
    (x - min(x)) / (max(x) - min(x))
  else
    2 * (x - min(x)) / (max(x) - min(x)) - 1
}

sumNA <- function(x) sum(is.na(x))

AdjMatNeighbors <- function(aMat, Order = 2) {
  if(Order < 0)
    stop("Not possible")
  else if(Order == 0)
    return(Matrix(0, nrow(aMat), ncol(aMat)))
  else if(Order == 1)
    return(aMat)

  if(!Matrix::isSymmetric(aMat))
    stop("Only symmetric matrices allowed.")

  AllMat = list(aMat)

  for(i in seq_len(Order)[-1]) {
    N = AllMat[[i - 1]] %*% aMat
    diag(N)  = 0
    if(is(N, "sparseMatrix"))
      N@x[N@x > 1] = 1
    else
      N[N > 1] = 1
    AllMat[[i]] = N
  }
  for(i in seq_along(AllMat)[-length(AllMat)])
    N = N - AllMat[[i]]
  if(is(N, "sparseMatrix"))
    N@x[N@x < 0] = 0
  else
    N[N < 0] = 0
  return(N)
}

AddInfoNetwork <- function() {
  N1 <- f <- AdjMat <- startNode <- ClaimID <- Party <- DtClaim <- Fraud <- i <- endNode <- NULL
  SumN1    = summary(N1)
  N1scores = tapply(SumN1$x, SumN1$j, function(x) do.call("cbind", lapply(f, function(f) f(x))))
  N1scores = do.call("rbind", N1scores)
  colnames(N1scores) = c("n1.q1", "n1.med", "n1.max", "n1.size")
  rownames(N1scores) = seq_len(ncol(N1))

  IDN2Claims   = unique(AdjMat[startNode %in% AdjMat[ClaimID %in% Dt$ClaimID][["startNode"]]][["ClaimID"]])
  AdjMatSubset = AdjMat[ClaimID %in% IDN2Claims]
  setorder(AdjMatSubset, ClaimID, Party)
  AdjMatSubset  = AdjMatSubset[, `:=` (
    startNode = as.numeric(as.factor(Party)),
    endNode   = as.numeric(as.factor(ClaimID))
  )]
  N1P   = tapply(AdjMatSubset$endNode, AdjMatSubset$startNode, unique)
  N1Cl  = tapply(AdjMatSubset$startNode, AdjMatSubset$endNode, unique)
  FrCl  = unique(AdjMatSubset[ClaimID %in% DtClaim[Fraud == T][["ClaimID"]]][["endNode"]])
  NFrCl = unique(AdjMatSubset[ClaimID %in% DtClaim[Fraud == F][["ClaimID"]]][["endNode"]])

  AdjMatSubsetOrig = AdjMatSubset
  W       = sparseMatrix(i = AdjMatSubset$endNode,
                         j = AdjMatSubset$startNode,
                         x = 1)
  AdjMatSubset = rbind(cbind(matrix(0, nrow(W), nrow(W)), W),
                       cbind(t(W), matrix(0, ncol(W), ncol(W))))
  NCl   = nrow(W)

  N2scores = foreach(i = seq_len(NCl)) %do%
    {
      x    = N1Cl[[i]]
      N2C  = unique(do.call("c", N1P[x]))
      N2C  = N2C[!N2C %in% i]
      if(length(N2C) != 0 ) {
        NFr  = sum(N2C %in% FrCl)
        NnFr = sum(N2C %in% NFrCl)
        n    = AdjMatSubsetOrig[endNode == i][1, ][["n2Size"]]
        Res = c(
          ratioFr = if (NFr == 0) 0 else NFr / n,
          ratioNonFr = if (NnFr == 0) 0 else NnFr / n,
          binFraud = as.numeric(NFr != 0)
        )
      } else {
        Res = rep(0, 3)
      }
      return(Res)
    }
  N2scores = data.table(do.call("rbind", N2scores))
  colnames(N2scores) = c("n2.ratioFraud", "n2.ratioNonFraud", "n2.binFraud")
  N2scores %>% class
  N2scores[, endNode := seq_len(nrow(N2scores))]
  N2scores[, ClaimID := AdjMatSubsetOrig$ClaimID[match(endNode, AdjMatSubsetOrig$endNode)]]
  N2scores[, endNode := NULL]

  Dt = merge(Dt, N2scores, by = "ClaimID")
  return(Dt)
}

AddInfoN2 <- function(Dt, AdjMat, Parallel = Parallel, std = stdize, DtCl = DtClaim) {
  stdize <- DtClaim <- startNode <- ClaimID <- Party <- Fraud <- i <- endNode <- NULL
  IDN2Claims   = unique(AdjMat[startNode %in% AdjMat[ClaimID %in% Dt$ClaimID][["startNode"]]][["ClaimID"]])
  AdjMatSubset = AdjMat[ClaimID %in% IDN2Claims]
  setorder(AdjMatSubset, ClaimID, Party)
  AdjMatSubset  = AdjMatSubset[, `:=` (
    startNode = as.numeric(as.factor(Party)),
    endNode   = as.numeric(as.factor(ClaimID))
  )]
  N1P   = tapply(AdjMatSubset$endNode, AdjMatSubset$startNode, unique)
  N1Cl  = tapply(AdjMatSubset$startNode, AdjMatSubset$endNode, unique)
  FrCl  = unique(AdjMatSubset[ClaimID %in% DtCl[Fraud == T][["ClaimID"]]][["endNode"]])
  NFrCl = unique(AdjMatSubset[ClaimID %in% DtCl[Fraud == F][["ClaimID"]]][["endNode"]])

  AdjMatSubsetOrig = AdjMatSubset
  W       = sparseMatrix(i = AdjMatSubset$endNode,
                         j = AdjMatSubset$startNode,
                         x = 1)
  AdjMatSubset = rbind(cbind(Matrix(0, nrow(W), nrow(W)), W),
                       cbind(t(W), Matrix(0, ncol(W), ncol(W))))
  NCl   = nrow(W)

  N2scores =
    if(Parallel) {
      foreach(i = seq_len(NCl), .packages = "data.table") %dopar%
        {
          x    = N1Cl[[i]]
          N2C  = unique(do.call("c", N1P[x]))
          N2C  = N2C[!N2C %in% i]
          if(length(N2C) != 0 ) {
            NFr  = sum(N2C %in% FrCl)
            NnFr = sum(N2C %in% NFrCl)
            n    = AdjMatSubsetOrig[endNode == i][1, ][["n2Size"]]
            Res = c(
              ratioFr = if (NFr == 0) 0 else NFr / n,
              ratioNonFr = if (NnFr == 0) 0 else NnFr / n,
              binFraud = as.numeric(NFr != 0)
            )
          } else {
            Res = rep(0, 3)
          }
          return(Res)
        }
    } else {
      foreach(i = seq_len(NCl)) %do%
        {
          x    = N1Cl[[i]]
          N2C  = unique(do.call("c", N1P[x]))
          N2C  = N2C[!N2C %in% i]
          if(length(N2C) != 0 ) {
            NFr  = sum(N2C %in% FrCl)
            NnFr = sum(N2C %in% NFrCl)
            n    = AdjMatSubsetOrig[endNode == i][1, ][["n2Size"]]
            Res = c(
              ratioFr = if (NFr == 0) 0 else NFr / n,
              ratioNonFr = if (NnFr == 0) 0 else NnFr / n,
              binFraud = as.numeric(NFr != 0)
            )
          } else {
            Res = rep(0, 3)
          }
          return(Res)
        }
    }
  N2scores = data.table(do.call("rbind", N2scores))
  colnames(N2scores) = c("n2.ratioFraud", "n2.ratioNonFraud", "n2.binFraud")
  stdVars = c("n2.ratioFraud", "n2.ratioNonFraud")
  N2scores[, (stdVars) := lapply(.SD, std), .SDcols = stdVars]
  N2scores[, endNode := seq_len(nrow(N2scores))]
  N2scores[, ClaimID := AdjMatSubsetOrig$ClaimID[match(endNode, AdjMatSubsetOrig$endNode)]]
  N2scores[, endNode := NULL]

  Dt = merge(Dt, N2scores, by = "ClaimID")
  return(Dt)
}

AddInfo <- function(DtClaim, DtSubset, SelectBatch, AdjMatOrig) {
  Fraud <- ClaimID <- startNode <- NULL
  Cols     = colnames(DtSubset)[!colnames(DtSubset) %in% c("Fraud", colnames(DtSubset)[grepl("n2.", colnames(DtSubset))])]
  PHIDs    = unique(DtClaim$IDPH)    %>% .[!. %in% SelectBatch]
  ClaimIDs = unique(DtClaim$ClaimID) %>% .[!. %in% DtSubset[["ClaimID"]]]
  FraudIDs = DtSubset[Fraud == 1,][["ClaimID"]]
  setkeyv(DtClaim, Cols)
  setkeyv(DtSubset, Cols)
  DtClaim[DtClaim[DtSubset[, .SD, .SDcols = Cols], nomatch = 0L, which = T], Fraud := DtSubset$Fraud]
  Neighbors = AdjMatOrig[ClaimID %in% FraudIDs, ][["startNode"]]
  n2Claims  = AdjMatOrig[startNode %in% Neighbors & ClaimID %in% ClaimIDs][["ClaimID"]]
  if(length(n2Claims) > 0.1 * NrUnique(DtClaim$ClaimID))
    n2Claims = sample(n2Claims, floor(0.01 * length(n2Claims)))
  ClaimIDs  = ClaimIDs[!ClaimIDs %in% n2Claims]
  Results =
    list(DtClaim  = DtClaim,
         ClaimIDs = ClaimIDs,
         n2Claims = n2Claims)
  Results
}

FindIntercept <- function(X, B, TargetPrev, yKnown, Range = c(-10, 10)) {
  RandGen <- function(N, p) rbinom(n = N, size = 1, prob = p)
  f <- function(B0, B, TargetPrev, X, Seed = 1, yKnown) {
    N = nrow(X)
    B = c(B0, B)
    p = binomial()$linkinv(as.vector(cbind(1, X) %*% B))
    set.seed(Seed)
    Y = do.call("RandGen", list(N = N, p = p))
    Prev = (sum(Y) + sum(yKnown)) / (N + length(yKnown))
    abs(Prev - TargetPrev)
  }
  Approx = F
  i = 1
  RangeB0 = T
  while(RangeB0) {
    B0 = optimise(f, Range, TargetPrev = TargetPrev, B = B, X = X, tol = .Machine$double.eps^0.5, Seed = i, yKnown = yKnown)$minimum
    B0 = round(B0, digits = 4)
    i  = i + 1
    RangeB0 <- B0 %in% c(-10, 10)
    if(i > 1e2)
      break
  }
  return(list(B0 = B0, Seed = i, Converged = i <= 1e2))
}

FactorToNumeric <- function(x) {
  if(!is.factor(x))
    stop("Provide variable of type factor.")
  as.numeric(as.character(x))
}

Prev <- function(x) {
  if(is.factor(x))
    x = FactorToNumeric(x)
  sum(x) / length(x)
}

unregister_dopar <- function() {
  ns              = getNamespace("foreach")
  .foreachGlobals = get(".foreachGlobals", envir = ns)
  rm(list = ls(name = .foreachGlobals), pos = .foreachGlobals)
}



AMH.sim <- function(u, v, alpha = -0.15) {
  # From package copula: https://github.com/cran/copula/blob/master/R/amhCopula.R
  w <- runif(length(u))
  b <- 1 - u
  A <- w * (alpha * b)^2 - alpha
  B <- (alpha + 1) - 2 * alpha * b * w
  C <- w - 1
  v <- (- B + sqrt(B^2 - 4 * A * C)) / 2 / A
  v <- 1 - v
  v
}

gender.cdf <- function(u, PropMale   = Gender$propMale, PropFemale = Gender$propFemale) {
  Gender <- NULL
  if (u <= PropFemale) {
    "female"
  }
  else if (u <= 1 - PropMale) {
    "non-binary"
  }
  else{
    "male"
  }
}


PlotGraph <- function(g, FrCl, FraudScores = NULL, PartyScores = NULL,
                      PosLegend = "topleft", vertex.label.dist = 2, vertex.label.degree = 90,
                      vertex.label.font = 2, vertex.label.cex = 1, ...) {
  igraph::V(g)$color = sapply(igraph::V(g)$name, function(x) {
    if (x %in% FrCl) {
      "red"
    } else if (grepl("C", x)) {
      "lightgreen"
    } else {
      "grey"
    }
  })
  igraph::V(g)$shape = sapply(igraph::V(g)$name, function(x) {
    if (grepl("C", x)) {
      "circle"
    } else {
      "square"
    }
  })
  if(!all(sapply(list(FraudScores, PartyScores), is.null))) {
    igraph::V(g)$score = sapply(igraph::V(g)$name, function(x) {
      if (grepl("C", x)) {
        round(FraudScores[which(names(FraudScores) == x)], 2)
      } else {
        round(PartyScores[which(names(PartyScores) == x)], 2)
      }
    })
    igraph::V(g)$scoreCol = sapply(igraph::V(g)$name, function(x) {
      if (grepl("C", x)) {
        "black"
      } else {
        "darkgrey"
      }
    })
    set.seed(1)
    plot.igraph(g, vertex.label = igraph::V(g)$score, vertex.label.dist = vertex.label.dist, vertex.label.degree = vertex.label.degree,
                vertex.label.font = vertex.label.font, vertex.label.cex = vertex.label.cex, vertex.label.color = igraph::V(g)$scoreCol,
                vertex.color = "white",
                vertex.frame.color = "white", edge.color="white", edge.label = NA, ...)
    plot.igraph(g, edge.label = NA, edge.color = 'black', add = T,
                vertex.label = igraph::V(g)$name, vertex.color = igraph::V(g)$color,
                vertex.label.color = 'black', main = "Social network",
                vertex.shape = igraph::V(g)$shape, ...)
  } else {
    plot.igraph(g, edge.label = NA, edge.color = 'black',
                vertex.label = igraph::V(g)$name, vertex.color = igraph::V(g)$color,
                vertex.label.color = 'black', main = "Social network",
                vertex.shape = igraph::V(g)$shape, ...)
  }
  legend(PosLegend, c("Unknown/ non-fraudulent claim", "Fraudulent claim", "Party (PH, broker, ...)"),
         pch = c(rep(21, 2), 22), col = "black", bty = "n", pt.bg = c("lightgreen", "red", "grey"),
         pt.cex = 2, y.intersp = 1.5)
}

AdjMatNeighbors <- function(aMat, Order = 2) {
  if(Order < 0)
    stop("Not possible")
  else if(Order == 0)
    return(Matrix(0, nrow(aMat), ncol(aMat)))
  else if(Order == 1)
    return(aMat)

  if(!Matrix::isSymmetric(aMat))
    stop("Only symmetric matrices allowed.")

  AllMat = list(aMat)

  for(i in seq_len(Order)[-1]) {
    N = AllMat[[i - 1]] %*% aMat
    diag(N)  = 0
    if(is(N, "sparseMatrix"))
      N@x[N@x > 1] = 1
    else
      N[N > 1] = 1
    AllMat[[i]] = N
  }
  for(i in seq_along(AllMat)[-length(AllMat)])
    N = N - AllMat[[i]]
  if(is(N, "sparseMatrix"))
    N@x[N@x < 0] = 0
  else
    N[N < 0] = 0
  return(N)
}


#' Function to examine whether the network exhibits signs of homophily
#'
#' This function take an object of class sfnData as input and computes the dyadicity and heterophilicity.
#'
#' @param SimObj an object of class \code{sfnData}, resulting from \code{\link{sfnGenerator}}
#'
#' @seealso [sfnGenerator()]
#' @return A vector with the dyadicity and heterophilicity.
#' @export
#'
#' @examples
#' \donttest{
#' data("SimObj")
#' Homophily(SimObj)
#' }
Homophily <- function(SimObj) {
  Fraud <- ClaimID <- Party <- startLabel <- startNode <- endLabel <- endNode <- Link <- NULL
  if(!is(SimObj, "sfnData"))
    stop("Only objects of type sfnData allowed.")
  Dt       = as.data.table(SimObj$Dt)
  AdjMat   = as.data.table(SimObj$AdjMatOrig)
  AdjMat[, Fraud := as.numeric(Dt$Fraud[match(ClaimID, Dt$ClaimID)])]
  setorder(Dt, ClaimID)
  setorder(AdjMat, ClaimID)

  ClaimLabel  = structure(as.numeric(Dt$Fraud), names = as.numeric(as.factor(Dt$ClaimID)))
  ClaimLabel  = ClaimLabel[order(as.numeric(names(ClaimLabel)))]
  FraudClaims = which(ClaimLabel == 1)
  if(length(FraudClaims) != sum(Dt$Fraud))
    stop("Something went wrong when making the vector with fraudulent claims.")

  AdjMat  = AdjMat[, `:=` (
    startNode = as.numeric(as.factor(Party)),
    endNode   = as.numeric(as.factor(ClaimID))
  )]
  AdjMatOrig  = copy(AdjMat)
  W       = sparseMatrix(i = AdjMat$endNode,
                         j = AdjMat$startNode,
                         x = 1)

  AdjMat = rbind(cbind(Matrix(0, nrow(W), nrow(W)), W),
                 cbind(t(W), Matrix(0, ncol(W), ncol(W))))
  N1     = W
  N2     = AdjMatNeighbors(AdjMat, 2)

  N2[seq_len(nrow(W)), seq_len(nrow(W))]

  CC       = N2[seq_len(nrow(W)), seq_len(nrow(W))]
  CCsum    = Matrix::summary(CC)
  CCsumAdj = tapply(CCsum$j, CCsum$i, unique)
  CCnetw   = data.table(startNode = rep(names(CCsumAdj), times = sapply(CCsumAdj, length)), endNode = unlist(CCsumAdj))
  CCnetw[, startLabel := fifelse(startNode %in% FraudClaims, "Fraud", "Legitimate")]
  CCnetw[, endLabel := fifelse(endNode %in% FraudClaims, "Fraud", "Legitimate")]
  CCnetw[, Link := paste0(startLabel, "-", endLabel)]
  setnames(CCnetw, c("startNode", "endNode"), c("from", "to"))

  n1  = length(FraudClaims)
  N   = NrUnique(CCnetw$from)
  n0  = N - n1
  m11 = sum(CCnetw$Link == "Fraud-Fraud")
  m00 = sum(CCnetw$Link == "Legitimate-Legitimate")
  m10 = sum(CCnetw$Link %in% c("Legitimate-Fraud", "Fraud-Legitimate"))
  M   = m11 + m00 + m10
  if(M != nrow(CCnetw))
    stop("Incorrect calculation of links.")

  p      = 2 * M / (N * (N - 1))
  m11bar = (n1 * (n1 - 1) * p) / 2
  m10bar = n1 * (N - n1) * p

  (Dyadicity = m11 / m11bar)
  (HeteroPh  = m10 / m10bar)
  return(c("Dyadicity" = Dyadicity, "Heterophily" = HeteroPh))
}


NrUnique <- function(x, na.rm = T) {
  if (na.rm)
    length(unique(na.omit(x)))
  else
    length(unique(x))
}


.checkTheta <- function(x) {
  defaultV = list(AgeVsGender = -0.15, AgeVsExposure = 0.15, AgeVsNrContracts = 0.95, AgeCarVsValueCar = -25)
  for(add in setdiff(names(defaultV), names(x)))
    x[[add]] = defaultV[[add]]
  if(!all(names(x) %in% names(defaultV)))
    stop("There are incorrect slot names for the theta list.")
  if(!between(x$AgeVsGender, -1, 1) | !between(x$AgeVsExposure, -1, 1) | !between(x$AgeVsNrContracts, -1, 1))
    stop("For the AMH copula, the theta value must lie in between -1 (negative dependence) and 1 (positive dependence).")
  return(x)
}


.checkExpert <- function(x) {
  defaultE = list(Sensitivity = 0.99, Specificity = 0.99)
  for(add in setdiff(names(defaultE), names(x)))
    x[[add]] = defaultE[[add]]
  if(!all(names(x) %in% names(defaultE)))
    stop("There are incorrect slot names for the ExpertJudgement list.")
  if(!between(x$Sensitivity, 0, 1) | !between(x$Specificity, 0, 1))
    stop("The sensitivity and specificity of the expert judgement must be a value in between 0 and 1.")
  return(x)
}


.checkCoefficients <- function(x) {
  defaultV = list(ClaimFrequency =
                    c(-2.17531,
                      c(log(0.85), log(0.75), log(0.7), log(0.6), log(0.55), log(0.6), log(0.7)), # Age policy holder
                      c(log(1.5), 0),                  # Gender = male & X
                      c(log(0.9), log(0.8), log(0.6)), # Age car
                      -0.12, -0.11,                    # Type of coverage
                      log(1.19),                        # Type of fuel = Diesel
                      c(0.12, 0.18, 0.34, 0.48, 0.54, 0.78) # Bonus-Malus scale
                    ),
                  ClaimSeverity =
                    c(
                      6.06,
                      c(log(0.85), log(0.75), log(0.85), log(0.85), log(1.15), log(1.25), log(1.5)), # Age policy holder
                      c(-0.16, 0.11),                      # Coverage
                      c(0.1, 0.15, 0.15, 0.15, 0.20, 0.30) # Bonus-Malus scale
                    ),
                  Fraud =
                    c(
                      -2.5,       # Intercept
                      0.2,        # Claim amount
                      -0.35,      # Claim age
                      2,          # Size first order nbh
                      -2,         # Size second order nbh
                      -1.5,       # Number of contracts
                      -2,         # Age policyholder
                      3           # n2.ratioFraud
                    )
                  )
  for(add in setdiff(names(defaultV), names(x)))
    x[[add]] = defaultV[[add]]
  if(!all(names(x) %in% names(defaultV)))
    stop("There are incorrect slot names for the Coefficients list.")
  return(x)
}


.checkFormulas <- function(x) {
  defaultV = list(ClaimFrequency = formula(~ AgePHBin + GenderPH + AgeCarBin + Coverage + Fuel + BonusMalusBin),
                  ClaimSeverity  = formula(~ AgePHBin + Coverage + BonusMalusBin),
                  Fraud = formula(~ ClaimAmount + ClaimAge + n1Size + n2Size + NrContractsPH + AgePH + n2.ratioFraud))
  for(add in setdiff(names(defaultV), names(x)))
    x[[add]] = defaultV[[add]]
  if(!all(names(x) %in% names(defaultV)))
    stop("There are incorrect slot names for the Formulas list.")
  return(x)
}

.fakeDt <- function() {
  DtE =
    structure(
      list(
        AgeCar = numeric(0),
        AgeCarBin = factor(levels = c("[0,5]", "(5,10]", "(10,16.1]", "(16.1,20]")),
        AgeCarScaled = numeric(0),
        AgePH = numeric(0),
        AgePHBin = factor(levels = c("[18,26]", "(26,30]", "(30,36]", "(36,50]", "(50,60]", "(60,65]", "(65,70]", "(70,79.8]")),
        AgePHOrig = numeric(0),
        AgePHScaled = numeric(0),
        BonusMalus = numeric(0),
        BonusMalusBin = factor(levels = c("[0,1)", "[1,2)", "[2,3)", "[3,7)", "[7,9)", "[9,11)", "[11,22]")),
        Broker = character(0),
        ClaimAge = numeric(0),
        ClaimAgeOrig = numeric(0),
        ClaimAmount = numeric(0),
        ClaimAmountOrig = numeric(0),
        ClaimDate = numeric(0),
        ClaimID = character(0),
        ContractID = character(0),
        Coverage = factor(levels = c("TPL", "PO", "FO")),
        Criminal = integer(0),
        Expert = character(0),
        ExpertJudgement = logical(0),
        ExpPH = numeric(0),
        ExpPHContracts = numeric(0),
        Fraud = logical(0),
        fraudScore = numeric(0),
        Fuel = factor(levels = c("Gasoline/LPG/Other", "Diesel")),
        Garage = character(0),
        GenderPH = factor(levels = c("female", "male", "non-binary")),
        IDPH = integer(0),
        Investigated = numeric(0),
        Lambda = numeric(0),
        n1.max = numeric(0),
        n1.med = numeric(0),
        n1.q1 = numeric(0),
        n1.size = numeric(0),
        n1Size = numeric(0),
        n2.binFraud = numeric(0),
        n2.max = numeric(0),
        n2.med = numeric(0),
        n2.q1 = numeric(0),
        n2.ratioFraud = numeric(0),
        n2.ratioFraudOrig = numeric(0),
        n2.ratioNonFraud = numeric(0),
        n2.ratioNonFraudOrig = numeric(0),
        n2.size = numeric(0),
        n2Size = numeric(0),
        NClaims = integer(0),
        nodeID = integer(0),
        nPersons = integer(0),
        NrContractsPH = numeric(0),
        OrigValueCar = numeric(0),
        Police = character(0),
        Policyholder = character(0),
        RateNrContracts = numeric(0),
        Rule1 = logical(0),
        Rule2 = logical(0),
        Rule3 = logical(0),
        TimeSinceClaim = numeric(0),
        ValueCar = numeric(0),
        ValueCarScaled = numeric(0)
      ),
      class = c("data.table", "data.frame")
    )
  DtE = DtE[rep(1, 5)]
  for(j in seq_along(DtE))
    if(is.numeric(DtE[[j]])) DtE[[j]] = runif(nrow(DtE))
  return(DtE)
}
