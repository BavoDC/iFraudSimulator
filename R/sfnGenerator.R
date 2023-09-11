#' Synthetic Insurance Fraud Network Data Generator
#'
#' Function to generate synthetic insurance fraud network data. This function corresponds to the simulation engine as described in
#' Campo, B.D.C. and Antonio, K. (2023). arXiv:2308.11659, available at https://arxiv.org/abs/2308.11659
#'
#' @param TargetPrev target class imbalance (i.e. the ratio of the number of fraudulent claims to the total number of claims).
#' @param NrPH total number of policyholders
#' @param NrExperts total number of experts in the network. Default is 1\% of NrPH.
#' @param NrBrokers total number of brokers in the network. Default is 1\% of NrPH.
#' @param NrGarages total number of garages in network. Default is 3\% of NrPH.
#' @param NrPersons total number of involved persons in network. Default is 150\% of NrPH.
#' @param ExcludeParties character vector to indicate which type of party (or parties) have to be excluded. By default, the expert is excluded. When no parties have to be excluded,
#' set \code{ExcludeParties = NULL}.
#' @param Age a list with the named parameters for generating the age of the policyholders. By default,
#' \code{Age = list(AvgAge = 40, SDAge = 15, RangeAge = c(18, 80))}. Hence, the average age is 40, the standard deviation 15, the minimum 18
#' and the maximum is 80.
#' @param Exposure a list with the named parameters for generating the exposure of the policyholders. By default,
#' \code{Exposure = list(AvgExp = 5, SDExp = 1.5, RangeExp = c(0, 20))}. Thus, the average exposure is 5, the standard deviation 1.5, the minimum 0
#' and the maximum 20.
#' @param Gender a list with the named parameters for generating the gender of the policyholders. By default,
#' \code{Gender = list(propMale   = 0.71, propFemale = 0.28, propNonBinary = 0.01)}. Hence, the proportion of males is 71\%, the proportion
#' of females is 28\% and the proportion of non-binaries is 1\%.
#' @param probBroker probability that a broker is involved in a policy contract. Default probability is 50\%.
#' @param stdize a function to standardize or normalize the covariates of the data-generating fraud model. This function should take vector
#' as input and return a vector. By default, the vector is normalized to the range [-1, 1] using the min-max feature scaling (see \code{\link{normalize}}).
#' @param Formulas named list with formulas for the \code{ClaimFrequency}, \code{ClaimSeverity} and \code{Fraud} model.
#' @param Coefficients named list with coefficients for the \code{ClaimFrequency}, \code{ClaimSeverity} and \code{Fraud} model.
#' @param zeta the parameter controlling the dependency between the claim frequency and claim severity. By default, we set \code{zeta = 0} so that the
#' claim frequency and claim severity are independent.
#' @param theta named list with the values for the dependence parameter in the copulas. The slot \code{AgeVsGender} specifies the dependency between the age of the policyholder
#' and the gender of the policyholder, the slot \code{AgeVsExposure} the dependency between the age of the policyholder and the exposure, the slot \code{AgeVsContracts}
#' the dependency between the age of the policyholder and the number of contracts and the slot \code{AgeCarVsValueCar} controls the dependency between the
#' age of the car and the value of the car. If the list contains just one entry, the remaining entries will be set to their default values.
#' @param BusinessRules logical expressions to flag suspicious claims
#' @param Parallel logical, indicates whether parallel computing has to be used.
#' @param NrCores the number of cores that are utilized for parallel computing .
#' @param printProgress logical, indicates whether the progress has to be printed.
#' @param tmpFiles logical, indicating whether temporary files have to be saved while generating the synthetic data set. Useful when generating
#' large data sets.
#' @param Seed the seed that is set at the beginning of the synthetic data generation. For reproducibility.
#'
#' @importFrom Hmisc capitalize
#'
#' @return An object of type \code{sfnData} with the following slots:
#' @return \item{call}{the matched call.}
#' @return \item{Dt}{the synthetic data set.}
#' @return \item{SummaryPlots}{a ggplot object with the summary plots of the synthetic data set.}
#' @return \item{TargetPrev}{the target class imbalance.}
#' @return \item{TruePrevalence}{the class imbalance in the synthetic data set.}
#' @return \item{AdjMatOrig}{a \code{data.table} which contains all information to construct the adjacency matrix.}
#' @return \item{AdjMat}{the adjacency matrix of the bipartite graph.}
#' @return \item{BiRank}{the object resulting from running the BiRank algorithm using the function \code{\link{BiRankFr}}.}
#' @export
#'
#' @examples
#' \dontrun{
#' SimObj = sfnGenerator(TargetPrev = 0.05, NrPH = 5000, Seed = 1,
#'  printProgress = FALSE, Parallel = FALSE)
#'  }
sfnGenerator <- function(TargetPrev = 0.01,
                   NrPH       = 1e4,
                   NrExperts  = floor(0.01 * NrPH),
                   NrBrokers  = floor(0.01 * NrPH),
                   NrGarages  = floor(0.03 * NrPH),
                   NrPersons  = NrPH * 1.5,
                   ExcludeParties = "Expert",
                   Age        = list(AvgAge = 40, SDAge = 15, RangeAge = c(18, 80)),
                   Exposure   = list(AvgExp = 5, SDExp = 1.5, RangeExp = c(0, 20)),
                   Gender     = list(propMale   = 0.71, propFemale = 0.28, propNonBinary = 0.01),
                   probBroker = 0.5,
                   stdize = normalize,
                   Formulas =
                     list(ClaimFrequency = formula(~ AgePHBin + GenderPH + AgeCarBin + Coverage + Fuel + BonusMalusBin),
                          ClaimSeverity  = formula(~ AgePHBin + Coverage + BonusMalusBin),
                          Fraud = formula(~ ClaimAmount + ClaimAge + n1Size + n2Size + NrContractsPH + AgePH + n2.ratioFraud)),
                   Coefficients =
                     list(ClaimFrequency =
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
                     ),
                   zeta = 0,
                   theta = list(AgeVsGender = -0.15, AgeVsExposure = 0.15, AgeVsNrContracts = 0.95),
                   BusinessRules =
                     list(
                       Expressions = list(
                         "cumsum(ClaimAmountOrig) > 2 * ValueCar",
                         "ClaimAmountOrig * 0.75 > ValueCar",
                         "c(F, diff(ClaimDate) <= 1)"),
                       byVariable = list(
                         "ContractID",
                         NULL,
                         "ContractID"
                       )),
                   Parallel = TRUE,
                   NrCores = detectCores() - 2,
                   printProgress = TRUE,
                   tmpFiles = FALSE,
                   Seed = 7032018) {

  #### 0.1 Packages and custom functions ####
  AvgAge <- SDAge <- RangeExp <- AvgExp <- SDExp <- IDPH <- ContractID <- AgeCar <- ExpPHContracts <- ClaimDate <- ClaimAge <- NULL
  Policyholder <- Expert <- ClaimAmount <- ClaimAmountOrig <- nPersons <- ClaimID <- Party <- variable <- endNode <- Criminal <- TimeSinceClaim <- NULL
  Fraud <- AgePHOrig <- ClaimAgeOrig <- n2Claims <- ClaimIDs <- ExpertJudgement <- Investigated <- nodeID <- fraudScore <- AgeClaim <- nodeID <- NULL
  call    = match.call()
  Formals = formals(sfnGenerator)
  if(is.null(ExcludeParties))
    call$ExcludeParties = "None"
  for(i in setdiff(names(Formals), names(call)))
    call[[i]] = Formals[[i]]

  if(!is.null(ExcludeParties)) {
    if(!all(ExcludeParties %in% c("Expert", "Broker", "Garage")))
      stop("Parties to be excluded must be one of 'expert', 'broker' or 'garage'.")
  }

  if(floor(TargetPrev * floor(1e-2 * NrPH)) < 1)
    stop("The number of policyholders is too low to reach the desired target prevalence. Please increase the NrPH.")

  list2env(Age, envir = environment())
  list2env(Exposure, envir = environment())

  if(!is.function(stdize))
    stop("Argument stdize has to be a function.")

  if(sum(unlist(Gender)) != 1)
    stop("The proportions specified in the argument Gender do not sum to 1.")

  if(printProgress) {
    cat("\f")
    for(i in 1:3) {
      cat("\rStarting the simulation", rep(".", i))
      Sys.sleep(0.25)
      flush.console()
    }
  }

  if(Parallel)
    if(NrCores >= detectCores()) {
      msg =
        paste("Number of cores is equal to or larger than the number of CPU cores.",
              "Will be set to one less than the number of cores to ensure sufficient memory.")
      warning(msg, immediate. = TRUE)
      NrCores = detectCores() - 1
    }

  theta = .checkTheta(theta)


  #### 0.2 Dashboard for simulation settings ####

  # definieer beschikbare netwerkpartijen (garages, tussenpersonen, experten en polishouders)
  garages.kandidaten <- c("Atelier Reparation","Auto Otto", "Garage Per Total", "Delessai", "Edcar", "Eliforp",
                          "Garage Ladeuze", "Garage Naamsestraat", "Georgette Voitures", "Happycar", "Garage Leeuwerik",
                          "Lady Limousine", "Bricoler & Tenter", "Occasie Olivier", "Quai Branly", "Karel Kilometerteller",
                          "Slow Fit", "Meneertje Trabant", "Wagens Rouleau", "Wagens Rousseau", "Vehicles Luxe Bene",
                          "Vehicles Bene Luxe", paste("Garage", 1:(NrGarages - 22)))

  tussenpersonen.kandidaten <- capitalize(c("handige harry","louche leontien", "madame jackpot", "monsieur malfait",
                                                   "sjacheraar jacque", "vlugge japie", "yvette integer", "vincent viable",
                                                   paste("Broker", 1:(NrBrokers - 8))))

  experten.kandidaten <- capitalize(c("simba", "mufasa", "christina", "britney", "gazoo", "barney", "bobientje",
                                             "chicolama", "pino", paste("Expert", 1:(NrExperts - 9))))

  Coverage   = data.frame(Type = c("TPL", "PO", "FO"), Prop = c(0.6, 0.25, 0.15))
  Fuel       = data.frame(Type = c("Gasoline/LPG/Other", "Diesel"), Prop = c(0.7, 0.3))
  BonusMalus = data.frame(Scale = 0:20, Prop = dgamma(0:20, 1, 1 / 3))


  #### 1. Simulation dataset ####
  if(printProgress)
    cat("\014Simulating policyholder characteristics")
  set.seed(Seed)

  #### 1.1 Simulation age PH ####
  RangeAge = sort(RangeAge)
  SeqAge   = min(RangeAge):max(RangeAge)
  AgePH    = sort(rnorm(NrPH, AvgAge, SDAge))
  OutsideRange = any(AgePH < RangeAge[1] | AgePH > RangeAge[2]) # Check if observations with age outside of range
  if(OutsideRange) {
    nRange    = sum(AgePH >= RangeAge[1] & AgePH <= RangeAge[2])
    probAge   = prop.table(table(ceiling(AgePH[AgePH >= RangeAge[1] & AgePH <= RangeAge[2]])))
    probAge   = probAge[!names(probAge) %in% RangeAge]
    probAge   = probAge / sum(probAge)
    AdjAge    = rmultinom(1, NrPH - nRange, probAge) # Redistribute age outside range
    AdjAge    = as.numeric(rep(rownames(AdjAge), AdjAge)) + runif(NrPH - nRange)
    AgePH = c(AgePH[AgePH >= RangeAge[1] & AgePH <= RangeAge[2]], AdjAge)
  }
  if(any(AgePH < RangeAge[1] | AgePH > RangeAge[2]))
    stop("Error in simulation age, out of range")

  Df      = data.frame(AgePH = AgePH)
  HistAge = ggplot(Df, aes(AgePH)) +
    geom_histogram(binwidth = 1,
                   color = "white",
                   fill = "steelblue")

  #### 1.2 Simulation gender ####
  # Proportion non-binary: https://www.nature.com/articles/s41598-021-81411-4

  AMHCopula = AMH.sim(runif(NrPH), pnorm(AgePH, AvgAge, SDAge), alpha = theta$AgeVsGender)
  GenderPH  = sapply(AMHCopula, gender.cdf, PropMale = Gender$propMale, PropFemale = Gender$propFemale)
  Df = cbind.data.frame(Df, data.frame(GenderPH = GenderPH))
  HistGender    = ggplot(Df, aes(GenderPH)) +
    geom_bar(color = "white", fill = "steelblue")
  HistAgeGender = ggplot(Df, aes(x = AgePH, fill = GenderPH)) + geom_density(alpha = .3)


  #### 1.3 Simulation duration policy using normal distribution ####
  MinExp    = min(RangeExp)
  MaxExp    = max(RangeExp)
  MaxExp    = pmin(AgePH - RangeAge[1], MaxExp)
  AMHCopula = AMH.sim(runif(NrPH), pnorm(AgePH, AvgAge, SDAge),alpha = theta$AgeVsExposure)
  ExpPH     = sapply(AMHCopula, function(x) qnorm(x, AvgExp, SDExp))
  if(any(ExpPH < MinExp | ExpPH > MaxExp)) {
    OutsideRange = which(ExpPH < MinExp | ExpPH > MaxExp)
    for(i in OutsideRange) {
      Conv = F
      while(!Conv) {
        ExpPH[i] = runif(1, MinExp, MaxExp[i])
        Conv = ExpPH[i] >= MinExp & ExpPH[i] <= MaxExp[i]
      }
    }
  }

  Df = cbind.data.frame(Df, data.frame(ExpPH = ExpPH))
  HistExp = ggplot(Df, aes(ExpPH)) +
    geom_histogram(binwidth = 1,
                   color = "white",
                   fill = "steelblue")
  ScatterExpAge = ggplot(Df, aes(ExpPH, AgePH)) + geom_point() + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))



  #### 1.4 Simulation number of contracts & cars using Poisson distribution with cubic budget rate via age ####
  ContractBudgetPar  = 0.25
  RateNrContracts = ContractBudgetPar * (1.05 + -0.0000025 * AgePH + 0.0025 * AgePH^2 - 0.0000265 * AgePH^3)
  RangeContracts  = c(1, 5)

  AMHCopula = AMH.sim(runif(NrPH), pnorm(AgePH, AvgAge, SDAge), alpha = theta$AgeVsNrContracts)
  NrContractsPH = qpois(AMHCopula, RateNrContracts)
  NrContractsPH = pmin(pmax(NrContractsPH, RangeContracts[1]), RangeContracts[2])
  Df = cbind.data.frame(Df, data.frame(NrContractsPH = as.factor(NrContractsPH), RateNrContracts = RateNrContracts))

  HistContracts    = ggplot(Df, aes(NrContractsPH)) + geom_bar(color = "white", fill = "steelblue")
  HistAgeContracts = ggplot(Df, aes(x = AgePH, fill = NrContractsPH)) + geom_density(alpha = 0.3)

  SummaryPlots =
    arrangeGrob(
      HistAge, HistGender,
      HistAgeGender, HistExp,
      ScatterExpAge, HistContracts,
      HistAgeContracts,
      ncol = 2, nrow = 4, layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 5), c(6, 7))
      )


  #### 1.5 Simulate claim frequency and severity ####
  if(printProgress)
    cat("\014Simulating claim frequency and severity")
  # Horsepower: https://www.hotcars.com/highest-and-lowest-horsepower-domestic-cars/
  #             https://www.cnet.com/roadshow/pictures/least-powerful-lowest-horsepower-cars-you-can-buy-pictures/
  # Car depreciation: https://goodcalculators.com/car-depreciation-calculator/
  # Cheapest cars: https://gocar.be/nl/autonieuws/autosalon/de-goedkoopste-modellen-autosalon-brussel-2020
  # Avg = 0.15. BMW 0.075

  Df$IDPH       = seq_len(NrPH)
  Df$NrContractsPH = as.numeric(as.character(Df$NrContractsPH))

  gc()
  Sys.sleep(1)
  if(Parallel & NrPH > 1e4) {
    Cl = parallel::makeCluster(NrCores)
    registerDoParallel(Cl)
  }
  if(printProgress)
    cat("\014Simulating claim frequency and severity:\n simulating contract-specific characteristics")
  Df = ddply(Df, .(IDPH), function(x) {
    Gender = x$GenderPH
    nContr = x$NrContractsPH
    if(x$NrContractsPH > 1) {
      ExpPHContracts = c(x$ExpPH, x$ExpPH - runif(x$NrContractsPH - 1, 0, x$ExpPH / 2))
      x = x[rep(1, nContr), ]
      x$ExpPHContracts = ExpPHContracts
      x$ContractID = paste0(x$IDPH, "_", seq_len(nrow(x)))
      uFrank = rCopula(nContr, frankCopula(theta$AgeCarVsValueCar, 2))
      uFrank[, 1] = uFrank[, 1] * (0.75 - 0.4) + 0.4 # Bound to range 0.4 to 0.75
      x$AgeCar    = pmax(qnorm(uFrank[, 2], 7.5, sqrt(5)), x$ExpPHContracts)
      ValueCar    = qexp(uFrank[, 1], x$RateNrContracts / nContr) * if(Gender == "male") 25e3 else if(Gender == "female") 20e3 else 22.5e3
      x$OrigValueCar = ValueCar
      x$ValueCar     = ValueCar * (1 - ifelse(ValueCar < 30e3, 0.15, 0.075))^x$AgeCar
    } else {
      x$ContractID = paste0(x$IDPH, "_", 1)
      x$ExpPHContracts = x$ExpPH
      uFrank = rCopula(nContr, frankCopula(-25, 2))
      uFrank[, 1]  = uFrank[, 1] * (0.75 - 0.4) + 0.4 # Bound to range 0.4 to 0.75
      x$AgeCar     = pmax(qnorm(uFrank[, 2], 7.5, sqrt(5)), x$ExpPHContracts)
      ValueCar     = qexp(uFrank[, 1], x$RateNrContracts / nContr) * if(Gender == "male") 25e3 else if(Gender == "female") 20e3 else 22.5e3
      x$OrigValueCar = ValueCar
      x$ValueCar   = ValueCar * (1 - ifelse(ValueCar < 30e3, 0.15, 0.075))^x$AgeCar
    }
    return(x)
  }, .parallel = Parallel & NrPH > 1e4)


  for(i in c("ValueCar", "AgeCar", "AgePH"))
    Df[[paste0(i, "Scaled")]] = stdize(Df[[i]])

  Df = ddply(Df, .(IDPH, ContractID), function(x) {
    X = model.matrix(~ ValueCarScaled + AgeCarScaled + AgePHScaled - 1, data = x)
    B = list(
      c(log(0.5), log(1.25), log(0.25)),
      c(log(1.25), log(0.75), log(1.05)),
      c(log(1.5), log(0.75), log(1.25))
    )
    P = do.call("cbind", lapply(B, function(b) exp(X %*% b)))
    Y = t(apply(P, 1, function(x) which(rmultinom(x, n = 1, size = 1) == 1)))
    x$Coverage   = Coverage$Type[Y]
    x$Fuel       = replicate(nrow(x), sample(Fuel$Type, 1, prob = Fuel$Prop))
    x$BonusMalus = rep(floor(rgamma(1, 1, 1 / 3)), nrow(x))
    x$BonusMalus = ifelse(x$BonusMalus > 22, 22, x$BonusMalus)
    return(x)
  }, .parallel = Parallel & NrPH > 1e4)
  if(Parallel & NrPH > 1e4)
    stopCluster(Cl)

  Dt = as.data.table(Df)
  Dt[, BonusMalus := rep(max(BonusMalus), .N), by = "IDPH"]

  Dt[,  `:=` (
    AgePHBin      = cut(AgePH, c(18, 26, 30, 36, 50, 60, 65, 70, max(AgePH)), include.lowest = T),
    AgeCarBin     = cut(AgeCar, c(0, 5, 10, 20, max(AgeCar)), include.lowest = T),
    BonusMalusBin = cut(BonusMalus, c(0, 1, 2, 3, 7, 9, 11, max(BonusMalus)), include.lowest = T, right = F)
  )]
  Dt[, `:=` (
    GenderPH = factor(GenderPH, levels = c("female", "male", "non-binary")),
    Coverage = factor(Coverage, levels = c("TPL", "PO", "FO")),
    Fuel     = factor(Fuel, levels = c("Gasoline/LPG/Other", "Diesel"))
  )]

  if(printProgress)
    cat("\014Simulating claim frequency and severity:\n claim frequency")
  CoefPois = Coefficients$ClaimFrequency
  X = as(model.matrix(Formulas$ClaimFrequency, data = Dt), "sparseMatrix")

  ## Check input ##
  if(ncol(X) != length(CoefPois)) {
    msg =
      c("Error claim frequency model: number of coefficients != number of columns model matrix based\n",
        paste0("\nFormula:\n", strwrap(capture.output(print(Formulas$ClaimFrequency)))),
        "\n\nCoefficients:\n",
        paste(capture.output(print(head(X))), collapse = "\n"))
    warning(msg, immediate. = T)
  }

  L = exp(X %*% CoefPois + log(Dt$ExpPHContracts))
  N = rpois(nrow(Dt), as.vector(L))
  rm(L)

  Dt$NClaims = N
  DtClaim    = Dt[rep(1:nrow(Dt), N), ]

  DtClaim[, `:=` (
    Garage    = sample(garages.kandidaten, .N, T),
    Expert    = sample(experten.kandidaten, .N, T),
    Broker    = ifelse(runif(1) < probBroker, sample(tussenpersonen.kandidaten, .N, T), rep(NA, .N)),
    nPersons  = replicate(.N, sample(0:5, 1, prob = c(0.025, 0.6, 0.2, 0.1, 0.1, 0.025))),
    ClaimID   = paste0(ContractID, seq_len(.N)),
    Police    = replicate(.N, sample(c("Yes", "No"), 1, prob = c(0.25, 0.75))),
    ClaimAge  = pmin(sort(floor(rexp(.N, 0.25))), floor(unique(ExpPHContracts) * 12)),
    ClaimDate = sort(runif(.N, 0, max = unique(ExpPHContracts)))
  ), by = ContractID]

  DtClaim[, ClaimDate := fifelse(ClaimDate <= ClaimAge / 12, ClaimAge / 12, ClaimDate)]
  DtClaim[, TimeSinceClaim := ExpPH - ClaimDate]
  # Sanity check
  if(any(DtClaim$TimeSinceClaim < 0))
    DtClaim[TimeSinceClaim < 0, TimeSinceClaim := 0]


  DtClaim[, Policyholder := rep(randomNames::randomNames(1, if(unique(GenderPH) == "non-binary") NA else if(unique(GenderPH) == "male") 0 else 1), .N),
          by = IDPH]

  ## Claim severity ##
  # Frees: Predictive Modeling Applications in Actuarial Science
  if(printProgress)
    cat("\014Simulating claim frequency and severity:\n Claim severity")
  CoefGamma = Coefficients$ClaimSeverity
  X = as(model.matrix(Formulas$ClaimSeverity, data = DtClaim), "sparseMatrix")

  ## Check input ##
  if(ncol(X) != length(CoefGamma)) {
    warning("Error claim severity model: number of coefficients != number of columns model matrix based", immediate. = T)
    cat("\nFormula:\n", strwrap(capture.output(print(Formulas$ClaimSeverity))))
    cat("\nCoefficients:\n")
    print(head(X))
  }

  muTrue   =
    if(zeta == 0) {
      exp(X %*% CoefGamma)
    } else {
      exp(X %*% CoefGamma + zeta * DtClaim$NClaims)
    }
  nu       = 0.25
  alpha1  = nu
  p       = (2 + alpha1) / (1 + alpha1)
  p
  Lambda  = as.vector(nu / muTrue)
  DtClaim[, Lambda := Lambda]
  DtClaim[, `:=` (
    ClaimAmount = max(rgamma(1, shape = nu, rate = Lambda), runif(1, 50, 150))
  ), by = "ClaimID"]

  # Set expert to NA if claim amount < 250
  # https://www.kbcbrussels.be/retail/en/products/insurance/vehicle/car-claims/accident-with-no-dispute
  # -or-claim-on-comprehensive-insurance.html#parsys1_chaptertitle_d573

  DtClaim[, Expert := ifelse(ClaimAmount < 250, NA, Expert), by = "ClaimID"]
  DtClaim[, ClaimAmountOrig := ClaimAmount]

  #### 1.5.1 Adjacency matrix #####
  if(printProgress)
    cat("\014Preparing for the simulation of fraudulent claims")
  DtG = DtClaim[!duplicated(IDPH)]
  if(sum(DtG$GenderPH == "non-binary") < 1e2)
    DtClaim = DtClaim[GenderPH != "non-binary"] # Creates problems if number is too low

  personNames = randomNames(NrPersons) %>% .[!. %in% unique(DtClaim$Policyholder)]
  AdjMat1 = melt(DtClaim, id.vars = "ClaimID", measure.vars = c("Garage", "Expert", "Broker", "Policyholder"),
                 value.name = "Party")
  AdjMat2 = DtClaim[nPersons > 0, .(variable = rep("Person", nPersons), Party = sample(personNames, nPersons, F)), by = ClaimID]
  AdjMat  = rbind(AdjMat1, AdjMat2)
  setorder(AdjMat, ClaimID, Party)
  AdjMat  = AdjMat[!is.na(Party)]
  if(!is.null(ExcludeParties)) {
    if(printProgress)
      cat("\nExcluding", ExcludeParties)
    AdjMat = AdjMat[!variable %in% ExcludeParties]
    AdjMat[, variable := droplevels(variable)]
  }
  AdjMat  = AdjMat[, `:=` (
    startNode = as.numeric(as.factor(Party)),
    endNode   = as.numeric(as.factor(ClaimID))
  )]
  AdjMatOrig = copy(AdjMat)
  W       = sparseMatrix(i = AdjMat$endNode,
                         j = AdjMat$startNode,
                         x = 1)
  AdjMat = rbind(cbind(Matrix(0, nrow(W), nrow(W)), W),
                 cbind(t(W), Matrix(0, ncol(W), ncol(W))))
  N1     = W
  N2     = AdjMatNeighbors(AdjMat, 2)
  n1Size = rowSums(N1)
  n2Size = colSums(N2[1:nrow(W), 1:nrow(W)])
  setorder(AdjMatOrig, endNode)
  AdjMatOrig[, `:=` (
    n1Size = n1Size[endNode],
    n2Size = n2Size[endNode]
  )]

  ## Sanity Check ##
  # head(AdjMatOrig)
  # Node = 2212
  # AdjMatOrig[endNode == Node]
  # AdjMatOrig[endNode %in% Node][["startNode"]] %>% NrUnique()
  # AdjMatOrig[startNode %in% AdjMatOrig[endNode == Node][["startNode"]]][endNode != Node][["endNode"]] %>% NrUnique

  DtClaim[, `:=` (
    n1Size = AdjMatOrig$n1Size[match(ClaimID, AdjMatOrig$ClaimID)] ,
    n2Size = AdjMatOrig$n2Size[match(ClaimID, AdjMatOrig$ClaimID)]
  )]

  CoefCriminal = c(
    -5,     # Intercept
    c(log(3), 0) # Gender = male & non-binary
  )

  X = model.matrix(~ GenderPH, data = DtClaim)
  p = binomial()$linkinv(X %*% CoefCriminal)
  Y = rbinom(nrow(DtClaim), 1, p)
  DtClaim[, Criminal := Y]
  # https://www.gov.uk/government/statistics/women-and-the-criminal-justice-system-2019/women-and-the-criminal-justice-system-2019


  #### 1.6 Simulate fraudulent claims ####
  if(tmpFiles) {
    save(list = ls(all.names = TRUE), file = "tempResults.RData")
    Sys.sleep(60)
  }

  if(Parallel) {
    unregister_dopar()
    Cl = parallel::makeCluster(NrCores)
    registerDoParallel(Cl)
  }


  if(printProgress)
    cat("\014Simulating fraudulent claims")
  ## Build network ##
  DtClaim[, Fraud := NA]

  DtClaim[, AgePHOrig := copy(AgePH)]
  DtClaim[, ClaimAgeOrig := copy(ClaimAge)]
  ScaleVars = c("ClaimAmount", "ClaimAge", "n1Size", "n2Size", "AgePH")
  DtClaim[, (ScaleVars) := lapply(.SD, function(x) stdize(x)), .SDcols = ScaleVars]
  ## First batch ##
  prevFraud = 1
  iter      = 0

  TempForm = all.vars(Formulas$Fraud) %>% .[!grepl("n1\\.|n2\\.", .)]
  TempForm = formula(paste("~ ", paste(TempForm, collapse = "+")))
  X        = as(model.matrix(TempForm, data = DtClaim), "sparseMatrix")
  CoefLog  = Coefficients$Fraud[1:ncol(X)]

  while(!(prevFraud < TargetPrev & prevFraud > 0)) {
    SelectBatch = sample(unique(DtClaim$IDPH),
                         floor(if(NrPH < 1e4) 0.05 * NrUnique(DtClaim$IDPH) else 0.01 * NrUnique(DtClaim$IDPH)), F)
    DtSubset    = DtClaim[IDPH %in% SelectBatch]
    X   = as(model.matrix(TempForm, data = DtSubset), "sparseMatrix")
    Res = FindIntercept(X[, -1], CoefLog[-1], TargetPrev, yKnown = vector())
    CoefLog[1] = Res$B0
    p = binomial()$linkinv(as.vector(X %*% CoefLog))
    if(mean(p) > TargetPrev) {
      Res = FindIntercept(X[, -1], CoefLog[-1], TargetPrev, yKnown = vector(), Range = c(-10, 0))
      CoefLog[1] = Res$B0
      p = binomial()$linkinv(as.vector(X %*% CoefLog))
    }
    set.seed(Res$Seed)
    Y = rbinom(nrow(DtSubset), 1, p)
    (prevFraud = Prev(Y))
    DtSubset[, Fraud := Y]
    iter  = iter + 1
    if(iter > 5e3)
      stop("Stuck in loop")
  }

  ResNetwork = AddInfo(DtClaim, DtSubset, SelectBatch, AdjMatOrig)
  list2env(ResNetwork, envir = environment())
  iter       = 0
  nSample    = floor(0.1 * nrow(DtClaim))
  NrLoop     = 0
  CoefLogistic = Coefficients$Fraud
  DtSubset = AddInfoN2(DtSubset, AdjMatOrig, Parallel, std = stdize, DtCl = DtClaim)
  X = as(model.matrix(Formulas$Fraud, data = DtSubset), "sparseMatrix")

  ## Check formula object ##
  if(ncol(X) != length(CoefLogistic)) {
    warning("Error fraud model: number of coefficients != number of columns model matrix", immediate. = T)
    cat("\nFormula:\n", strwrap(capture.output(print(Formulas$Fraud))))
    cat("\nCoefficients:\n")
    print(head(X))
  }

  if(printProgress)
    cat("\014\nProgress simulation fraudulent claims:\n")
  while(anyNA(DtClaim$Fraud)) {
    iter  = iter + 1
    RemCl = sumNA(DtClaim$Fraud) / nrow(DtClaim)
    if(printProgress) {
      cat("\r", round(RemCl * 1e2), "% of claims unlabeled")
    }
    SelectBatch = c(n2Claims, if(RemCl <= 0.1) unique(ClaimIDs) else sample(ClaimIDs, nSample, F))
    DtSubset    = DtClaim[ClaimID %in% SelectBatch]
    DtSubset    = DtSubset[is.na(DtSubset$Fraud)]
    DtSubset    = AddInfoN2(DtSubset, AdjMatOrig, Parallel, std = stdize, DtCl = DtClaim)


    X = as(model.matrix(Formulas$Fraud, data = DtSubset), "sparseMatrix")
    yKnown = as.numeric(DtClaim[!is.na(Fraud), get("Fraud")])
    Res = FindIntercept(X[, -1], CoefLogistic[-1], TargetPrev, yKnown)
    if(!Res$Converged) {
      NrLoop = NrLoop + 1
      if(NrLoop < 5)
        next
    } else {
      NrLoop = 0
    }
    if(NrLoop == 0)
      CoefLogistic[1] = Res$B0
    p = binomial()$linkinv(as.vector(X %*% CoefLogistic))
    set.seed(Res$Seed)
    Y = rbinom(nrow(DtSubset), 1, p)
    (prevFraud = Prev(Y))
    DtSubset[, Fraud := Y]
    ResNetwork = AddInfo(DtClaim, DtSubset, SelectBatch, AdjMatOrig)
    list2env(ResNetwork, envir = environment())
  }

  DtClaim = AddInfoN2(DtClaim, AdjMatOrig, Parallel, std = stdize, DtCl = DtClaim)

  #### 1.7 Claims that are investigated ####
  # Business rule 1: if cumulative claim amount > 200% value car, investigation
  # Business rule 2: if individual claim amount > 75% value car, investigation
  # Business rule 3: if new claim within a year, investigation
  # DtClaim[, Rule1 := cumsum(ClaimAmountOrig) > 2 * ValueCar, by = "ContractID"]
  # DtClaim[, Rule2 := ClaimAmountOrig * 0.75 > ValueCar]
  # DtClaim[, Rule3 := c(F, diff(ClaimDate) <= 1), by = ContractID]
  # DtClaim[, Investigated := as.numeric(Rule1 | Rule2 | Rule3)]
  exprTmp = BusinessRules$Expressions
  exprTmp = lapply(seq_along(exprTmp), function(i) paste0("Rule", i, " := ", exprTmp[[i]]))

  for(i in seq_along(exprTmp))
    DtClaim[, eval(parse(text = exprTmp[[i]])), by = eval(BusinessRules$byVariable[[i]])]

  varTmp = paste0("Investigated := as.numeric(", paste0("Rule", seq_along(exprTmp), collapse = "|"), ")")
  DtClaim[, eval(parse(text = varTmp))]


  DtClaim[, ExpertJudgement := if(all(Investigated == 0)) rep(NA, .N) else if(all(Fraud == T)) rbinom(.N, 1, 0.99) else rbinom(.N, 1, 0.01),
          by = c("Investigated", "Fraud")]




  #### 1.8 Birank algorithm ####
  if(printProgress)
    cat("\014Running Birank algorithm")
  setorder(DtClaim, ClaimID)
  DtClaim[, nodeID := as.numeric(as.factor(ClaimID))]
  setorder(DtClaim, nodeID)
  c0 = DtClaim[, get("ExpertJudgement")]
  c0 = ifelse(is.na(c0), 0, c0)

  ResultsNetw = BiRankFr(AdjMatOrig, data.frame(FraudInd = c0))

  #### 1.8.1 Features based on fraud ####
  if(printProgress)
    cat("\014Computing fraud based features.")
  aMat = ResultsNetw$AdjacencyMatrix
  ResP = ResultsNetw$ResultsParties
  ResP = ResP[order(ResP$ID), "ScaledScore"]
  ResC = ResultsNetw$ResultsClaims
  ResC = ResC[order(ResC$ID), "ScaledScore"]

  N1   = aMat * ResP
  N2   = t(t(aMat) * ResC)

  ## Functions ##
  Q1 <- function(x) quantile(x, 0.25)
  f  = list(Q1, median, max, length)

  #### 1.8.1.1 First order neighborhood ####
  SumN1    = Matrix::summary(N1)
  N1scores = tapply(SumN1$x, SumN1$j, function(x) do.call("cbind", lapply(f, function(f) f(x))))
  N1scores = do.call("rbind", N1scores)
  colnames(N1scores) = c("n1.q1", "n1.med", "n1.max", "n1.size")
  rownames(N1scores) = seq_len(ncol(N1))


  #### 1.8.1.2 Second order neighborhood ####
  SumN2 = Matrix::summary(N2)
  N1P   = tapply(AdjMatOrig$endNode, AdjMatOrig$startNode, unique)
  N1Cl  = tapply(AdjMatOrig$startNode, AdjMatOrig$endNode, unique)
  FrCl  = DtClaim$nodeID[which(DtClaim$ExpertJudgement == 1)]
  NFrCl = DtClaim$nodeID[which(DtClaim$ExpertJudgement == 0)]
  if(length(c(FrCl, NFrCl)) != sum(DtClaim$Investigated == 1))
    stop("Something went wrong when creating the variable FraudInd!")

  NCl   = ncol(aMat)

  rm(aMat, N1, N2, SumN1)
  gc()

  if(Parallel) {
    stopCluster(Cl)
    Cl = makeCluster(NrCores, type = "SOCK")
    registerDoSNOW(Cl)
  }
  if(printProgress)
    pb       = txtProgressBar(max = NCl, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts     = list(progress = progress)
  N2scores = foreach(i = seq_len(NCl), .options.snow = if(printProgress) opts else NULL) %dopar%
    {
      x    = N1Cl[[i]]
      N2C  = unique(do.call("c", N1P[x]))
      N2C  = N2C[!N2C %in% i]
      if(length(N2C) != 0 ) {
        ClSc = ResC[N2C]
        NFr  = sum(N2C %in% FrCl)
        NnFr = sum(N2C %in% NFrCl)
        n    = length(ClSc)

        Res1 = do.call("c", lapply(f, function(f) f(ClSc)))
        Res2 = c(
          ratioFr = if (NFr == 0) 0 else NFr / n,
          ratioNonFr = if (NnFr == 0) 0 else NnFr / n,
          binFraud = as.numeric(NFr != 0)
        )
      } else {
        Res1 = rep(0, 4)
        Res2 = rep(0, 3)
      }
      return(c(Res1, Res2))
    }
  if(printProgress)
    close(pb)
  if(Parallel)
    stopCluster(Cl)

  N2scores = do.call("rbind", N2scores)
  colnames(N2scores) = c("n2.q1", "n2.med", "n2.max", "n2.size", "n2.ratioFraud",
                         "n2.ratioNonFraud", "n2.binFraud")
  rownames(N2scores) = seq_len(nrow(N2scores))
  N2scores = as.data.frame(N2scores)
  N2scores$n2.ratioFraudOrig    = N2scores$n2.ratioFraud
  N2scores$n2.ratioNonFraudOrig = N2scores$n2.ratioNonFraud
  N2scores[, colnames(N2scores) %in% c("n2.ratioFraud", "n2.ratioNonFraud")] =
    apply(N2scores[, colnames(N2scores) %in% c("n2.ratioFraud", "n2.ratioNonFraud")], 2, stdize)

  #### 1.8.2 Combine results ####
  if(printProgress) {
    for(i in 1:3) {
      cat("\rFinishing the simulation", rep(".", i))
      Sys.sleep(0.25)
    }
  }
  unregister_dopar()

  NetwScores = data.table(cbind.data.frame(N1scores, N2scores))
  NetwScores[, nodeID := seq_len(nrow(NetwScores))]

  RmCols = intersect(colnames(DtClaim), colnames(NetwScores)) %>% .[. != "nodeID"]
  if(length(RmCols) > 0)
    DtClaim[, (RmCols) := NULL]
  DtClaimAll = merge(DtClaim, NetwScores, by = "nodeID")
  DtClaimAll[, fraudScore := ResultsNetw$ResultsClaims$ScaledScore[match(nodeID, ResultsNetw$ResultsClaims$ID)]]

  ResultSimulation =
    list(call = call,
         Dt = DtClaimAll,
         SummaryPlots = SummaryPlots,
         TargetPrev = TargetPrev,
         TruePrevalence = sum(DtClaimAll$Fraud == 1) / nrow(DtClaimAll),
         AdjMatOrig = AdjMatOrig,
         AdjMat = AdjMat,
         BiRank = ResultsNetw)
  class(ResultSimulation) = "sfnData"
  return(ResultSimulation)
}

#' Print function for a sfnData object
#'
#' Prints the call, target class imbalance, class imbalance in the synthetic data set, the first part of the data set and
#' the summary ggplot object.
#'
#' @param x an object of class \code{sfnData}, resulting from \code{\link{sfnGenerator}}
#' @param ... arguments passed to \code{\link{print}}
#'
#' @seealso [sfnGenerator()]
#' @method print sfnData
#'
#' @return The original \code{sfnData} object is returned.
#' @export
print.sfnData <- function(x, ...) {
  cat("\nSynthetic social network data set for fraud detection:\n")
  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(" - Target prevalence:", round(x$TargetPrev * 1e2), "%")
  cat("\n - Prevalence in data set:", round(x$TruePrevalence * 1e2), "%\n\n")
  cat("Simulated data set:\n\n")
  print(head(x$Dt), ...)
  plot(x$SummaryPlots)
  invisible(x)
}


