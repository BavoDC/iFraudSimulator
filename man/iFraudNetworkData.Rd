\name{iFraudNetworkData}
\alias{iFraudNetworkData}
\docType{data}
\title{
  Description of a synthetically generated insurance fraud network data set
}
\description{
 Here, we provide a short description of the variables in the synthetic data set generated using the function \code{\link{sfnGenerator}}. For a detailed description  on how the variables are generated, we refer to the paper on arXiv.
}
\usage{data("iFraudNetworkData")}
\format{
  A data frame with 6948 observations on the following 60 variables.
  \describe{
    \item{\code{nodeID}}{unique ID for the node}
    \item{\code{ClaimID}}{unique ID to identify the claim}
    \item{\code{AgePH}}{age of the policyholder in years}
    \item{\code{GenderPH}}{the gender of the policyholder, encoded as a factor with levels \code{female}, \code{male} and \code{non-binary}}
    \item{\code{ExpPH}}{time since inception of the first contract of the policyholder, in years}
    \item{\code{NrContractsPH}}{number of contracts of the policyholder}
    \item{\code{RateNrContracts}}{rate parameter for generating the number of contracts}
    \item{\code{IDPH}}{unique ID to identify the policyholder}
    \item{\code{ContractID}}{unique ID to identify the contract of the policyholder}
    \item{\code{ExpPHContracts}}{duration or exposure of the contract, in years}
    \item{\code{AgeCar}}{age of the vehicle, in years}
    \item{\code{OrigValueCar}}{original value of the car}
    \item{\code{ValueCar}}{current value of the car}
    \item{\code{ValueCarScaled}}{normalized version of \code{ValueCar}}
    \item{\code{AgeCarScaled}}{normalized version of \code{AgeCarg}}
    \item{\code{AgePHScaled}}{normalized version of \code{AgePH}}
    \item{\code{Coverage}}{type of coverage provided by the insurance company: \code{TPL} (only third party liability),
    \code{PO} (partial omnium = TPL + limited material damage),  \code{FO} (full omnium, TPL + comprehensive material damage)}
    \item{\code{Fuel}}{type of fuel of the vehicle, encoded as a factor with levels \code{Gasoline/LPG/Other} and \code{Diesel}}
    \item{\code{BonusMalus}}{level occupied in bonus-malus scale of the insurance company}
    \item{\code{AgePHBin}}{the binned version of the variable \code{AgePH}: a factor with levels \code{[18,26]} \code{(26,30]} \code{(30,36]} \code{(36,50]} \code{(50,60]} \code{(60,65]} \code{(65,70]} \code{(70,79.9]}}
    \item{\code{AgeCarBin}}{the binned version of the variable \code{AgeCar}: a factor with levels \code{[0,5]} \code{(5,10]} \code{(10,15.3]} \code{(15.3,20]}}
    \item{\code{BonusMalusBin}}{the binned version of the variable \code{BonusMalus}: a factor with levels \code{[0,1)} \code{[1,2)} \code{[2,3)} \code{[3,7)} \code{[7,9)} \code{[9,11)} \code{[11,22]}}
    \item{\code{NClaims}}{the number of claims}
    \item{\code{Garage}}{the garage to which the claim is connected}
    \item{\code{Expert}}{the expert, if any, to which the claim is connected}
    \item{\code{Broker}}{the broker to which the claim is connected}
    \item{\code{nPersons}}{the number of persons, other than the policyholder, to which the claim is connected}
    \item{\code{Police}}{a character vector indicating whether the police was called when the incident}
    \item{\code{ClaimAge}}{the normalized version of the number of months from beginning of contract to the date of the incident}
    \item{\code{ClaimDate}}{number of years between the start of the contract and the claim's filing date}
    \item{\code{TimeSinceClaim}}{the number of years elapsed since the claim was filed}
    \item{\code{Policyholder}}{name of the policyholder}
    \item{\code{Lambda}}{the rate parameter in the gamma distribution for generating the individual claim amounts}
    \item{\code{ClaimAmount}}{the normalized version of the claim amount}
    \item{\code{ClaimAmountOrig}}{the claim amount}
    \item{\code{n1Size}}{the normalized version of the size of the first order neighborhood}
    \item{\code{n2Size}}{the normalized version of the size of the second order neighborhood}
    \item{\code{Criminal}}{a binary variable indicating whether the policyholder is a criminal or not}
    \item{\code{Fraud}}{the ground truth label, indicating whether the claim is fraudulent or non-fraudulent}
    \item{\code{AgePHOrig}}{the age of the policyholder, on its original scale}
    \item{\code{ClaimAgeOrig}}{number of months from beginning of contract to the date of the incident}
    \item{\code{Rule1}}{a logical vector, indicating whether the claim is flagged according to the first business rule}
    \item{\code{Rule2}}{a logical vector, indicating whether the claim is flagged according to the second business rule}
    \item{\code{Rule3}}{a logical vector, indicating whether the claim is flagged according to the third business rule}
    \item{\code{Investigated}}{a binary variable, indicating if the claim is investigated}
    \item{\code{ExpertJudgement}}{a logical vector, indicating if the claim was judged to be fraudulent according to the expert}
    \item{\code{n1.q1}}{the first quartile of the empirical distribution of the fraud scores in the node's first order neighborhood}
    \item{\code{n1.med}}{the median of the empirical distribution of the fraud scores in the node's first order neighborhood}
    \item{\code{n1.max}}{the maximum of the empirical distribution of the fraud scores in the node's first order neighborhood}
    \item{\code{n1.size}}{the number of nodes in the node's first order neighborhood}
    \item{\code{n2.q1}}{the first quartile of the empirical distribution of the fraud scores in the node's second order neighborhood}
    \item{\code{n2.med}}{the median of the empirical distribution of the fraud scores in the node's second order neighborhood}
    \item{\code{n2.max}}{the maximum of the empirical distribution of the fraud scores in the node's second order neighborhood}
    \item{\code{n2.size}}{the number of nodes in the node's second order neighborhood}
    \item{\code{n2.ratioFraud}}{the normalized version of the number of known fraudulent claims in the node's second order neighborhood divided by \code{n2.size}}
    \item{\code{n2.ratioNonFraud}}{the normalized version of the number of known non-fraudulent claims in the node's second order neighborhood divided by \code{n2.size}}
    \item{\code{n2.binFraud}}{a binary value indicating whether there is a known fraudulent claim in the node's second order neighborhood}
    \item{\code{n2.ratioFraudOrig}}{the number of known fraudulent claims in the node's second order neighborhood divided by \code{n2.size}}
    \item{\code{n2.ratioNonFraudOrig}}{the number of known non-fraudulent claims in the node's second order neighborhood divided by \code{n2.size}}
    \item{\code{fraudScore}}{the fraud score as computed by the BiRank algorithm}
  }
}
%%\details{
%%}
%%\source{
%%}
\references{
  Campo, Bavo D.C., and Antonio, Katrien (2023). An engine to simulate insurance fraud network data.
}
\examples{
data(iFraudNetworkData)
str(iFraudNetworkData)
# Class imbalance
sum(iFraudNetworkData$Fraud) / nrow(iFraudNetworkData)
}
\seealso{
\code{\link{sfnGenerator}}
}
\keyword{datasets}
