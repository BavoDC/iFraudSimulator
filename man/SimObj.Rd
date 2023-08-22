\name{SimObj}
\alias{SimObj}
\docType{data}
\title{
 An example of a synthetic insurance fraud network data set
}
\description{
 The \code{SimObj} object is an example of a synthetically generated data set using the function \code{\link{sfnGenerator}}. The object was created using the
 code \code{SimObj = sfnGenerator(TargetPrev = 0.05, NrPH = 10000, Seed = 3, printProgress = TRUE, Parallel = TRUE)}.
}
\usage{data("SimObj")}
\value{
The \code{SimObj} object has the following slots:

\item{call}{the matched call.}

\item{Dt}{the synthetic data set.}

\item{SummaryPlots}{a ggplot object with the summary plots of the synthetic data set.}

\item{TargetPrev}{the target class imbalance.}

\item{TruePrevalence}{the class imbalance in the synthetic data set.}

\item{AdjMatOrig}{a \code{data.table} which contains all information to construct the adjacency matrix.}

\item{AdjMat}{the adjacency matrix of the bipartite graph.}

\item{BiRank}{the object resulting from running the BiRank algorithm using the function \code{\link{BiRankFr}}.}


}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
  Campo, Bavo D.C., and Antonio, Katrien (2023). An engine to simulate insurance fraud network data.
}
\examples{
data(SimObj)
str(SimObj)
Homophily(SimObj)
}
\seealso{
\code{\link{sfnGenerator}}
}
\keyword{datasets}
