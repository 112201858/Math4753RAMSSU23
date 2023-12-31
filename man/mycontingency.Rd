% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mycontingency.R
\name{mycontingency}
\alias{mycontingency}
\title{mycontingency - Chi-squared test for contingency tables}
\usage{
mycontingency(table, alpha = 0.05, yates = TRUE, verbose = TRUE, visual = TRUE)
}
\arguments{
\item{table}{an rXc contingency table}

\item{alpha}{A confidence coefficient (default = .05)}

\item{yates}{A Boolean for whether or not to use Yates Correction (default = TRUE)}

\item{verbose}{A Boolean for whether or not to print extra information (default = TRUE)}

\item{visual}{A Boolean for whether or not to construct plots (default = TRUE)}
}
\value{
an invisible list called 'mc' that contains various tables,
an exhaustive list of all the values related to the chi-squared test
}
\description{
mycontingency - Chi-squared test for contingency tables
}
\examples{
mycontingency(table = table(c(1,2),c(1,2)))
}
