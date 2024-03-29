% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_information4.R
\name{score-information}
\alias{score-information}
\alias{getScore}
\alias{getScore.rcon}
\alias{getScore.rcor}
\title{Get score and information for RCOX models.}
\usage{
getScore(m, K, scale = "original")

\method{getScore}{rcon}(m, K, scale = "original")

\method{getScore}{rcor}(m, K, scale = "original")
}
\arguments{
\item{m}{Model obect.}

\item{K}{Concentration matrix}

\item{scale}{On which scales are the parameters computed.}
}
\description{
Get score and information for RCOX models.
}
\author{
Søren Højsgaard, \email{sorenh@math.aau.dk}
}
