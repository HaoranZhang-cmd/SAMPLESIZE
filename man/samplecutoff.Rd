% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.4.R
\name{samplecutoff}
\alias{samplecutoff}
\title{SAMPLE SIZE FOR DETERMINING A SUITABLE CUTOFF VALUE}
\usage{
samplecutoff(SP, SE_con, SE_min, b, alpha, beta)
}
\arguments{
\item{SP}{the minimum specificity required at the cutoff}

\item{SE_con}{the conjectured sensitivity at the cutoff}

\item{SE_min}{the sensitivity at which the cutoff value for the test would no longer be useful}

\item{b}{the binormal parameter b for the test result}

\item{alpha}{The choice of alpha determines an overall confidence probability of 1-alpha for the statement "the specificity at the chosen cutoff is at least SP"}

\item{beta}{The choice of 1-beta determines the probability of finding such a cutoff, when such a cutoff truly exists}
}
\value{
N the rough total sample size when the cutoff point is not determined based on a prespecified specificity and sensitivity
        N_pre the precise total sample size when the cutoff point is determined based on a prespecified specificity and sensitivity
}
\description{
SAMPLE SIZE FOR DETERMINING A SUITABLE CUTOFF VALUE
}
