% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.1.6.R
\name{samplepartialarea}
\alias{samplepartialarea}
\title{Sample Size for Estimating the Partial Area Under the ROC Curve}
\usage{
samplepartialarea(a, b, e1, e2, alpha, beta, L, R)
}
\arguments{
\item{a}{binormal parameters}

\item{b}{binormal parameters}

\item{e1}{lower bound for partial area}

\item{e2}{upper bound for partial area}

\item{alpha}{significance level}

\item{beta}{1-beta is the power of the test}

\item{L}{the desired width of one-half of the CI}

\item{R}{ratio of undiseased subjects and diseased subjects}
}
\value{
n_one sample size of diseased subjects in a one-sided hypothesis testing
        N_one total sample size in a one-sided hypothesis testing
        n_two sample size of diseased subjects in a two-sided hypothesis testing
        N_two total sample size in a two-sided hypothesis testing
}
\description{
Sample Size for Estimating the Partial Area Under the ROC Curve
}
