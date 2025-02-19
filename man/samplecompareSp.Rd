% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.2.2.R
\name{samplecompareSp}
\alias{samplecompareSp}
\title{Sample Size for Comparing Tests Sensitivity and/or Specificity}
\usage{
samplecompareSp(Sp1, Sp2, p, alpha, beta, R)
}
\arguments{
\item{Sp1}{the conjectured value of sensitivity under the alternative hypothesis of test 1}

\item{Sp2}{the conjectured value of sensitivity under the alternative hypothesis of test 2}

\item{p}{P(T1=0|T2=0)}

\item{alpha}{significance level}

\item{beta}{1-beta is the power of the test}

\item{R}{the ratio of undiseased subjects and diseased subjects}
}
\value{
n_one sample size of diseased subjects in a one-sided hypothesis testing
        N_one total sample size in a one-sided hypothesis testing
        n_two sample size of diseased subjects in a two-sided hypothesis testing
        N_two total sample size in a two-sided hypothesis testing
}
\description{
Sample Size for Comparing Tests Sensitivity and/or Specificity
}
