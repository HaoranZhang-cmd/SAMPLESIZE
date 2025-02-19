% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.1.2.R
\name{samplearea}
\alias{samplearea}
\title{Sample Size for Estimating the Area Under the ROC Curve}
\usage{
samplearea(A, R, alpha, beta, L)
}
\arguments{
\item{A}{the conjectured AUC of the test}

\item{R}{the ratio of the number of patients without the condition to patients with the condition in the study sample}

\item{alpha}{1-alpha is the size of the desired confidence interval}

\item{beta}{1-beta is the desired power}

\item{L}{the desired width of one-half of the CI}
}
\value{
n.normal is a more precise estimation of the number of patients with condition
        n.rough is a rough estimation of the number of patients with condition
        N.normal is a more precise estimation of the total number of patients
        N.rough is a rough estimation of the total number of patients
}
\description{
Sample Size for Estimating the Area Under the ROC Curve
}
