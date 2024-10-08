% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm-utils.R
\name{gextract.local_pwm}
\alias{gextract.local_pwm}
\title{Extracts local position weight matrix (PWM) scores for given intervals and a PWM.}
\usage{
gextract.local_pwm(
  intervals,
  pssm,
  spat = NULL,
  spat_min = 0,
  spat_max = NULL,
  bidirect = TRUE,
  prior = 0.01
)
}
\arguments{
\item{intervals}{The intervals to extract}

\item{pssm}{a PSSM matrix or data frame. The columns of the matrix or data frame should be named with the nucleotides ('A', 'C', 'G' and 'T').}

\item{spat}{a data frame with the spatial model (as returned from the \code{$spat} slot from the regression). Should contain a column called 'bin' and a column called 'spat_factor'.}

\item{spat_min}{the minimum position to use from the sequences. The default is 1.}

\item{spat_max}{the maximum position to use from the sequences. The default is the length of the sequences.}

\item{bidirect}{is the motif bi-directional. If TRUE, the reverse-complement of the motif will be used as well.}

\item{prior}{a prior probability for each nucleotide.}
}
\value{
A matrix with \code{nrow(intervals)} rows and \code{ncol(pssm)} columns with the local PWM for each sequence in each position.
}
\description{
Extracts local position weight matrix (PWM) scores for given intervals and a PWM.
}
