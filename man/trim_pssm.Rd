% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm-utils.R
\name{trim_pssm}
\alias{trim_pssm}
\alias{pssm_trim}
\title{Trim PSSM}
\usage{
trim_pssm(pssm, bits_thresh = 0.1)

pssm_trim(pssm, bits_thresh = 0.1)
}
\arguments{
\item{pssm}{A data frame representing the PSSM, with columns for position (pos) and bits per position (bits).}

\item{bits_thresh}{The threshold value for bits per position. Positions with bits above this threshold will be kept, while positions with bits below this threshold at the beginning and the end of the motif will be removed. The default value is 0.1.}
}
\value{
A trimmed PSSM data frame, with positions filtered based on the bits threshold.
}
\description{
This function trims a Position-Specific Scoring Matrix (PSSM) by removing positions with low information content at the beginning and end of the motif.
}
