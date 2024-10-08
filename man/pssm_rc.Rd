% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm-utils.R
\name{pssm_rc}
\alias{pssm_rc}
\title{Reverse complement a PSSM}
\usage{
pssm_rc(pssm)
}
\arguments{
\item{pssm}{A PSSM. Data frame with columns 'A', 'C', 'G', 'T' and 'pos' or a matrix with columns 'A', 'C', 'G', 'T'}
}
\value{
A PSSM with the same format, but reverse complemented.
}
\description{
Reverse complement a PSSM
}
\examples{
# Create simulated PSSM data frame
pssm <- data.frame(
    pos = 1:4,
    A = c(0.1, 0.2, 0.3, 0.1),
    C = c(0.1, 0.3, 0.2, 0.1),
    G = c(0.1, 0.3, 0.3, 0.7),
    T = c(0.7, 0.2, 0.2, 0.1)
)

# Reverse complement the PSSM
rc_pssm <- pssm_rc(pssm)

}
