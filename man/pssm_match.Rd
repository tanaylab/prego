% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm-utils.R
\name{pssm_match}
\alias{pssm_match}
\title{Match PSSM to a directory of motifs}
\usage{
pssm_match(
  pssm,
  motifs,
  best = FALSE,
  method = "spearman",
  parallel = getOption("prego.parallel", TRUE)
)
}
\arguments{
\item{pssm}{PSSM matrix or data frame}

\item{motifs}{a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name}

\item{best}{return the best match only}

\item{method}{method to use for computing the correlation. See \code{\link[stats]{cor}} for details.}

\item{parallel}{use parallel processing. Set the number of cores using \code{set_parallel}.}
}
\value{
if \code{best} is \code{TRUE}, a string with the best match. Otherwise, a data frame with a row per motif and a column named 'cor' with its correlation to \code{pssm}. The data frame is sorted by descreasing correlation.
}
\description{
Match a PSSM to a directory of motifs. The PSSM is matched to each motif in the directory by computing the correlation between the two PSSMs.
}
\examples{
\dontrun{
res1 <- regress_pwm(cluster_sequences_example, cluster_mat_example[, 1])
head(pssm_match(res1$pssm, JASPAR_motifs))
pssm_match(res1$pssm, JASPAR_motifs, best = TRUE)
}

}
