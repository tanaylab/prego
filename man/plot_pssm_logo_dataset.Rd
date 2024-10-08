% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pssm-utils.R
\name{plot_pssm_logo_dataset}
\alias{plot_pssm_logo_dataset}
\title{Plot LOGO of pssm from dataset (e.g. "HOMER" or "JASPAR")}
\usage{
plot_pssm_logo_dataset(
  motif,
  dataset = all_motif_datasets(),
  title = motif,
  subtitle = ggplot2::waiver(),
  pos_bits_thresh = NULL,
  revcomp = FALSE,
  method = "bits"
)
}
\arguments{
\item{motif}{the motif name (e.g. "GATA4")}

\item{dataset}{a data frame with PSSMs ('A', 'C', 'G' and 'T' columns), with an additional column 'motif' containing the motif name, for example \code{HOMER_motifs}, \code{JASPAR_motifs} or all_motif_datasets()}

\item{title}{title of the plot}

\item{subtitle}{subtitle of the plot}

\item{pos_bits_thresh}{Positions with bits above this threshold would be highlighted in red. If \code{NULL}, no positions would be highlighted.}

\item{revcomp}{whether to plot the reverse complement of the PSSM}

\item{method}{Height method, can be one of "bits" or "probability" (default:"bits")}
}
\value{
a ggplot object
}
\description{
Plot LOGO of pssm from dataset (e.g. "HOMER" or "JASPAR")
}
\examples{

plot_pssm_logo_dataset("JASPAR.Brachyury")

plot_pssm_logo_dataset("GATA5", JASPAR_motifs)

}
