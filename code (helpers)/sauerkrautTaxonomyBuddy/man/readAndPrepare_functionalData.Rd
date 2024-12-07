% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readAndPrepare_functionalData.R
\name{readAndPrepare_functionalData}
\alias{readAndPrepare_functionalData}
\title{Prepare the functional datasets}
\usage{
readAndPrepare_functionalData(
  path_KEGG_KO,
  path_KEGG_pathways,
  path_sampleLookupXlsx,
  path_krautLookupXlsx,
  abundance_unit = "TPM",
  exclude_medicatedParticipants = TRUE,
  exclude_sickParticipants
)
}
\arguments{
\item{path_KEGG_KO}{Path to the file \code{KO-sum-wide.txt}}

\item{path_KEGG_pathways}{Path to the file \code{kegg_otu_table.txt}}

\item{path_sampleLookupXlsx}{Path to the \code{Metabolon_Blutproben_Stuhlproben_intern.xlsx}
file which contains one sheet \code{Sample Manifest} summarizing all relevant
sample information, including the information which stool and blood samples were
analyzed twice in the lab for us to measure the measurement variability.}

\item{path_krautLookupXlsx}{Path to the \code{Sauerkrautproben für Sequenzierung.xlsx}
file which contains lookup information on the sauerkraut samples.}

\item{abundance_unit}{One of \code{c("TPM", "relAbd")}, specifying if the
returned abundances should be in unit TPM (transcripts per million reads
mapped) or returned as relative abundances (i.e., TPM divided by one million).}

\item{exclude_medicatedParticipants}{Indicator if the two participants SK005
and SK093 that were on medication during the whole study should be removed
from the data. Defaults to TRUE.}

\item{exclude_sickParticipants}{One of \code{c("no exclusion", "stool type analysis", "taxonomy paper", "blood paper")}.
If one of the latter is specified, all measurements biased by a recent sickness
are excluded. Since the taxonomy paper (mainly stool measurements) and the blood paper
(blood measurements) analyze different biosamples, the exclusion list differs
between the two papers. For the stool type analysis, a few more measurements
are excluded compared to the taxonomy paper.

Note that medicated or sick participants' measurements' information is not
excluded from the returned blood Boston duplicate data since this is only
a secondary lookup table used for matching original with duplicate measurements
measured by the Boston laboratory, and not used for filtering other data.}
}
\value{
List of data.frames
}
\description{
This function reads the KEGG-KO and KEGG-pathway datasets and prepares them.
}
