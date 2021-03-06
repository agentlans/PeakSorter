% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\docType{package}
\name{PeakSorter-package}
\alias{PeakSorter}
\alias{PeakSorter-package}
\title{Tools for processing LC-MS peak lists.}
\description{
PeakSorter provides functions for sorting LC-MS peak lists
to prioritize peaks of interest.
}
\details{
Use \code{\link{read.peak.file}} to read peak lists from files
and \code{\link{bin.prioritize}} to sort the peak list.
If you have a list of true positive peaks, you can use \code{\link{roc.auc}}
to check the accuracy of the sort.
}
