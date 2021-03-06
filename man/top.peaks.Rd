% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{top.peaks}
\alias{top.peaks}
\title{Returns the top peaks from sorted data frame}
\usage{
top.peaks(peak.df, npeaks = min(100, nrow(peak.df)))
}
\arguments{
\item{peak.df}{Peak data frame}

\item{npeaks}{Number of peaks to return (default 100)}
}
\value{
A subset of the peak data frame with the specified number of peaks
}
\description{
Returns a specified number of peaks from the top of the peak data frame.
If data frame doesn't have that many rows, then returns entire data frame.
}
