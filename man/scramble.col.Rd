% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{scramble.col}
\alias{scramble.col}
\title{Scrambles the values in each column of peak data frame}
\usage{
scramble.col(peak.df, column.name)
}
\arguments{
\item{peak.df}{Peak data frame}

\item{column.name}{The name of column to scramble}
}
\value{
The data frame but with the scrambled column scrambled
}
\description{
Given a peak data frame, scrambles specified columns separately
and returns the scrambled version
}
