% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{roc.auc}
\alias{roc.auc}
\title{Receiver operating curve for peak classification}
\usage{
roc.auc(peak.df, true.pos.peaks, make.plot = TRUE, title.text = "")
}
\arguments{
\item{peak.df}{a sorted peak data frame
Note: peak data frame must be arranged so that most important peaks are on top
 of the data frame. Least important peaks on the bottom.}

\item{true.pos.peaks}{a vector of peak IDs considered to be true positive peaks}

\item{make.plot}{if true, plots ROC curve. Otherwise, returns area under ROC curve.}

\item{title.text}{the title of the the ROC curve. Only affects result if make.plot is TRUE.}
}
\value{
Plot of ROC curve or area under ROC curve
}
\description{
Plots accuracy of peak ranking against known true positive peaks as
receiver operating curve (ROC) plot.
Alternatively, can return area under the ROC curve (AUC).
}
