#' Tools for processing LC-MS peak lists.
#'
#' PeakSorter provides functions for sorting LC-MS peak lists
#' to prioritize peaks of interest.
#'
#' Use \code{\link{read.peak.file}} to read peak lists from files
#' and \code{\link{bin.prioritize}} to sort the peak list.
#' If you have a list of true positive peaks, you can use \code{\link{roc.auc}}
#' to check the accuracy of the sort.
"_PACKAGE"

#' Idh1 knock-in mouse dataset
#'
#' Peak data frame containing peaks from untargeted LC-MS profiling
#' of Idh1 knock-in mouse plasma on C18 column negative mode.
#' For more details, please see (unpublished paper/thesis)
#'
#' @format A data frame with 5000 rows and 4 variables:
#' \describe{
#' \item{Rt}{Retention time in minutes}
#' \item{Mz}{Mass to charge ratio}
#' \item{Intensity}{Normalized intensity}
#' \item{Name}{Unique peak ID}
#' }
#' @seealso \code{\link{idh1.true.pos}}
"idh1"

#' True positive peaks in the Idh1 knock-in metabolomics dataset
#'
#' A list of peak IDs from the Idh1 knock-in metabolomics experiment
#' whose m/z correspond to metabolites in the HMDB database.
#' For more details, please see (unpublished paper/thesis)
#'
#' @seealso \code{\link{idh1}}
"idh1.true.pos"

#' Read peak data from file
#'
#' Reads peak data from given text file
#'
#' @param filename The file containing peak data
#' (exactly four tab-separated columns:
#' retention time (min.), m/z, intensity, peak name)
#' No header row.
#' @return data frame of peaks
#' @export
#' @importFrom utils read.table
read.peak.file <- function(filename) {
  # Read the sorted peaks
  sorted.peaks.df <- utils::read.table(
    filename, sep="\t", header=FALSE)
  colnames(sorted.peaks.df) <-
    c("Rt", "Mz", "Intensity", "Name")
  check.peak.df(sorted.peaks.df)
  sorted.peaks.df
}

#' Validate peak data frame
#'
#' Checks peak data frame to make sure it's valid
#'
#' @param peak.df Peak data frame
#' @return TRUE if OK otherwise FALSE
#' @export
check.peak.df <- function(peak.df) {
  if (!is.data.frame(peak.df)) {
    stop("Need valid peak data frame.")
    return(FALSE)
  } else if (!all(c("Rt","Mz","Intensity","Name") %in% colnames(peak.df))) {
    stop("Peak data frame must include columns named Rt, Mz, Intensity, and Name.")
    return(FALSE)
  }
  if (nrow(peak.df) == 0) {
    warning("There are no peaks in the dataset.")
  }
  if (nrow(peak.df) == 1) {
    warning("There is only one peak in the dataset.")
  }
  if (any(duplicated(peak.df$Name))) {
    warning("Duplicate names in peak data frame.")
    return(FALSE)
  }
  return(TRUE)
}

#' Scrambles the values in each column of peak data frame
#'
#' Given a peak data frame, scrambles specified columns separately
#' and returns the scrambled version
#'
#' @param peak.df Peak data frame
#' @param column.name The name of column to scramble
#' @return The data frame but with the scrambled column scrambled
#' @export
scramble.col <- function(peak.df, column.name) {
  check.peak.df(peak.df)
  # Case where there's only 0 or 1 peak
  if (nrow(peak.df) <= 1) {
    return(peak.df)
  }
  temp <- peak.df
  if (!all(column.name %in% colnames(peak.df))) {
    stop("Column names are not in original peak data frame.")
  }
  if (length(column.name) == 1) {
    # Scramble single column
    temp[,column.name] <- sample(temp[,column.name])
  } else if (length(column.name) > 1) {
    # Scramble each column separately as specified in input
    for (col.id in column.name) {
      temp[,col.id] <- sample(temp[,col.id])
    }
  } else {
    stop("Error. Must have column to scramble.")
  }
  temp
}

#' Sorts peak data frame by intensity
#'
#' Given a peak data frame, reorders data frame
#' by peak intensity (maximum to minimum intensity). Note: rownames will be removed in output.
#'
#' @param peak.df Peak data frame
#' @return The peak data frame arranged in descending order of intensity
#' @seealso \code{\link{bin.prioritize}}, \code{\link{random.prioritize}}
#' @export
intensity.prioritize <- function(peak.df) {
  check.peak.df(peak.df)
  # Case where there's only 0 or 1 peak in dataset
  if (nrow(peak.df) <= 1) {
    return(peak.df)
  }
  temp <- peak.df[order(-peak.df$Intensity),]
  rownames(temp) <- NULL
  temp
}

#' Randomly change order of peaks in data frame
#'
#' Returns peaks in random order. Note: rownames will be removed in output.
#'
#' @param peak.df Peak data frame
#' @return Data frame with rows in random order
#' @seealso \code{\link{bin.prioritize}}, \code{\link{intensity.prioritize}}
#' @export
random.prioritize <- function(peak.df) {
  check.peak.df(peak.df)
  # Case where there's only 0 or 1 peak in dataset
  if (nrow(peak.df) <= 1) {
    return(peak.df)
  }
  temp <- peak.df[sample(1:nrow(peak.df)),]
  rownames(temp) <- NULL
  temp
}

#' Finds rank of values in descending order
#'
#' Ranking function called by various prioritization methods
#'
#' Orders values in descending order
#'
#' @param x a vector to be sorted
#' @return sorted vector
rank.fun <- function(x) rank(-x, ties.method="first")

#' Binned sorting method
#'
#' Reorders peaks using the binned method.
#'
#' Binned method works as follows:
#' \enumerate{
#' \item Retention time is divided into bins of specified binwidth.
#' \item Most intense peaks from each bin are arranged in order of decreasing intensity.
#' \item Second most intense peaks from each bin are arranged in order of decreasing intensity.
#' \item Process repeated until all peaks are sorted.
#' }
#' Note: in this implementation, the rownames of original peak data frame will be removed in output.
#'
#' @param peak.df Peak data frame to sort
#' @param binwidth Width of the bins from which to pick peaks
#' @return Data frame of peaks
#' @seealso \code{\link{intensity.prioritize}}, \code{\link{random.prioritize}}
#' @export
#' @importFrom stats ave
bin.prioritize <- function(peak.df, binwidth=0.01) {
  check.peak.df(peak.df)
  # Case where there's only 0 or 1 peak in dataset
  if (nrow(peak.df) <= 1) {
    return(peak.df)
  }
  # Save copy of peak data frame
  temp <- peak.df
  # and copy of rownames
  #rowname.vec <- rownames(peak.df)
  #names(rowname.vec) <- peak.df$Name
  # Bin the retention time.
  if (binwidth <= 0) {
    stop("Bin width for the binned method must be positive.")
  }
  rt.range <- max(peak.df$Rt) - min(peak.df$Rt)
  if (rt.range / binwidth >= 2) {
    peak.df$RtGroup <- cut(peak.df$Rt, floor(rt.range/binwidth))
  } else {
    # If binwidth is too large, then put every peak in one bin
    peak.df$RtGroup = 1
  }
  # Assign ranking within each group
  peak.df$BinRank <- ave(peak.df$Intensity, peak.df$RtGroup, FUN=rank.fun)
  # Assign ranking for peaks with same rank in each group
  peak.df$IntensityRank <- ave(peak.df$Intensity, peak.df$BinRank, FUN=rank.fun)
  # Select the peak columns and remove row names
  temp3 <- peak.df[order(peak.df$BinRank, peak.df$IntensityRank),
                   c("Rt", "Mz", "Intensity", "Name")]
  rownames(temp3) <- NULL
  temp3
}

#' Returns the top peaks from sorted data frame
#'
#' Returns a specified number of peaks from the top of the peak data frame.
#' If data frame doesn't have that many rows, then returns entire data frame.
#'
#' @param peak.df Peak data frame
#' @param npeaks Number of peaks to return (default 100)
#' @return A subset of the peak data frame with the specified number of peaks
#' @export
top.peaks <- function(peak.df, npeaks=min(100, nrow(peak.df))) {
  check.peak.df(peak.df)
  if (nrow(peak.df) <= 1) {
    return(peak.df)
  } else {
    peak.df[1:npeaks,]
  }
}

#' Dwell times of peaks in data frame
#'
#' Finds the dwell time of a data frame containing peak data
#' (defined here as the retention time difference
#' between successive peaks in the data frame).
#'
#' @param peak.df Peak data frame
#' @return A vector of retention time differences
#' @export
dwell.times <- function(peak.df) {
  check.peak.df(peak.df)
  if (nrow(peak.df) >= 2) {
    # Find difference between successive retention times
    rt <- sort(peak.df$Rt)
    rt[2:length(rt)] - rt[1:(length(rt)-1)]
  } else {
    warning("Should have at least 2 peaks to calculate retention time. Returning NA.")
    NA
  }
}

#' Receiver operating curve for peak classification
#'
#' Plots accuracy of peak ranking against known true positive peaks as
#' receiver operating curve (ROC) plot.
#' Alternatively, can return area under the ROC curve (AUC).
#'
#' @param peak.df a sorted peak data frame
#' Note: peak data frame must be arranged so that most important peaks are on top
#'  of the data frame. Least important peaks on the bottom.
#' @param true.pos.peaks a vector of peak IDs considered to be true positive peaks
#' @param make.plot if true, plots ROC curve. Otherwise, returns area under ROC curve.
#' @param title.text the title of the the ROC curve. Only affects result if make.plot is TRUE.
#' @return Plot of ROC curve or area under ROC curve
#' @export
#' @importFrom ROCR prediction performance
#' @importFrom grDevices rainbow
#' @importFrom graphics abline mtext plot
roc.auc <- function(peak.df, true.pos.peaks, make.plot=TRUE, title.text="") {
  check.peak.df(peak.df)
  # Validate the true positive peaks
  if (!all(true.pos.peaks %in% peak.df$Name)) {
    stop("Some true positive peaks aren't in original dataset!")
  }
  if (any(duplicated(true.pos.peaks))) {
    warning("Some true positive peaks are duplicated. Those will be removed from dataset.")
    true.pos.peaks <- unique(true.pos.peaks)
  }
  # Set up the ROC data
  pred <- ROCR::prediction(
    # So peaks in beginning of data frame have higher value than peaks at end
    -(1:nrow(peak.df)),
    peak.df$Name %in% true.pos.peaks)
  # Find AUC
  auc <- ROCR::performance(pred, measure="auc")@y.values[[1]]
  # Make ROC plot
  if (make.plot) {
    perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
    ROCR::plot(perf, col=grDevices::rainbow(10), main=title.text)
    graphics::abline(0,1)
    graphics::mtext(paste("AUC =", signif(auc, 3)))
  } else {
    # If not, then make AUC
    auc
  }
}
