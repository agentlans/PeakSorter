
#' Validate a list of peak data frames
#'
#' Checks whether all objects in list are valid peak data frames
#'
#' @param peak.df.list list of peak data frames to test
#' @return TRUE if all data frames are valid. Otherwise, FALSE
#' @seealso \code{\link{check.peak.df}}
validate.peak.df.list <- function(peak.df.list) {
  if (!is.list(peak.df.list)) {
    stop("Can only plot a list of peak data frames.")
  }
  if (!all(sapply(peak.df.list, check.peak.df))) {
    stop("Not all the peak data frames in list have right format.")
  }
  return(TRUE)
}

#' Wrapper for boxplot function that also prints P values
#'
#' @param y.list a named list of vectors indicating values for each group
#' @param ... additional arguments to the boxplot function
#' @return a boxplot with P values
#' @examples
#'     boxplot2(list(`Foo`=1:3, `Bar`=10:13))
#' @export
boxplot2 <- function(y.list, ...) {
  if (length(y.list) == 0) {
    stop("Must have at least one group to make boxplot.")
  }
  if (is.null(names(y.list))) {
    names(y.list) <- paste("Group", 1:length(y.list))
  }
  # Convert the list into a data frame of (x, y) pairs
  plottable <- do.call(rbind, lapply(1:length(y.list), function(i) {
    data.frame(
      x=rep(names(y.list)[i], length(y.list[[i]])),
      y=y.list[[i]])
  }))
  # Create boxplot
  boxplot(y~x, data=plottable, ...)
  if (length(y.list) >= 2) {
    pval <- kruskal.test(y~x, data=plottable)$p.value
    mtext(paste("P = ", signif(pval, 2), ", Kruskal-Wallis rank sum test", sep=""))
  }
}

#' Boxplots of dwell times of multiple peak lists
#'
#' Plots a boxplot of dwell times of different peak lists
#'
#' @param peak.df.list a list of peak data frames
#' @param xlab label for X axis
#' @param ylab label for y axis
#' @param ... additional parameters for boxplot
#' @return boxplot of dwell times where each box corresponds to each peak data frame in list
#' @export
dwell.time.boxplot <- function(peak.df.list,
                               xlab="Method", ylab="Dwell time (min.)",
                               ...) {
  validate.peak.df.list(peak.df.list)
  dwell.times <- lapply(peak.df.list, dwell.times)
  boxplot2(dwell.times, xlab=xlab, ylab=ylab, ...)
}

#' Generates graphable data from a list of peak data frames
#'
#' Note: intended to be called by other graphing functions
#'
#' @param peak.df.list a list of peak data frames. Names indicate the method to be labelled on plot.
#' @return a data frame for graphing
peak.plot.df <- function(peak.df.list) {
  # Create names if there aren't any
  if (is.null(names(peak.df.list))) {
    names(peak.df.list) <- paste("Method", 1:length(peak.df.list))
  }
  # Validate list of peak data frames
  validate.peak.df.list(peak.df.list)
  # For each peak data frame, collect information as graphable data frame
  do.call(rbind, lapply(names(peak.df.list), function(df.name) {
    peak.df <- peak.df.list[[df.name]]
    cbind(data.frame(Method=rep(df.name, nrow(peak.df))),
          peak.df)
  }))
}

#' Creates plots of peaks across multiple sorted peak data frames
#'
#' Generates plots of the peaks in peak lists
#'
#' @param peak.df.list a list of peak data frames
#' @return plots showing intensity vs. retention time of each peak
#' @export
#' @importFrom ggplot2 ggplot geom_segment facet_grid scale_x_continuous scale_y_log10 aes
#' @importFrom graphics boxplot
#' @importFrom stats kruskal.test
peak.plot <- function(peak.df.list) {
  # Format data properly
  peak.data.df <- peak.plot.df(peak.df.list)
  # Make the plots
  ggplot(peak.data.df, aes(peak.data.df$Rt, peak.data.df$Intensity)) +
    geom_segment(aes(x=peak.data.df$Rt,
                     y=min(peak.data.df$Intensity),
                     xend=peak.data.df$Rt,
                     yend=peak.data.df$Intensity)) +
    facet_grid(Method~.) +
    scale_x_continuous("Retention time (min.)") +
    scale_y_log10("Normalized intensity (normalized counts)")
}

#' Receiver operating characteristic curves of multiple sorted peak data frames
#'
#' Make multiple receiver operating characteristic curves
#'
#' @param peak.df.list a list of peak data frames
#' @param true.pos.peaks a vector containing names of true positive peaks
#' @return plot of multiple ROC curves
#' @seealso \code{\link{roc.auc}}
#' @export
roc.plot <- function(peak.df.list, true.pos.peaks) {
  validate.peak.df.list(peak.df.list)
  # Make new names if there aren't any
  if (is.null(names(peak.df.list))) {
    names(peak.df.list) <- paste("Method", 1:length(peak.df.list))
  }
  # Plot ROC curves for each peak list
  lapply(1:length(peak.df.list), function(i) {
    roc.auc(peak.df.list[[i]], true.pos.peaks,
            title.text=names(peak.df.list)[i])
  })
  NULL
}

