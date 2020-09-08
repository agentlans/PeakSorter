library(PeakSorter)
context("Sorting methods")

# Toy data
peak.df <- data.frame(
  Rt=c(0.1, 0.2, 0.3, 0.6, 0.7, 1.1),
  Mz=c(150, 151, 152, 140, 141, 133),
  Intensity=c(1,3,2,3,2,4),
  Name=c("P1","P2","P3","P4","P5","P6")
)

# No peaks
zero.peak.df <- peak.df[0,]

# Just one peak
one.peak.df <- peak.df[1,]

# Expected output if peak.df sorted with bin size 0.5
# With this bin size, there will be 2 retention time bins:
# one from (0 to 0.6] and one from (0.6, 1.1]
# P1-P4 will be in first bin, P5-P6 will be in second.
peak.df.sorted <- data.frame(
  Rt=c(1.1, 0.2, 0.6, 0.7, 0.3, 0.1),
  Mz=c(133, 151, 140, 141, 152, 150),
  Intensity=c(4,3,3,2,2,1),
  Name=c("P6","P2","P4","P5","P3","P1")
)

# Expected output if peak.df sorted with bin size 0.3
# For this bin size, retention time bins are
# (0, 0.433], (0.433, 0.767], (0.767, 1.1]
# Peaks will be binned as
# {P1, P2, P3}, {P4, P5}, {P6}
peak.df.sorted2 <- data.frame(
  Rt=c(1.1, 0.2, 0.6, 0.3, 0.7, 0.1),
  Mz=c(133, 151, 140, 152, 141, 150),
  Intensity=c(4,3,3,2,2,1),
  Name=c("P6","P2","P4","P3","P5","P1")
)

test_that("Methods works with 0 peaks", {
  expect_equal(random.prioritize(zero.peak.df), zero.peak.df)
  expect_equal(intensity.prioritize(zero.peak.df), zero.peak.df)
  expect_equal(bin.prioritize(zero.peak.df), zero.peak.df)
})

test_that("Methods work with 1 peak", {
  expect_equal(random.prioritize(one.peak.df), one.peak.df)
  expect_equal(intensity.prioritize(one.peak.df), one.peak.df)
  expect_equal(bin.prioritize(one.peak.df), one.peak.df)
})

# Function to reorder peak list by peak name
reorder.by.name <- function(peak.df) {
  temp <- peak.df[order(peak.df$Name),]
  rownames(temp) <- NULL # Remove row names
  temp
}
reordered.peak.df <- reorder.by.name(peak.df)
test_that("Sorted peak lists can be reordered to match original data", {
  expect_equal(reorder.by.name(random.prioritize(peak.df)), reordered.peak.df)
  expect_equal(reorder.by.name(intensity.prioritize(peak.df)), reordered.peak.df)
  expect_equal(reorder.by.name(bin.prioritize(peak.df)), reordered.peak.df)
})

test_that("Intensity method works", {
  expect_equal(intensity.prioritize(peak.df), peak.df.sorted2)
})

# For the following tests, we compare the sequence of intensities only.
# It's possible that the top peaks in different bins can have the same intensity.
# In that case, there is no a priori reason to rank one peak over the other.
test_that("Binned method works with different bin widths", {
  # Default settings
  expect_equal(bin.prioritize(peak.df, 0.5)$Intensity, peak.df.sorted$Intensity)
  # Smaller bin width
  expect_equal(bin.prioritize(peak.df, 0.3)$Intensity, peak.df.sorted2$Intensity)
  # One big bin. Row names don't match.
  expect_equal(bin.prioritize(peak.df, 1)$Intensity, 
               intensity.prioritize(peak.df)$Intensity)
})

test_that("Area under ROC curve calculation works", {
  expect_equal(roc.auc(peak.df, c("P1", "P2", "P3"), make.plot=FALSE), 1)
  expect_equal(roc.auc(peak.df, c("P4", "P5", "P6"), make.plot=FALSE), 0)
})

test_that("Dwell time calculation is right", {
  expect_equal(dwell.times(zero.peak.df), NA)
  expect_equal(dwell.times(one.peak.df), NA)
  expect_equal(dwell.times(peak.df[1:2,]), 0.1)
  expect_equal(dwell.times(random.prioritize(peak.df)), c(0.1, 0.1, 0.3, 0.1, 0.4))
})

test_that("Top peaks are chosen correctly", {
  expect_equal(top.peaks(zero.peak.df), zero.peak.df)
  expect_equal(top.peaks(one.peak.df), one.peak.df)
  expect_equal(top.peaks(peak.df, 3), peak.df[1:3,])
  expect_equal(top.peaks(peak.df), peak.df)
})


# Test whether x is permutation of y (and that x != y)
is.perm <- function(x,y) {
  (!identical(x,y)) && identical(sort(x), sort(y))
}
set.seed(12345)
# Single column shuffling
test_that("Single column shuffling works", {
  shuffled <- scramble.col(peak.df, "Mz")
  # Every other column should be same
  expect_equal(peak.df[,colnames(peak.df) != "Mz"],
               shuffled[,colnames(peak.df) != "Mz"])
  # Only Mz column should be scrambled
  expect_true(is.perm(peak.df$Mz, shuffled$Mz))
})

# Multiple column shuffling
test_that("Multiple column shuffling works", {
  shuffled <- scramble.col(peak.df, c("Mz", "Intensity"))
  # Every other column should be the same
  expect_equal(peak.df[,!colnames(peak.df) %in% c("Mz", "Intensity")],
               shuffled[,!colnames(shuffled) %in% c("Mz", "Intensity")])
  # The columns should be permuted
  expect_true(is.perm(peak.df$Mz, shuffled$Mz))
  expect_true(is.perm(peak.df$Intensity, shuffled$Intensity))
})

