# PeakSorter
Peaksorter is a R package for sorting LC-MS peak lists to prioritize peaks of interest.

Example uses:
- to reduce datasets for sample classification
- to select peaks for further investigation
- to evaluate different peak sorting methods

# Installation
1. Download and install R version 2.10 or higher. [Link](https://www.r-project.org/)
2. In R, run the following code to install prerequisites and this package:
```r
# Install prerequisite packages
install.packages(c('ggplot2', 'ROCR'))

# Install PeakSorter from GitHub
library(devtools)
install_github("agentlans/PeakSorter")
```

# Usage
Let `peak.data.frame` be a data frame where:
- each row represents a LC-MS peak
- and columns are Rt, Mz, Intensity, Name
    - representing retention time (min.), m/z, intensity, and peak ID, respectively

Then to sort the peaks using the binned method, run the following code:
```r
sorted.list <- bin.prioritize(peak.data.frame)
```

To view the top 10 peaks, run:
```r
top.peaks(sorted.list, 10)
```

# License
The GNU General Public License v3.0

Software maintainer: Alan Tseng
