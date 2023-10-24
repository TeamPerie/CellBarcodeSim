#!/usr/bin/env Rscript

#' # This script will show you some examples of how to simulate a barcode library.

#' ## Load the required packages
library(CellBarcodeSim)

temp_dir = tempdir()

#' ## Example A. Simulate a random barcode library with 10^5 barcode sequences of length 15 by uniform distribution
simu_barcode_random_uniform(
    length = 15,    ## the length of the barcode sequences
    n = 1e5,        ## the number of barcode sequences
    output = paste0(temp_dir, "/random_barcodes_uniform.tsv") ## the output file
)

#' ## Example B. Simulate a random barcode library with 10^5 barcode sequences of length 15 by normal distribution
## The mean of the normal distribution is 0.5 and the standard deviation is 0.1
simu_barcode_random_norm(
    length = 15,    ## the length of the barcode sequences
    n = 1e5,        ## the number of barcode sequences
    mean = 0.5,     ## the mean of the normal distribution
    sd = 0.1,       ## the standard deviation of the normal distribution
    output = paste0(temp_dir, "/random_barcodes_norm.tsv") ## the output file
)

#' ## Example C. Simulate a random barcode library with 10^5 barcode sequences of length 15 by lognormal distribution
## The mean of the lognormal distribution is 0.5 and the standard deviation is 0.1
simu_barcode_random_lnorm(
    length = 15,    ## the length of the barcode sequences
    n = 1e5,        ## the number of barcode sequences
    log_mean= 0.5,  ## the mean of the lognormal distribution
    log_sd= 0.1,    ## the standard deviation of the lognormal distribution
    output = paste0(temp_dir, "/random_barcodes_lnorm.tsv") ## the output file
)

#' ## Example D. Simulate a random barcode library with 10^5 barcode sequences of length 15 by exponential distribution
## The rate of the exponential distribution is 0.1
simu_barcode_random_exp(
    length = 15,    ## the length of the barcode sequences
    n = 1e5,        ## the number of barcode sequences
    rate = 0.1,     ## the rate of the exponential distribution
    output = paste0(temp_dir, "/random_barcodes_exp.tsv") ## the output file
)

#' ## Example E. Simulate a hamming barcode library of 10bp barcodes with minimum hamming distance of 2.
simu_barcode_hamming_uniform(
    length = 10,    ## the length of the barcode sequences
    dist = 2,        ## the minimum hamming distance
    output = paste0(temp_dir, "/hamming_barcodes_uniform.tsv") ## the output file
)
