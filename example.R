library(data.table)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(DNABarcodes)

source("lib/lib_simulation.R")

simulate_main(
    barcode_type         = "random",
    barcode_library_file = NULL,
    clone_size_dist      = "uniform",
    clone_n              = 20,
    clone_size_dist_par  = list(size_max      = 1000, size_min = 1),
    cycle                = 30,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_cell    = 50,
    output_prefix        = "./tmp/simu_seq",
    ngs_profile          = "MSv1",
    reads_length         = 35,
    is_replicate         = F,
    top_seq              = "AAAAAAAAAAGGGGG",
    bottom_seq           = "TTTTTTTTTT",
    sequence_trunk       = 10
)
