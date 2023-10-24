#!/usr/bin/env Rscript

#' # Example 1: Simulate a Barcode sequencing experiment

library(CellBarcodeSim)

#' Load barcode library

barcode_library_file = system.file("data", "random_barcodes.tsv", package = "CellBarcodeSim")

#' Simulate a Barcode sequencing experiment

## Simulate a Barcode sequencing experiment
simulate_main(
    barcode_library_file = barcode_library_file,  ## Define the barcode library
    clone_n              = 20,  ## Define the number of clones
    clone_size_dist      = "uniform", ## Define the clone size distribution
    clone_size_dist_par  = list(size_max = 1000, size_min = 1),  ## Define the parameters of the clone size distribution
    cycle                = 30,  ## Define the number of cell cycles
    efficiency           = 0.705, ## Define the efficiency of the cell cycle
    error                = 1e-6,  ## Define the error rate of the PCR, mutation per base per cycle
    pcr_read_per_cell    = 50,  ## Define the number of PCR reads per cell (clone_n)
    output_prefix        = "./tmp/barcode_sim" ,  ## Define the output prefix
    ngs_profile          = "MSv1",  ## Define the NGS profile (refer to ART sequencing simulator manual)
    reads_length         = 35,  ## Define the length of the reads
    is_replicate         = F,  ## Define whether to divide the reads into two replicates
    top_seq              = "AAAAAAAAAAGGGGG",  ## Define the fixed sequence at the 5' end of the reads to be added
    bottom_seq           = "TTTTTTTTTT",  ## Define the fixed sequence at the 3' end of the reads to be added
    sequence_trunk       = 10, ## Define the length of the fixed sequence to be added
    art_bin              = NULL ## Use the default ART binary file
)

