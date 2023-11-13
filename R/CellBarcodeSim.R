#' CellBarcodeSim 
#'
#' This package aim to simulate DNA barcode sequencing results,
#' including both UMI (unique molecular identifier) sequencing and non-UMI
#' sequencing techniques.
#'
#'
#' @name CellBarcodeSim
#' @docType package
#' @importFrom magrittr %>% %<>% extract extract2
#' @importFrom data.table data.table rbindlist .N := fread
#' @importFrom plyr . count
#' @importFrom stringr str_glue str_match str_split_fixed str_c
#' @import readr
#' @import methods
#' @import Rcpp
#' @useDynLib CellBarcodeSim
NULL
