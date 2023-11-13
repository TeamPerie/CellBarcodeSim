######################
#  Clone simulation  #
######################

#' Simulate clone size with log normal distribution
#' 
#' @param n clone number
#' @param size_mean meanlog of log normal distribution with default 1.2
#' @param size_variant meansd of log normal distribution with default 2
#' @return a vector of clone size
#' @export
simu_clone_size_lognormal = function(
    n,                                 # clone number
    size_mean = 1.2,                   # meanlog
    size_variant = 2                   # meansd
    ) {
    size_v = rlnorm(n, meanlog = size_mean, sdlog = size_variant) %>% ceiling()
    size_v[size_v <= 1] = 1
    size_v
}

#' Simulate clone size with power law distribution
#'
#' @param n clone number
#' @param constant constant of power law distribution with default 10
#' @param scale scale of power law distribution with default 1
#' @param alpha alpha of power law distribution with default 2
#' @return a vector of clone size
#' @export
simu_clone_size_powerlaw = function(
    n,
    constant = 10,
    scale = 1,
    alpha = 2
    ) {
    size_v = constant * (scale * runif(n))^-alpha
    size_v[size_v <= 1] = 1
    size_v
}

#' Simulate clone size with exponential distribution
#'
#' @param n clone number
#' @param size_max max clone size with default 1000
#' @param size_min min clone size with default 1
#' @return a vector of clone size
#' @export
simu_clone_size_uniform = function(
    n,                                 # clone number
    size_max = 1000,
    size_min = 1
    ) {
    seq.int(size_min, size_max, length.out = n)
}

########################
#  Barcode simulation  #
########################

#' Simulate hamming distance barcode library with uniform distribution
#' 
#' @param length barcode length with default 10
#' @param dist hamming distance with default 3
#' @param output output file name with default "hamming_barcodes.tsv"
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @return a barcode library file or a data frame
#' @export
simu_barcode_hamming_uniform = function(
    length = 10,
    dist = 3,
    top_fix = "",
    bottom_fix = "",
    output = "hamming_barcodes.tsv"
    ) {
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    ## save as two column: seq, freq
    x = table(b_v)
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}

#' Simulate hamming distance barcode library with normal distribution
#'
#' @param length barcode length with default 10
#' @param dist hamming distance with default 3
#' @param n number of barcodes with default 1e6
#' @param mean mean of normal distribution with default 1
#' @param sd sd of normal distribution with default 1
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "hamming_barcodes_norm.tsv"
#' @return a barcode library file or a data frame
#' @export
simu_barcode_hamming_norm = function(
    length = 10,
    dist = 3,
    n = 1e6,
    mean = 1,
    sd = 1,
    top_fix = "",
    bottom_fix = "",
    output = "hamming_barcodes_norm.tsv"
    ) {
    b_freq = rnorm(n, mean = mean, sd = sd) %>% ceiling()
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}

#' Simulate hamming distance barcode library with log normal distribution
#'
#' @param length barcode length with default 10
#' @param dist hamming distance with default 3
#' @param n number of barcodes with default 1e6
#' @param log_mean meanlog of log normal distribution with default 1
#' @param log_sd sdlog of log normal distribution with default 1
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "hamming_barcodes_lnorm.tsv"
#' @return a barcode library file or a data frame
#' @export
simu_barcode_hamming_lnorm = function(
    length = 10,
    dist = 3,
    n = 1e6,
    log_mean = 1,
    log_sd = 1,
    top_fix = "",
    bottom_fix = "",
    output = "hamming_barcodes_lnorm.tsv"
    ) {
    b_freq = rlnorm(n, meanlog = log_mean, sdlog = log_sd) %>% ceiling()
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}

#' Simulate hamming distance barcode library with exponential distribution
#' 
#' @param length barcode length with default 10
#' @param dist hamming distance with default 3
#' @param n number of barcodes with default 1e6
#' @param rate rate of exponential distribution with default 1
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "hamming_barcodes_exp.tsv"
#' @return a barcode library file or a data frame
simu_barcode_hamming_exp = function(
    length = 10,
    dist = 3,
    n = 1e6,
    rate = 1,
    top_fix = "",
    bottom_fix = "",
    output = "hamming_barcodes_exp.tsv"
    ) {
    b_freq = rexp(n, rate) %>% ceiling()
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}


#' Simulate random barcode library with uniform distribution
#'
#' @param length barcode length with default 14
#' @param n number of barcodes with default 1e6
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "random_barcodes_uniform.tsv"
#' @return a barcode library file or a data frame
#' @export
simu_barcode_random_uniform = function(
    length = 14,
    n = 1e6,
    top_fix = "",
    bottom_fix = "",
    output = "random_barcodes_uniform.tsv"
    ) {
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(b_v)
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}

#' Simulate random barcode library with normal distribution
#'
#' @param length barcode length with default 14
#' @param n number of barcodes with default 1e6
#' @param mean mean of normal distribution with default 1
#' @param sd sd of normal distribution with default 1
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "random_barcodes_norm.tsv"
#' @return a barcode library file or a data frame
#' @export
simu_barcode_random_norm = function(
    length = 14,
    n = 1e6,
    mean = 1,
    sd = 1,
    top_fix = "",
    bottom_fix = "",
    output = "random_barcodes_norm.tsv"
    ) {
    b_freq = rnorm(n, mean = mean, sd = sd) %>% ceiling()
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}

#' Simulate random barcode library with log normal distribution
#'
#' @param length barcode length with default 14
#' @param n number of barcodes with default 1e6
#' @param log_mean meanlog of log normal distribution with default 1
#' @param log_sd sdlog of log normal distribution with default 1
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "random_barcodes_lnorm.tsv"
#' @return a barcode library file or a data frame
#' @export
simu_barcode_random_lnorm = function(
    length = 14,
    n = 1e6,
    log_mean = 1,
    log_sd = 1,
    top_fix = "",
    bottom_fix = "",
    output = "random_barcodes_lnorm.tsv"
    ) {
    b_freq = rlnorm(n, meanlog = log_mean, sdlog = log_sd) %>% ceiling()
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}


#' Simulate random barcode library with exponential distribution
#'
#' @param length barcode length with default 14
#' @param n number of barcodes with default 1e6
#' @param rate rate of exponential distribution with default 1
#' @param top_fix 5 end fix region sequence of the barcode with default "".
#' @param bottom_fix 3 end fix region sequence of the barcode with default "".
#' @param output output file name with default "random_barcodes_exp.tsv"
#' @return a barcode library file or a data frame
#' @export
simu_barcode_random_exp = function(
    length = 14,
    n = 1e6,
    rate = 1,
    top_fix = "",
    bottom_fix = "",
    output = "random_barcodes_exp.tsv"
    ) {
    b_freq = rexp(n, rate) %>% ceiling()
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    d$x = paste0(top_fix, d$seq, bottom_fix)
    if (is.null(output)) {
        return(d)
    } else {
        readr::write_tsv(d, output)
    }
}


#' Function sample barcode library
#' 
#' @param n number of barcodes to sample
#' @param d_barcode_lib a data.frame with two columns: `seq` and `freq`.
#' @return a data frame
#' @export
sample_barcode = function(
    n,
    d_barcode_lib = NULL
    ) {

    res = sample(
        d_barcode_lib$seq, 
        n, 
        replace=TRUE, 
        prob = d_barcode_lib$freq)
    res
}


####################
#  PCR simulation  #
####################
# Get the directory of the currently executing script
# current_dir <- dirname(parent.frame(2)$ofile)
# Rcpp::sourceCpp(paste0(current_dir, "/lib_pcr_simulation.cpp"))

#' PCR amplification
#'
#' @param temp a list with two items: `seq` is character vector keeps sequences; `freq` is a integer vector keeps the frequency of each sequences.
#' @param cycle PCR cycle
#' @param efficiency pcr efficiency.
#' @param error pcr error per base per cycle.
#' @param reads reads number sampled from PCR results for sequencing.
#' @return a list with two items: `seq` is character vector keeps sequences; `freq` is a integer vector keeps the frequency of each sequences.
#' @export
do_pcr = function(temp, cycle = 30, efficiency = 0.705, error = 1e-6, reads = 5e5) {
    ## do the PCR amplification
    pcr_res = pcr_amplify(temp, cycle, efficiency, error)
    if (reads > 0) {
        sample(pcr_res$seq, reads, replace = T, prob = pcr_res$freq / sum(pcr_res$freq))
    } else {
        pcr_res
    }
}


###########################
#  Sequencing simulation  #
###########################

# Need to download the ART to the `lib/ngs_simu` fold.
# MiSeq sequencing model

#' Run the ART sequencing simulator
#'
#' @param input input fasta file
#' @param profile sequencing profile
#' @param output output prefix of the simulated sequencing result
#' @param reads_length reads length of the simulated sequencing result with default 110
#' @param qc_shift quality score shift, default is 0
#' @param art_bin path to the ART bin directory, default is NULL, which will use the package buildin ART bin directory.
simu_sequence_run_command = function(input, profile, output, reads_length = 110, qc_shift = 0, art_bin = NULL) {

    ## configure the art binary
    if (is.null(art_bin)) {
        os <- tolower(Sys.info()["sysname"])

        if (grepl("darwin", os)) {
            art_dir <- system.file("data", "art_bin_MountRainier_Mac", package = "CellBarcodeSim")
        } else if (grepl("mingw|win", os)) {
            art_dir <- system.file("data", "art_bin_MountRainier_Win", package = "CellBarcodeSim")
        } else if (grepl("linux", os)) {
            art_dir <- system.file("data", "art_bin_MountRainier_Linux", package = "CellBarcodeSim")
        } else {
            stop("Unknown OS, please install ART sequencing simulator manually.\n
                https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm\n
                And set the path to the ART bin directory.")
        }
        art_bin = str_glue("{art_dir}/art_illumina")
    } else {
        if (!dir.exists(art_bin)) {
            stop("The ART bin directory is not exist!")
        }
    }

    buildin_profile = c("MSv1", "HS20")

    ## configure the profile
    if (file.exists(profile)) {
        cmd = str_glue("{art_bin} -qs {qc_shift} -1 {profile} -i {input} -amp -o {output} -l {reads_length} -f 1")
    } else if (profile %in% buildin_profile) {
        cmd = str_glue("{art_bin} -ss {profile}  -qs {qc_shift} -i {input} -amp -o {output} -l {reads_length} -f 1")
    } else {
        stop(str_glue("The profile/file {profile} is not exist."))
    }

    system(cmd)
}

##################################
#  Main function for simulation  #
##################################

#' Run the non-UMI barcode simulation
#'
#' @param barcode_library_file barcode library file, if not provided, the barcode_library should be provided.
#' @param barcode_library barcode library data frame with two columns: `seq` and `freq`.
#' @param clone_size_dist clone size distribution, should be one of `uniform` and `lognormal`.
#' @param clone_n number of clones.
#' @param clone_size_dist_par parameters for clone size distribution, check the function \code{\link{simu_clone_size_uniform}}, \code{\link{simu_clone_size_lognormal}}, and \code{\link{simu_clone_size_powerlaw}}.
#' @param cycle PCR cycle count. 
#' @param efficiency pcr efficiency, default is 0.705.
#' @param error pcr error per base per cycle with default 1e-6.
#' @param pcr_read_per_cell reads number sampled from PCR results for sequencing with default 50.
#' @param output_prefix output prefix of the simulated sequencing result with default "./tmp/simu_seq".
#' @param ngs_profile sequencing profile with default "MSv1".
#' @param reads_length reads length of the simulated sequencing result with default 100.
#' @param top_seq 5 end fix region sequence added after PCR default "".
#' @param bottom_seq 3 end fix region sequence added after PCR default "".
#' @param sequence_trunk the length of the barcode sequence to be used for the simulation with default 10.
#' @param qc_shift quality score shift, default is 0.
#' @param art_bin path to the ART bin directory, default is NULL, which will use the package buildin ART bin directory.
#' @return a list with two items: `library_fasta` is the fasta file of the simulated library; `sequencing_result` is the prefix of the simulated sequencing result.
#' @export
simulate_main = function(
    barcode_library_file = NULL,
    barcode_library      = NULL,
    clone_size_dist      = "uniform",
    clone_n              = 300,
    clone_size_dist_par  = list(size_max = 1000, size_min = 1),
    cycle                = 30,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_cell    = 50,
    output_prefix        = "./tmp/simu_seq",
    ngs_profile          = "MSv1",
    reads_length         = 100,
    top_seq              = "",
    bottom_seq           = "",
    sequence_trunk       = 10,
    qc_shift             = 0,
    art_bin              = NULL
    ) {

    dir.create("./tmp/", showWarnings=F)

    if (!is.null(barcode_library_file)) {
        d_barcoce_library = fread(barcode_library_file)
        d_barcode_label = sample_barcode(
            clone_n,
            d_barcode_lib = d_barcoce_library
        )
    } else if (!is.null(barcode_library)) {
        d_barcode_label = sample_barcode(
            clone_n,
            d_barcode_lib = barcode_library
        )
    } else {
        stop("No barcode library provided!")
    }

    ## set the barcode length for random barcode
    d_barcode_label = substring(d_barcode_label, 1, sequence_trunk)
    
    ## label the cells
    if (clone_size_dist == "uniform") {
        clone_size_v = simu_clone_size_uniform(
            clone_n, 
            size_max = clone_size_dist_par$size_max,
            size_min = clone_size_dist_par$size_min
            )
    } else if (clone_size_dist == "lognormal") {
        clone_size_v = simu_clone_size_lognormal(
            clone_n,
            size_mean = clone_size_dist_par$size_mean,
            size_variant = clone_size_dist_par$size_variant
            )
    } else if (clone_size_dist == "powerlaw") {
        clone_size_v = simu_clone_size_powerlaw(
            clone_n,
            constant = clone_size_dist_par$constant,
            scale = clone_size_dist_par$scale,
            alpha = clone_size_dist_par$alpha
            )
    } else {
        stop(str_glue("The clone size distribution {clone_size_dist} is not exist."))
    }

    ## barcoded labeled cells
    d_cell = data.table(seq = d_barcode_label, freq = clone_size_v)
    d_cell = d_cell[, .(freq = sum(freq)), by = seq]
    cell_barcode = str_glue("{output_prefix}_ref.tsv")
    readr::write_tsv(d_cell, cell_barcode)

    ## calc how many reads to seq
    pcr_read_n = pcr_read_per_cell * sum(d_cell$freq)

    ## PCR process
    d_library = do_pcr(
        temp = d_cell,
        cycle = cycle,
        efficiency = efficiency,
        error = error,
        reads = pcr_read_n
    )

    ## Sequence
    library_fasta = str_glue("{output_prefix}_library.fasta")
    x = paste0(">", 1:length(d_library), "\n", paste0(top_seq, d_library, bottom_seq))
    readr::write_lines(x, library_fasta)

    simu_sequence_run_command(
        input = library_fasta,
        profile = ngs_profile,
        reads_length = reads_length,
        output = output_prefix,
        qc_shift = qc_shift,
        art_bin = art_bin
        )

    if (FALSE) {
        output_prefix2 = paste0(output_prefix, "_2")
        library_fasta2 = str_glue("{output_prefix2}_library.fasta")

        d_library = do_pcr(
            temp = d_cell,
            cycle = cycle,
            efficiency = efficiency,
            error = error,
            reads = pcr_read_n
        )
        x = paste0(">", 1:length(d_library), "\n", substring(paste0(top_seq, d_library, bottom_seq)), 1, reads_length)
        readr::write_lines(x, library_fasta2)

        simu_sequence_run_command(
            input = library_fasta2,
            profile = ngs_profile,
            reads_length = reads_length,
            output = output_prefix2,
            qc_shift = qc_shift,
            art_bin = art_bin
        )

        library_fasta = c(library_fasta, library_fasta2)
        output_prefix = c(output_prefix, output_prefix2)
    }

    list(
        library_fasta = library_fasta,
        sequencing_result = output_prefix
    )
}

#######################################
#  UMI_barcode sequencing simulation  #
#######################################

#' Run the simulation for UMI_barcode sequencing
#'
#' @param barcode_library_file barcode library file, if not provided, the barcode_library should be provided.
#' @param barcode_library barcode library data frame with two columns: `seq` and `freq`.
#' @param clone_size_dist clone size distribution, should be `uniform`, `lognormal` or `powerlaw`.
#' @param clone_n number of clones.
#' @param clone_size_dist_par parameters for clone size distribution, check the function \code{\link{simu_clone_size_uniform}}, \code{\link{simu_clone_size_lognormal}}, and \code{\link{simu_clone_size_powerlaw}}.
#' @param cycle total PCR cycle count with default 40.
#' @param efficiency pcr efficiency, default is 0.705.
#' @param error pcr error per base per cycle with default 1e-6.
#' @param pcr_read_per_umi reads number sampled from PCR results for sequencing with default 50.
#' @param output_prefix output prefix of the simulated sequencing result with default "./tmp/simu_umi_seq".
#' @param ngs_profile sequencing profile with default "MSv1".
#' @param reads_length reads length of the simulated sequencing result with default 100.
#' @param top_seq 5 end fix region sequence added after PCR default "".
#' @param bottom_seq 3 end fix region sequence added after PCR default "".
#' @param sequence_trunk the max length of the barcode sequence to be used for the simulation with default 10.
#' @param preamp_n pre-amplification cycle number with default 0.
#' @param umi_length UMI length with default 8.
#' @param umi_tagging_efficiency UMI tagging efficiency with default 0.25.
#' @param qc_shift quality score shift, default is 0.
#' @param art_bin path to the ART bin directory, default is NULL, which will use the package buildin ART bin directory.
#' @return a list with two items: `library_fasta` is the fasta file of the simulated library; `sequencing_result` is the prefix of the simulated sequencing result.
#' @export
simulate_main_umi = function(
    barcode_library_file = NULL,
    barcode_library      = NULL,
    clone_size_dist      = "uniform",
    clone_n              = 300,
    clone_size_dist_par  = list(size_max = 1000, size_min = 1),
    cycle                = 40,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_umi     = 50,
    output_prefix        = "./tmp/simu_umi_seq",
    ngs_profile          = "MSv1",
    reads_length         = 100,
    top_seq              = "",
    bottom_seq           = "",
    sequence_trunk       = 10,
    preamp_n             = 0,
    umi_length           = 8,
    umi_tagging_efficiency = 0.25,
    qc_shift             = 0,
    art_bin              = NULL
    ) {

    dir.create("./tmp/", showWarnings=F)

    if (!is.null(barcode_library_file)) {
        d_barcoce_library = fread(barcode_library_file)
        d_barcode_label = sample_barcode(
            clone_n,
            d_barcode_lib = d_barcoce_library
        )
    } else if (!is.null(barcode_library)) {
        d_barcode_label = sample_barcode(
            clone_n,
            d_barcode_lib = barcode_library
        )
    } else {
        stop("No barcode library provided!")
    }

    ## set the barcode length for random barcode
    d_barcode_label = substring(d_barcode_label, 1, sequence_trunk)
    
    ## label the cells
    if (clone_size_dist == "uniform") {
        clone_size_v = simu_clone_size_uniform(
            clone_n, 
            size_max = clone_size_dist_par$size_max,
            size_min = clone_size_dist_par$size_min
            )
    } else if (clone_size_dist == "lognormal") {
        clone_size_v = simu_clone_size_lognormal(
            clone_n,
            size_mean = clone_size_dist_par$size_mean,
            size_variant = clone_size_dist_par$size_variant
            )
    } else if (clone_size_dist == "powerlaw") {
        clone_size_v = simu_clone_size_powerlaw(
            clone_n,
            constant = clone_size_dist_par$constant,
            scale = clone_size_dist_par$scale,
            alpha = clone_size_dist_par$alpha
            )
    } else {
        stop(str_glue("The clone size distribution {clone_size_dist} is not exist."))
    }

    ## barcoded labeled cells
    d_cell = data.table(seq = d_barcode_label, freq = clone_size_v)
    d_cell = d_cell[, .(freq = sum(freq)), by = seq]
    cell_barcode = str_glue("{output_prefix}_ref.tsv")
    readr::write_tsv(d_cell, cell_barcode)

    ## PCR
    
    if (preamp_n == 0) {

        umi_seq = stringi::stri_rand_strings(n = sum(d_cell$freq), length = umi_length, pattern = '[ATCG]')

        ## UMI tagging
        x = paste0(umi_seq, rep(d_cell$seq, d_cell$freq))
        x = sample(x, round(length(x) * umi_tagging_efficiency), replace = F)
        x = table(x)
        d_cell = data.table(seq = names(x), freq = as.integer(x))

        pcr_read_n = pcr_read_per_umi * sum(d_cell$freq)

        d_library = do_pcr(
            temp = d_cell,
            cycle = cycle,
            efficiency = efficiency,
            error = error,
            reads = pcr_read_n
        )
    } else {
        d_preamp = do_pcr(
            temp = d_cell,
            cycle = preamp_n,
            efficiency = efficiency,
            error = error,
            reads = -1
        )
        umi_seq = stringi::stri_rand_strings(n = sum(d_preamp$freq), length = umi_length, pattern = '[ATCG]')

        ## UMI tagging
        x = paste0(umi_seq, rep(d_preamp$seq, d_preamp$freq))
        x = sample(x, round(length(x) * umi_tagging_efficiency), replace = F)
        x = table(x)
        d_cell = data.table(seq = names(x), freq = as.integer(x))
        d_cell = d_cell[, .(freq = sum(freq)), by = seq]

        pcr_read_n = pcr_read_per_umi * sum(d_cell$freq)

        d_library= do_pcr(
            temp = d_cell,
            cycle = cycle - preamp_n,
            efficiency = efficiency,
            error = error,
            reads = pcr_read_n 
        )
    }
    
    ## Sequence
    library_fasta = str_glue("{output_prefix}_library.fasta")
    x = paste0(">", 1:length(d_library), "\n", paste0(top_seq, d_library, bottom_seq))
    readr::write_lines(x, library_fasta)

    simu_sequence_run_command(
        input = library_fasta,
        profile = ngs_profile,
        reads_length = reads_length,
        output = output_prefix,
        qc_shift = qc_shift,
        art_bin = art_bin
        )

    
    list(
        library_fasta = library_fasta,
        sequencing_result = output_prefix
    )
}
