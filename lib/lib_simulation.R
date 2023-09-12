## Creator: Wenjie SUN
## Email: sunwjie at gmail.com
## Date: 2023-06-19
## Purpose: simulation functions for single cell barcode sequencing data

######################
#  Clone simulation  #
######################

#' log normal: 
#'  parameters: 
#'      clone size variation
#'      clone number
#'      average clone size
simu_clone_size_lognormal = function(
    n,                                 # clone number
    size_mean = 1.2,                   # meanlog
    size_variant = 2                   # meansd
    ) {
    size_v = rlnorm(n, meanlog = size_mean, sdlog = size_variant) %>% ceiling()
    size_v[size_v <= 1] = 1
    size_v
}

#' uniform:
#'  parameters:
#'      max clone size
#'      min clone size
#'      clone number
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
simu_barcode_hamming_uniform = function(
    length = 10,
    dist = 3,
    output = "hamming_barcodes.tsv"
    ) {
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    ## save as two column: seq, freq
    x = table(b_v)
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}

#' Simulate hamming distance barcode library with normal distribution
simu_barcode_hamming_norm = function(
    length = 10,
    dist = 3,
    n = 1e6,
    mean = 1,
    sd = 1,
    output = "hamming_barcodes_norm.tsv"
    ) {
    b_freq = rnorm(n, mean = mean, sd = sd) %>% ceiling()
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}

#' Simulate hamming distance barcode library with log normal distribution
simu_barcode_hamming_lnorm = function(
    length = 10,
    dist = 3,
    n = 1e6,
    log_mean = 1,
    log_sd = 1,
    output = "hamming_barcodes_lnorm.tsv"
    ) {
    b_freq = rlnorm(n, meanlog = log_mean, sdlog = log_sd) %>% ceiling()
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}

#' Simulate hamming distance barcode library with exponential distribution
simu_barcode_hamming_exp = function(
    length = 10,
    dist = 3,
    n = 1e6,
    rate = 1,
    output = "hamming_barcodes_exp.tsv"
    ) {
    b_freq = rexp(n, rate) %>% ceiling()
    b_v = DNABarcodes::create.dnabarcodes(n = length, dist = 3)
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}


#' Simulate random barcode library with uniform distribution
simu_barcode_random_uniform = function(
    length = 14,
    n = 1e6,
    output = "random_barcodes_uniform.tsv"
    ) {
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(b_v)
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}

#' Simulate random barcode library with normal distribution
simu_barcode_random_norm = function(
    length = 14,
    n = 1e6,
    mean = 1,
    sd = 1,
    output = "random_barcodes_norm.tsv"
    ) {
    b_freq = rnorm(n, mean = mean, sd = sd) %>% ceiling()
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}

#' Simulate random barcode library with log normal distribution
simu_barcode_random_lnorm = function(
    length = 14,
    n = 1e6,
    log_mean = 1,
    log_sd = 1,
    output = "random_barcodes_lnorm.tsv"
    ) {
    b_freq = rlnorm(n, meanlog = log_mean, sdlog = log_sd) %>% ceiling()
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}


#' Simulate random barcode library with exponential distribution
simu_barcode_random_exp = function(
    length = 14,
    n = 1e6,
    rate = 1,
    output = "random_barcodes_exp.tsv"
    ) {
    b_freq = rexp(n, rate) %>% ceiling()
    b_v = stringi::stri_rand_strings(n = n, length = length, pattern = '[ATCG]')
    x = table(rep(b_v, b_freq))
    d = data.frame(seq = names(x), freq = as.integer(x))
    readr::write_tsv(d, file)
}



#' Function loading barcode library file
simu_barcode = function(
    n,
    file = NULL
    ) {
    if (is.null(file)) {
        stop("no barcode library file!")
    }

    d_barcode_lib = fread(file)

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
current_dir <- dirname(parent.frame(2)$ofile)
Rcpp::sourceCpp(paste0(current_dir, "/lib_pcr_simulation.cpp"))

#' Input
#'  temp: 
#'      A list with two items. `seq` is character vector keeps sequences;
#'      `freq` is a integer vector keeps the frequency of each sequences.
#'  cycle: PCR cycle
#'  efficiency: pcr efficiency.
#'  error: pcr error per base per cycle.
#'  reads: reads number sampled from PCR results for sequencing.
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

#' Need to download the ART to the `lib/ngs_simu` fold.

#' MiSeq sequencing model

#' parameters
#'  input: fasta file
#'  profile: it can be MSv1, HS20 ...
#'  reads_length: default is 110
#'  output: the output of the simulated sequencing reads
simu_sequence_run_command = function(input, profile, reads_length = 110, output, art_bin = "./lib/art_bin_MountRainier/art_illumina") {
    cmd = str_glue("{art_bin} -ss {profile} -i {input} -amp -o {output} -l {reads_length} -f 1")
    system(cmd)
}

##################################
#  Main function for simulation  #
##################################

#' Run the simulation
simulate_main = function(
    barcode_library_file = NULL,
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
    is_replicate         = F,
    top_seq              = "",
    bottom_seq           = "",
    sequence_trunk       = 10,
    art_bin              = "./lib/art_bin_MountRainier/art_illumina"
    ) {

    dir.create("./tmp/", showWarnings=F)

    ## simulate barcode
    d_barcode_label = simu_barcode(
        clone_n,
        file = barcode_library_file
    )
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
    } else {
        stop(str_glue("The clone size distribution {clone_size_dist} is not exist."))
    }

    ## barcoded labeled cells
    d_cell = data.frame(seq = d_barcode_label, freq = clone_size_v)
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
        art_bin = art_bin
        )

    if (is_replicate) {
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

simulate_main_umi = function(
    barcode_library_file = NULL,
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
    art_bin              = "./lib/art_bin_MountRainier/art_illumina"
    ) {

    dir.create("./tmp/", showWarnings=F)
    ## simulate barcode
    d_barcode_label = simu_barcode(
        clone_n,
        file = barcode_library_file
    )
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
    } else {
        stop(str_glue("The clone size distribution {clone_size_dist} is not exist."))
    }

    ## barcoded labeled cells
    d_cell = data.frame(seq = d_barcode_label, freq = clone_size_v)
    cell_barcode = str_glue("{output_prefix}_ref.tsv")
    readr::write_tsv(d_cell, cell_barcode)

    ## PCR
    
    if (preamp_n == 0) {

        umi_seq = stringi::stri_rand_strings(n = sum(d_cell$freq), length = umi_length, pattern = '[ATCG]')

        ## UMI tagging
        x = paste0(umi_seq, rep(d_cell$seq, d_cell$freq))
        x = sample(x, round(length(x) * umi_tagging_efficiency), replace = F)
        x = table(x)
        d_cell = data.frame(seq = names(x), freq = as.integer(x))

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
        d_cell = data.frame(seq = names(x), freq = as.integer(x))

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
        art_bin = art_bin
        )

    
    list(
        library_fasta = library_fasta,
        sequencing_result = output_prefix
    )
}
