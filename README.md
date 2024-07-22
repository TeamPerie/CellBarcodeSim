# DNA barcode sequencing simulation kit

The scripts in this repo were used to simulate DNA barcode sequencing results, including both UMI (unique molecular identifier) sequencing and non-UMI sequencing techniques.

## Configuration

### Install ART (optional)

The ART sequencing simulator is included in this repo, but you can also install it by yourself.
The building ART simulator can work on 64-bit Linux, MacOS and Windows.

In case to download the ART NGS simulator, please check out following link:

<https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>

To use the customed ART simulator, please make the `art_bin` variable point to the `art_illumina` file.
Otherwise, please keep the `art_bin` variable as default, which is `NULL`, to use the built-in ART simulator.

### Install R dependencies

The following R packages are required for running the simulation:

Script for installing the simulator:

```R
## install remotes package if you do not have it
install.packages("remotes")

## install CellBarcodeSim from github
remotes::install_github("TeamPerie/CellBarcodeSim")
```

## Run the simulation

### Example

The `example1_simulate_sequencing.R` script in the viggnettes fold will give you a demo of barcode sequencing simulation, you can find in the root dir.
Normally, it takes less than 1 min to run it.

Though we provided some barcode library in `viggnettes/` folder, if you need to generate more barcode libraries with customed parameters you can explore the example in `example2_simulate_barcode_library.R`.

### The built-in barcode library

The barcode library file is needed for the sequencing simulation function described below.

The file should be a tab-separated file (TSV) and consist of two columns, 

1. `seq` column: barcode sequence.
2. `freq` column: barcode frequency. It is used as barcode probability while sampling the barcode library to label the cells.

Three example barcode libraries are in the `inst/data` fold. They are:

1. Random barcode: `random_barcodes.`tsv` contains a 14bp random barcode library consisting of 1e6 simulated sequences. Exists some duplicated sequences.
2. Hamming distance barcode: `hamming_barcodes.tsv` contains 9155 14bp sequences that have been identified with a minimum Hamming distance of 3 between each sequence.
3. Vdj barcode: `vdj_barcodes.tsv` contains $10^7$ simulated VDJ barcodes that contain around $1.4\times10^5$ unique barcodes.

Once the barcode library is loaded, you can get the barcode library by using:

```r
barcode_library = system.file("data", "random_barcodes.tsv", package = "CellBarcodeSim")
```

The simulation methods were described in the paper (TBD: add paper link).

### The barcode library simulation function

There are also 8 functions to simulate barcode libraries.

They simulated 2 types of barcodes library X 4 barcode library distribution due to potential barcode amplification bias.

The 2 types of barcodes library are:

- Random barcode: The barcode sequences are random.
- Hamming distance barcode: The barcode sequences have a minimum Hamming distance between each other.

The 4 barcode library distribution are:

- Uniform distribution.
- Normal distribution.
- Lognormal distribution.
- Exponential distribution.

#### Hamming distance barcode

The four functions used to simulate the hamming distance barcode library are:

- `simu_barcode_hamming_uniform()`: Simulate the hamming distance barcode library with uniform distribution.
- `simu_barcode_hamming_norm()`: Simulate the hamming distance barcode library with normal distribution.
- `simu_barcode_hamming_lnorm()`: Simulate the hamming distance barcode library with lognormal distribution.
- `simu_barcode_hamming_exp()`: Simulate the hamming distance barcode library with exponential distribution.

Those wrap the barcode simulator in the [DNAbarcodes](https://bioconductor.org/packages/release/bioc/html/DNABarcodes.html) package to simulate the barcode library given hamming distance.
It accepts parameters:

1. `length`:  barcode length.
2. `dist`: hamming distance.
3. `output`: output file path. 
4. `mean`: mean of the distribution, only for `simu_barcode_hamming_norm()`.
5. `sd`: standard deviation of the distribution, only for `simu_barcode_hamming_norm()`.
6. `log_mean`: log mean of the distribution, only for `simu_barcode_hamming_lnorm()`.
7. `log_sd`: log standard deviation of the distribution, only for `simu_barcode_hamming_lnorm()`.
8. `rate`: rate of the distribution, only for `simu_barcode_hamming_exp()`.
9. `top_fix`: 5 end fix flanking sequence.
10. `bottom_fix`: 3 end fix flanking sequence.

The output file dir is given, the output file is a two columns TSV file which is described in the previous section.

Otherwise, the function will return a `data.table` with `seq` and `freq` column.

For more complex cases please simulate directly with the [DNABarcodes](https://bioconductor.org/packages/release/bioc/html/DNABarcodes.html) package.

#### Random barcode

The four functions used to simulate the random barcode library are:

- `simu_barcode_random_uniform()`: Simulate the random barcode library with uniform distribution.
- `simu_barcode_random_norm()`: Simulate the random barcode library with normal distribution.
- `simu_barcode_random_lnorm()`: Simulate the random barcode library with lognormal distribution.
- `simu_barcode_random_exp()`: Simulate the random barcode library with exponential distribution.

They accept parameters:

1. `length`: barcode length.
2. `n`: number of sequences to be simulated. The unique barcodes can be less than the total sequence number due to the potential duplications.
3. `output`: output file path.
4. `mean`: mean of the distribution, only for `simu_barcode_random_norm()`.
5. `sd`: standard deviation of the distribution, only for `simu_barcode_random_norm()`.
6. `log_mean`: log mean of the distribution, only for `simu_barcode_random_lnorm()`.
7. `log_sd`: log standard deviation of the distribution, only for `simu_barcode_random_lnorm()`.
8. `rate`: rate of the distribution, only for `simu_barcode_random_exp()`.

### Non UMI simulation

We can run the `simulate_main()` function to do the simulation without UMI.
You can find an example in `example.R` file.

Following are the parameters and the default:

```r
simulate_main(
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
    reads_length         = 10,
    is_replicate         = F,s
    top_seq              = "",
    bottom_seq           = "",
    sequence_trunk       = 10,
    art_bin              = NULL
)
```

The parameters:

- `barcode_library_file`: The location of the barcode frequency list.
This parameter is necessary and has no default value. We sample the list to have the barcode for labeling.
Please refer to the sections above to know how to prepare it.
- `barcode_library`: If the `barcode_library_file` is empty, a `data.frame` of barcode library is used here. The `data.frame` should contain two columns `seq` and `freq`.
- `clone_size_dist`: clone size distribution, it can be `uniform` or `lognormal`.
The `clone_size_dist_par` should match the chosen model (see more detail below).
- `clone_n`: Number of labeled cells (progenitors) at initiation, default 300.
- `clone_size_dist_par`: a list.
If `uniform` clone size distribution is chosen, the list should contain two items `size_max` and `size_min` for the range of the uniform distribution. 
If `lognormal` is chosen, the list should contain `size_mean` which is log mean, and `size_variant` which is log sd.
If `powerlaw` clone size distribution is chosen, the list should contain three items `constant`, `scale` and `alpha`. They are corresponding to the parameters in following equation:
$$f(scale \cdot x) = constant (scale \cdot x)^{-\alpha}$$
Where $x \in [0,1]$.
- `cycle`: PCR cycle number.
The run time and output file size increase will increase exponentially with the `cycle` value, because of the exponential nature of PCR amplification.
- `efficiency`: PCR efficiency with a default value of 0.705 for each cycle.
- `error`: PCR error per base with a default value of $10^{-6}$.
- `pcr_read_per_cell`: Expacted sequencing reads per each progeny cell.
- `output_prefix`: Output prefix.
- `ngs_profile`: NGS profile, this option will be transmitted to ART sequencing simulator.
The building profile is, "MSv1", "HS20" etc.
For more detail, please check the ART simulator.
- `reads_length`: The NGS sequencing reads length.
**It should be shorter than the PCR results.**
You can use `top_seq` or `bottom_seq` (described below) to add the constant sequences to increase the PCR sequence length.
- `is_replicate`: Generate technical replicates or not. The technical replicates will generate by dividing the progeny cells into two samples and do the PCR for each. The sequencing data will be labeled by pending `_2` to label the replicate.
- `top_seq`: Add constant sequences to the left (5') of the PCR result.
- `bottom_seq`: Add constant sequences to the right (3') of the PCR result.
- `sequence_trunk`: Substring the barcode sequence by choosing the first n base, default 10.
If the value is bigger than the barcode length, the barcode sequence will not be truncated.
- `art_bin`: The executable ART simulator location. The default is `NULL`, which makes use of the building ART simulator.

### UMI sequencing simulation

We can simulate the UMI sequencing result with the function `simulate_main_umi()`.

```R
simulate_main_umi(
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
    art_bin              = NULL
)
```

For most of the parameters, you can find explanations in the documents of `simu_main`.
The `simulate_main_umi()` specific parameters:

- `preamp_n`: PCR cycle before catenating UMI, default is 0.
- `umi_length`: The UMI base pair length, default is 8.
- `umi_tagging_efficiency`: The UMI tagging efficiency, default is 0.25.
It means 25% of sequences in the pre-PCR pool will be tagged by UMI.

## About the output

The output location and file name are defined by the `output_prefix` option.
For example, by running the `example.R`, the option is `./tmp/simu_seq`. The output files are in `./tmp/` fold:

- `simu_seq_ref.tsv`: The induced barcodes.
- `simu_seq_library.fasta`: The sequencing library, which will only carry the PCR error.
- `simu_seq.aln`: The intermediate files from ART.
- `simu_seq.fq` : Fianl barcode sequencing fastq file.
