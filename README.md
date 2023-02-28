# DNA barcode sequencing simulation kit



The scrips in this repo were used to simulate the DNA barcode sequencing results, which can be UMI and non-UMI sequencing.

## Configuration

### Install ART 

Download the ART NGS simulator from here:

https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm 

And put it to the `lib` fold. The default simulator executable file should be `./lib/art_bin_MountRainier/art_illumina`. Otherwise, the `art_bin` variable should be configured to point the `art_illumina` file.

### Install R dependencies

The scrips depends on the following packages:

-   plyr
-   magrittr
-   readr
-   stringr
-   data.table

## Run

## Non UMI simulation

We can run the `simulate_main()` function to do the simulation without UMI. The executable example is `example.R` file.

Following is the all parameters and default values:

```r
simulate_main(
    barcode_type         = "hamming",
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
    reads_length         = 10,
    is_replicate         = F,
    top_seq              = "",
    bottom_seq           = "",
    sequence_trunk       = 10,
    art_bin              = "./lib/art_bin_MountRainier/art_illumina"
)
```

The parameters:

-   `barcode_type`: the building barcode simulator is `hamming`, `random` barcode, `vdj` barcode. This option can be `custom` to allow using a given barcode frequency table as barcode library (see more detail in `barcode_library_file par`).
-   `barcode_library_file`: the location of barcode frequency list. We sample the list to have the barcode used for labeling. If the file exists, the program will use it without redo the simulation. The default value is `NULL`, in this case, the new generated barcode library file will be named `./tmp/random_barcodes.tsv`, `./tmp/hamming_barcodes.tsv`. And when `barcode_type` is `custom` , this parameter is mandatory. When a costum file is given, it should have two columns with first column named `seq` and second column named `freq`.
-   `clone_size_dist`: clone size distribution, it can be `uniform` or `lognormal`. The `clone_size_dist_par` should match the this parameter's value (see more detail in below).
-   `clone_n`: Number of progenitors, default 300.
-   `clone_size_dist_par`: a list. When `clone_size_dist` is `uniform`, the list should contain two item `size_max` and `size_min`, when `clone_size_dist` is `lognormal`, the list should contain `size_mean` - log mean and `size_variant` - log sd.
-   `cycle`: PCR cycle number. Because the nature of the exponential sequence of PCR process, the run time and output file size increase will increase exponentially with the `cycle` value.
-   `efficiency`: PCR efficiency with default value 0.705.
-   `error`: PCR error per base with default 1e-6.
-   `pcr_read_per_cell`: Reads number needed for sequence, how much per progeny cell.
-   `output_prefix`: Output prefix.
-   `ngs_profile`: NGS profile, this option will be transmit to ART sequencing simulator. The building profile is, "MSv1", "HS20" etc. For more detail, please check the ART simulator.
-   `reads_length`: The NGS sequencing reads length. It should be shorter than the PCR results. Otherwise there will be an error. In that case, you can add `top_seq` or `bottom_seq` to increase the PCR results.
-   `is_replicate`: Generate technical replicates or not. The technical replicates will generate by dividing the progeny cells into two sample and do the PCR for each of them. The two sequencing output will be labeled by pendix `_2` of the output file.
-   `top_seq`: Add 5 end constant sequences to the barcodes for sequencing simulation.
-   `bottom_seq`: Add 3 end constant sequences to the barcodes sequencing simulation.
-   `sequence_trunk`: To make the barocde sequence short by choosing the first n base, default 10.
-   `art_bin`: the executable ART simulator location. The default is "./lib/art_bin_MountRainier/art_illumina".

## UMI sequencing simulation

Similarly we can simulate the UMI sequencing result with function `simulate_main_umi()`.

Following shows all the parameters and default values.

```R
simulate_main_umi(
    barcode_type         = "hamming",
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
)
```

For the `simulate_main_umi()` specific parameters:

-   `preamp_n`: PCR cycle before conjuncting UMI, default is 0.
-   `umi_length`: The UMI base pair length, defautl is 8.
-   `umi_tagging_efficiency`: The UMI tagging efficiency, default is 0.25.

## Output

The output location and file name is defined by the `output_prefix` option. For example, by running the `example.R`, the option is `./tmp/simu_seq`. The output files are in `./tmp/` fold:

-   `simu_seq.fq` : The fastq of barcode sequencing.
-   `simu_seq.aln`: The middle file of ART.
-   `simu_seq_ref.tsv`: The induced barcodes.
-   `simu_seq_library.fasta`: The sequencing library loaded to the sequencing machine, which will carry the PCR error.

