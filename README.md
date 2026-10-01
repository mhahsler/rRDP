
<img src="https://bioconductor.org/images/logo/svg/Logo.svg" align="right" />

# R package rRDP - Interface to the RDP Classifier

[![r-universe
status](https://mhahsler.r-universe.dev/badges/rRDP)](https://mhahsler.r-universe.dev/rRDP)
[![Package on
Bioc](https://img.shields.io/badge/Bioconductor-blue)](https://bioconductor.org/packages/rRDP)

This package installs and interfaces the naive Bayesian classifier for
16S rRNA sequences developed by the Ribosomal Database Project (RDP).
With this package the classifier trained with the standard training set
can be used or a custom classifier can be trained.

Other R interfaces for bioinformatics are also available:

- [rBLAST](https://bioconductor.org/packages/rBLAST): Interface to the
  Basic Local Alignment Search Tool (BLAST) to search genetic sequence
  data bases
- [rMSA](https://mhahsler.r-universe.dev/ui#package:rMSA): Interface for
  Popular Multiple Sequence Alignment Tools including ClustalW, MAFFT,
  MUSCLE, and Kalign

## Installation

``` r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("rRDP")
```

## Usage

``` r
library(rRDP)
```

Load a few sequences and classify them with RDP. On first use, `rdp()` downloads
the default 16S classifier to the user data directory; the download may take a while.

``` r
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta", package = "rRDP"))
seq
```

    ## RNAStringSet object of length 5:
    ##     width seq                                               names               
    ## [1]  1481 AGAGUUUGAUCCUGGCUCAGAAC...GGUGAAGUCGUAACAAGGUAACC 1675 AB015560.1 d...
    ## [2]  1404 GCUGGCGGCAGGCCUAACACAUG...CACGGUAAGGUCAGCGACUGGGG 4399 D14432.1 Rho...
    ## [3]  1426 GGAAUGCUNAACACAUGCAAGUC...AACAAGGUAGCCGUAGGGGAACC 4403 X72908.1 Ros...
    ## [4]  1362 GCUGGCGGAAUGCUUAACACAUG...UACCUUAGGUGUCUAGGCUAACC 4404 AF173825.1 A...
    ## [5]  1458 AGAGUUUGAUUAUGGCUCAGAGC...UGAAGUCGUAACAAGGUAACCGU 4411 Y07647.2 Dre...

``` r
pred <- predict(rdp(), seq)
pred
```

    ##        domain         phylum               class            order
    ## 1675 Bacteria   Nitrospinota         Nitrospinia    Nitrospinales
    ## 4399 Bacteria Pseudomonadota Alphaproteobacteria Rhodospirillales
    ## 4403 Bacteria Pseudomonadota Alphaproteobacteria Rhodospirillales
    ## 4404 Bacteria Pseudomonadota Alphaproteobacteria Rhodospirillales
    ## 4411 Bacteria Pseudomonadota Alphaproteobacteria Rhodospirillales
    ##                 family          genus
    ## 1675    Nitrospinaceae     Nitrospina
    ## 4399 Rhodovibrionaceae    Rhodovibrio
    ## 4403  Acetobacteraceae    Roseococcus
    ## 4404  Acetobacteraceae Sediminicoccus
    ## 4411  Acetobacteraceae           <NA>

``` r
attr(pred, "confidence")
```

    ##      domain phylum class order family genus
    ## 1675      1      1     1     1      1  1.00
    ## 4399      1      1     1     1      1  1.00
    ## 4403      1      1     1     1      1  1.00
    ## 4404      1      1     1     1      1  1.00
    ## 4411      1      1     1     1      1  0.29

## Citation Request

To cite package ‘rRDP’ in publications use:

> Hahsler M, Nagar A (2020). “rRDP: Interface to the RDP Classifier.”
> Bioconductor version: Release (3.19). <doi:10.18129/B9.bioc.rRDP>
> <https://doi.org/10.18129/B9.bioc.rRDP>, R package version 1.23.3.

    @Misc{,
      title = {{rRDP:} Interface to the {RDP} Classifier},
      author = {Michael Hahsler and Annurag Nagar},
      year = {2020},
      doi = {10.18129/B9.bioc.rRDP},
      note = {R package version 1.23.3},
      howpublished = {Bioconductor version: Release (3.19)},
    }

## Acknowledgments

This work was partially supported by grant no. R21HG005912 from the
[National Human Genome Research Institute](https://www.genome.gov/).
