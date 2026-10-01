# rRDP: Interface to the RDP Classifier

Abstract

This package installs and interfaces the naive Bayesian classifier for
16S rRNA sequences developed by the Ribosomal Database Project (RDP).
With this package the classifier trained with the standard training set
can be used or a custom classifier can be trained.

``` r

library("rRDP")
set.seed(1234)
```

## Installation and System Requirements

rRDP requires the Bioconductor package Biostrings and R to be configured
with Java.

``` r

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("Biostrings")
```

Install rRDP and the database used by the default RDP classifier.

``` r

install.packages(c('rRDP'), repos = c('https://mhahsler.r-universe.dev'))
```

RDP uses Java and you need a working installation of the `rJava`
package. You need to have a Java JDK installed. On Linux, you can
install Open JDK and run in your shell `R CMD javareconf` to configure R
for using Jave. On Windows, you can install the latest version of the
JDK from <https://www.oracle.com/java/technologies/downloads/> and set
the `JAVA_HOME` environment variable in R using (make sure to use the
correct location). An example would look like this:

``` r

Sys.setenv(JAVA_HOME = "C:\\Program Files\\Java\\jdk-20")
```

Note the double backslashes (i.e. escaped slashes) used in the path.

Details can be found at <https://www.rforge.net/rJava/index.html>. To
configure R for Java,

### How to cite this package

``` r
citation("rRDP")
To cite package 'rRDP' in publications use:

  Hahsler M, Nagar A (2020). "rRDP: Interface to the RDP Classifier."
  Bioconductor version: Release (3.19). doi:10.18129/B9.bioc.rRDP
  <https://doi.org/10.18129/B9.bioc.rRDP>. R package version 1.23.3.

A BibTeX entry for LaTeX users is

  @Misc{,
    title = {{rRDP:} Interface to the {RDP} Classifier},
    author = {Michael Hahsler and Annurag Nagar},
    year = {2020},
    doi = {10.18129/B9.bioc.rRDP},
    note = {R package version 1.23.3},
    howpublished = {Bioconductor version: Release (3.19)},
  }
```

## Classification with RDP

The RDP classifier was developed by the Ribosomal Database Project (Cole
et al. 2003) which provides various tools and services to the scientific
community for data related to 16S rRNA sequences. The classifier uses a
Naive Bayesian approach to quickly and accurately classify sequences.
The classifier uses 8-mer counts as features (Wang et al. 2007).

### Training a custom RDP classifier

RDP can be trained using
[`trainRDP()`](http://michael.hahsler.net/rRDP/reference/rdp.md). We use
an example of a tiny training data set that is shipped with the package.

``` r

trainingSequences <- readDNAStringSet(
    system.file("examples/trainingSequences.fasta", package = "rRDP")
)
trainingSequences
#> DNAStringSet object of length 20:
#>      width seq                                              names               
#>  [1]  1384 TAGTGGCGGACGGGTGAGTAACG...GAAGTTCGAATTTGGGTCAAGT 13652 Root;Bacter...
#>  [2]  1386 ATCTCACCTCTCAATAGCGGCGG...GCCGTCGAAGGTGGGGTTGGTG 13655 Root;Bacter...
#>  [3]  1440 ATCTCACCTCTCAATAGCGGCGG...GTGCGGCTGGATCACCTCCTTA 13661 Root;Bacter...
#>  [4]  1421 AATAGCGGCGGACGGGTGAGTAA...GCCGTATCGGAAGGTGCGGCTG 13671 Root;Bacter...
#>  [5]  1439 ATCTCACCTCTCAATANCGGCGG...CGGAAGGTGCGGCTGGATCACC 13677 Root;Bacter...
#>  ...   ... ...
#> [16]  1478 TGGCTCAGGACGAACGCTGGCGG...GTAGCCGTATCGGAAGGTGCGG 13763 Root;Bacter...
#> [17]  1507 CCTGGCTCAGGACGAACGCTGGC...AGCCGTATCGGAAGGTGCGGCT 13781 Root;Bacter...
#> [18]  1481 TGGAGAGTTTGATCCTGGCTCAG...NAACCGCAAGGATATAGCCGTC 13797 Root;Bacter...
#> [19]  1463 CGGCGTGCTTGGACCCACCCAAA...AGGAGGGTCCTAAGGTGGGGGC 13799 Root;Bacter...
#> [20]  1389 CGAGTGGCAAACGGGTGAGTAAC...AAACCGCAAGGATGCAGCCGTC 13800 Root;Bacter...
```

Note that the training data needs to have names in a specific RDP
format:

`"<ID> <Kingdom>;<Phylum>;<Class>;<Order>;<Family>;<Genus>"`

In the following we show the name for the first sequence. We use here
`sprintf` to display only the first 65~characters so the it fits into a
single line.

``` r

sprintf(names(trainingSequences[1]), fmt = "%.65s...")
#> [1] "13652 Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Peptococc..."
```

Now, we can train a the classifier. The model is stored in a directory
specified by the parameter `dir`.

``` r

customRDP <- trainRDP(trainingSequences, dir = "myRDP")
customRDP
#> RDPClassifier
#> Location: /home/runner/work/rRDP/rRDP/vignettes/myRDP
```

As test sequences, we sample 5 sequences form the training data and
extract a subsequences of length 500.

``` r

testSequences <- sample(trainingSequences, 5)

testSequences <- subseq(testSequences, start = 500, width = 500)
testSequences
#> DNAStringSet object of length 5:
#>     width seq                                               names               
#> [1]   500 CGTTGTCCGGAATTACTGGGCGT...TTCTCTGAAGGAGCTGTGAGACA 13763 Root;Bacter...
#> [2]   500 GGGTGAAATACCGCAGCTCAACT...TCGTGAGATGTTGGGTTAAGTCC 13677 Root;Bacter...
#> [3]   500 AAAACCTGGGCTCAACCGAGGGT...TTGGGTTAAGTCCCGCAACGAGC 13757 Root;Bacter...
#> [4]   500 CGTAGGCGGCTATAAAAGTCAGA...GGTTGTCGTCAGCTCGTGTCGTG 13762 Root;Bacter...
#> [5]   500 ATCAGCTCAACTGATAGCCTGCT...GGGTTAAGTCCCGCAACGAGCGC 13682 Root;Bacter...
```

``` r

pred <- predict(customRDP, testSequences)
pred
#>           domain     Phylum      Class         Order
#> 13763 Firmicutes Firmicutes Clostridia Clostridiales
#> 13677 Firmicutes Firmicutes Clostridia Clostridiales
#> 13757 Firmicutes Firmicutes Clostridia Clostridiales
#> 13762 Firmicutes Firmicutes Clostridia Clostridiales
#> 13682 Firmicutes Firmicutes Clostridia Clostridiales
#>                                                  Family                 Genus
#> 13763 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 13677                                    Peptococcaceae      Desulfotomaculum
#> 13757 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 13762 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 13682                                    Peptococcaceae      Desulfotomaculum
```

The prediction confidence is supplied as the attribute `"confidence"`.

``` r

attr(pred, "confidence")
#>       domain Phylum Class Order Family Genus
#> 13763      1      1     1     1      1     1
#> 13677      1      1     1     1      1     1
#> 13757      1      1     1     1      1     1
#> 13762      1      1     1     1      1     1
#> 13682      1      1     1     1      1     1
```

To evaluate the classification accuracy we can compare the known
classification with the predictions. The known classification was stored
in the FASTA file and encoded in RDP format. We can decode the
annotation using
[`decode_RDP()`](http://michael.hahsler.net/rRDP/reference/classification.md).
Greengenes format can be decoded using
[`decode_Greengenes()`](http://michael.hahsler.net/rRDP/reference/classification.md).

``` r

annotation <- names(testSequences)
actual <- decode_RDP(annotation)
actual
#>    Kingdom     Phylum      Class         Order
#> 1 Bacteria Firmicutes Clostridia Clostridiales
#> 2 Bacteria Firmicutes Clostridia Clostridiales
#> 3 Bacteria Firmicutes Clostridia Clostridiales
#> 4 Bacteria Firmicutes Clostridia Clostridiales
#> 5 Bacteria Firmicutes Clostridia Clostridiales
#>                                              Family                 Genus
#> 1 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 2                                    Peptococcaceae      Desulfotomaculum
#> 3 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 4 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 5                                    Peptococcaceae      Desulfotomaculum
#>   Species  Otu Org_name    Id
#> 1    <NA> <NA>     <NA> 13763
#> 2    <NA> <NA>     <NA> 13677
#> 3    <NA> <NA>     <NA> 13757
#> 4    <NA> <NA>     <NA> 13762
#> 5    <NA> <NA>     <NA> 13682
```

Now we can compare the prediction with the actual classification by
creating a confusion table and calculating the classification accuracy.
Here we do this at the Genus level.

``` r

confusionTable(actual, pred, rank = "genus")
#>                        predicted
#> actual                  Desulfotomaculum Thermoanaerobacterium
#>   Desulfotomaculum                     2                     0
#>   Thermoanaerobacterium                0                     3
accuracy(actual, pred, rank = "genus")
#> [1] 1
```

Note that the confidence and accuracy values are not very informative
since we trained the classifier with a tiny data set so the vignette can
be created fast. You can use a larger training set or a pretrained
classifier to get more informative results.

Since the custom classifier is stored on disc it can be recalled anytime
using [`rdp()`](http://michael.hahsler.net/rRDP/reference/rdp.md).

``` r

customRDP <- rdp(dir = "myRDP")
customRDP
#> RDPClassifier
#> Location: /home/runner/work/rRDP/rRDP/vignettes/myRDP
```

To permanently remove the classifier use
[`removeRDP()`](http://michael.hahsler.net/rRDP/reference/rdp.md). This
will delete the directory containing the classifier files.

``` r

removeRDP(customRDP)
```

RDP is trained with a 16S rRNA training set. The default classifier is
downloaded on first use and cached in the user’s R data directory. It
uses RDP Classifier 2.14 released in August 2023, with the bacterial and
archaeal taxonomy training set No. 19 (Wang and Cole 2024).

Since the download take a while, it is not shown executed in this
vignette.

``` r

default_16S_classifier <- rdp()
```

## Session Info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8          LC_NUMERIC=C             
#>  [3] LC_TIME=C.UTF-8           LC_COLLATE=C.UTF-8       
#>  [5] LC_MONETARY=C.UTF-8       LC_MESSAGES=C.UTF-8      
#>  [7] LC_PAPER=C.UTF-8          LC_NAME=C.UTF-8          
#>  [9] LC_ADDRESS=C.UTF-8        LC_TELEPHONE=C.UTF-8     
#> [11] LC_MEASUREMENT=C.UTF-8    LC_IDENTIFICATION=C.UTF-8
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] rRDP_2.0.0          Biostrings_2.80.2   Seqinfo_1.2.0      
#> [4] XVector_0.52.0      IRanges_2.46.0      S4Vectors_0.50.3   
#> [7] BiocGenerics_0.58.1 generics_0.1.4     
#> 
#> loaded via a namespace (and not attached):
#>  [1] crayon_1.5.3      cli_3.6.6         knitr_1.52        rlang_1.3.0      
#>  [5] xfun_0.61         otel_0.2.0        textshaping_1.0.5 rJava_1.0-18     
#>  [9] jsonlite_2.0.0    htmltools_0.5.9   ragg_1.5.2        sass_0.4.10      
#> [13] rmarkdown_2.32    evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
#> [17] yaml_2.3.12       lifecycle_1.0.5   compiler_4.6.1    fs_2.1.0         
#> [21] systemfonts_1.3.2 digest_0.6.39     R6_2.6.1          bslib_0.12.0     
#> [25] tools_4.6.1       pkgdown_2.2.1     cachem_1.1.0      desc_1.4.3
```

## Acknowledgments

This research is supported by research grant no. R21HG005912 from the
National Human Genome Research Institute (NHGRI / NIH).

## References

Cole, J. R., B. Chai, T. L. Marsh, et al. 2003. “The Ribosomal Database
Project (RDP-II): Previewing a New Autoaligner That Allows Regular
Updates and the New Prokaryotic Taxonomy.” *Nucleic Acids Research* 31
(1): 442–43. <https://doi.org/10.1093/nar/gkg039>.

Wang, Qiong, and James R. Cole. 2024. “Updated RDP Taxonomy and RDP
Classifier for More Accurate Taxonomic Classification.” *Microbiology
Resource Announcements* 0 (0): e01063–23.
<https://doi.org/10.1128/mra.01063-23>.

Wang, Qiong, George M Garrity, James M Tiedje, and James R Cole. 2007.
“Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into
the New Bacterial Taxonomy.” *Applied and Environmental Microbiology* 73
(16): 5261–67.
