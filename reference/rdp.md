# Ribosomal Database Project (RDP) Classifier for 16S rRNA

RDP is a naive Bayes classifier for biological sequences using 8-mers as
features. Use the RDP classifier (Wang et al, 2007) to classify 16S rRNA
sequences. This package contains currently RDP version 2.14 released in
August 2023.

## Usage

``` r
rdp(dir = NULL)

# S3 method for class 'RDPClassifier'
predict(object, newdata, confidence = 0.8, rdp_args = "", verbose = FALSE, ...)

trainRDP(x, dir, rank = "genus", verbose = FALSE)

removeRDP(object)
```

## Arguments

- dir:

  directory where the classifier information is stored.

- object:

  a RDPClassifier object.

- newdata:

  new data to be classified as a
  [Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html).

- confidence:

  numeric; minimum confidence level for classification. Results with
  lower confidence are replaced by NAs. Set to 0 to disable.

- rdp_args:

  additional RDP arguments for classification (e.g., `"-minWords 5"` to
  set the minimum number of words for each bootstrap trial.). See RDP
  documentation.

- verbose:

  logical; print additional information.

- ...:

  additional arguments (currently unused).

- x:

  an object of class
  [Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  with the 16S rRNA sequences for training.

- rank:

  Taxonomic rank at which the classification is learned.

## Value

`rdp()` and `trainRDP()` return a `RDPClassifier` object.

`predict()` returns a data.frame containing the classification results
for each sequence (rows). The data.frame has an attribute called
`"confidence"` with a matrix containing the confidence values.

## Details

### Train a new classifier

`trainRDP()` creates a new classifier for the data in `x` and stores the
classifier information in `dir`. The data in `x` needs to have
annotations in the following format:

`"<ID> <Kingdom>;<Phylum>;<Class>;<Order>;<Family>;<Genus>"`

The trained classifier can be used as long as the directory is
available. A created classifier can be removed with `removeRDP()`. This
will remove the directory which stores the classifier information
permanently.

### Use a pretrained classifier

Calling `rdp()` without a directory downloads the default 16S classifier
from the RDP project and caches it in the user data directory. The first
download may take a while, and you will be asked before it starts. A
custom classifier can be loaded by passing its directory to `rdp(dir)`.

### Classify sequences

`rdp()` loads a classifier model from a specified directory. `predict()`
can then be used to create a prediction.

## References

Hahsler M, Nagar A (2020). "rRDP: Interface to the RDP Classifier." R
Package, Bioconductor.
[doi:10.18129/B9.bioc.rRDP](https://doi.org/10.18129/B9.bioc.rRDP) .

RDP classifier software:
<https://sourceforge.net/projects/rdp-classifier/>

Qiong Wang, George M. Garrity, James M. Tiedje and James R. Cole. Naive
Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New
Bacterial Taxonomy, Appl. Environ. Microbiol. August 2007 vol. 73 no. 16
5261-5267.
[doi:10.1128/AEM.00062-07](https://doi.org/10.1128/AEM.00062-07)

Qiong W. and Cole J.R. Updated RDP taxonomy and RDP Classifier for more
accurate taxonomic classification, Microbial Ecology, Announcement, 4
March 2024.
[doi:10.1128/mra.01063-23](https://doi.org/10.1128/mra.01063-23)

## Examples

``` r
## Train a RDP classifier on new data
trainingSequences <- readDNAStringSet(
    system.file("examples/trainingSequences.fasta", package = "rRDP")
)

customRDP <- trainRDP(trainingSequences, dir = "test_classifier")
customRDP
#> RDPClassifier
#> Location: /home/runner/work/rRDP/rRDP/docs/reference/test_classifier 

testSequences <- readDNAStringSet(
    system.file("examples/testSequences.fasta", package = "rRDP")
)
predict(customRDP, testSequences)
#>           domain     Phylum      Class         Order
#> 13811 Firmicutes Firmicutes Clostridia Clostridiales
#> 13813 Firmicutes Firmicutes Clostridia Clostridiales
#> 13678 Firmicutes Firmicutes Clostridia Clostridiales
#> 13755 Firmicutes Firmicutes Clostridia Clostridiales
#> 13661 Firmicutes Firmicutes Clostridia Clostridiales
#>                                                  Family                 Genus
#> 13811                                   Veillonellaceae           Selenomonas
#> 13813                                   Veillonellaceae           Selenomonas
#> 13678                                    Peptococcaceae      Desulfotomaculum
#> 13755 Thermoanaerobacterales Family III. Incertae Sedis Thermoanaerobacterium
#> 13661                                    Peptococcaceae      Desulfotomaculum

## clean up (don't do this if you want to use the classifier later on again)
removeRDP(customRDP)

## Download a pretrained model manually
if (FALSE) { # \dontrun{
# the download may take a while

if (!file.exists("data")) {
  if (!file.exists("data.tgz"))
    download.file(paste0("https://downloads.sourceforge.net/",
      "project/rdp-classifier/rdp-classifier/data.tgz"), 
      "data.tgz", mode = "wb")
  untar("data.tgz")
}
  
classifier <- rdp("data/classifier/16srrna")
classifier

seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
    package = "rRDP"
))

# decode the actual classification
actual <- decode_Greengenes(names(seq))

# use RDP to predict the classification
pred <- predict(classifier, seq)

# calculate accuracy
confusionTable(actual, pred, "genus")
accuracy(actual, pred, "genus")

# Don't forget to delete the models if you do not need them anymore
} # }
```
