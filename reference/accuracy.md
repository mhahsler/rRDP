# Calculate Classification Accuracy

Calculate the classification accuracy at a given phylogenetic level.

## Usage

``` r
accuracy(actual, predicted, rank)

confusionTable(actual, predicted, rank)
```

## Arguments

- actual:

  data.frame with the actual classification hierarchy.

- predicted:

  data.frame with the predicted classification hierarchy.

- rank:

  rank at which the accuracy should be evaluated.

## Value

The accuracy or a confusion table.

## Examples

``` r

if (FALSE) { # \dontrun{
# the download may take a while
if (!file.exists("data")) {
  if (!file.exists("data.tgz"))
    download.file(paste0("https://master.dl.sourceforge.net/project/",
      "rdp-classifier/rdp-classifier/data.tgz"), "data.tgz", mode = "wb")
  untar("data.tgz")
}

classifier <- rdp("data/classifier/16srrna")

seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
    package = "rRDP"
))

### decode the actual classification
actual <- decode_Greengenes(names(seq))

### use RDP to predict the classification
pred <- predict(classifier, seq)

### calculate accuracy
confusionTable(actual, pred, "genus")
accuracy(actual, pred, "genus")
} # }
```
