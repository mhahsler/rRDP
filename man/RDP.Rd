\name{RDP}
\alias{rdp}
\alias{predict}
\alias{predict.RDPClassifier}
\alias{print.RDPClassifier}
\alias{trainRDP}
\alias{removeRDP}
\title{
Ribosomal Database Project (RDP) Classifier for 16S rRNA
}
\description{
Use the RDP classifier to classify 16S rRNA sequences. This package contains
currently RDP version 2.9.
}
\usage{
rdp(dir = NULL)
\method{predict}{RDPClassifier}(object, newdata,
  confidence=.8, rdp_args="", java_args="-Xmx1g", ...)
trainRDP(x, dir="classifier", rank="genus", java_args="-Xmx1g")
removeRDP(object)
}
\arguments{
  \item{dir}{directory where the classifier information is stored.} 
  \item{object}{a RDPClassifier object.}
  \item{newdata}{new data to be classified as a DNAStringSet.}
  \item{confidence}{numeric; minimum confidence level for classification.
    Results with lower confidence are replaced by NAs. Set to 0 to disable.}
  \item{rdp_args}{additional RDP arguments for classification
    (e.g., "-minWords 5" to set the 
    minimum number of words for each bootstrap trial.). 
    See RDP documentation.}
  \item{java_args}{additional arguments for java (default sets the 
    max. heap memory to 1GB).}
  \item{x}{an object of class DNAStringSet with the 16S rRNA sequences 
    for training.}
  \item{rank}{Taxonomic rank at which the classification is learned.}
  \item{...}{additional arguments (currently unused).}
}
\details{
RDP is a naive Bayes classifier using 8-mers as features.

\code{rdp()} creates a default classifier trained with the data shipped with RDP. 
Alternatively, a directory with the data for an existing classifier 
(created with \code{trainRDP()}) can be supplied.

\code{trainRDP()} creates a new classifier for the data in \code{x} and stores
the classifier information in \code{dir}. The data in
\code{x} needs to have annotations in the following format:

"<ID> <Kingdom>;<Phylum>;<Class>;<Order>;<Family>;<Genus>"

A created classifier can be removed with \code{removeRDP()}. This will remove
the directory which stores the classifier information.

The data for the default 16S rRNA classifier can be found in package 
\pkg{rRDPData}.
}
\value{
\code{rdp()} and \code{trainRDP()} return a RDPClassifier object.

\code{predict()} returns a data.frame containing the classification 
results for each sequence
(rows). The data.frame has an
attribure called "confidence" with a matrix containing the confidence values.
}
\references{
RDP Classifier \url{http://sourceforge.net/projects/rdp-classifier/}

Qiong Wang, George M. Garrity, James M. Tiedje and James R. Cole. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy, Appl. Environ. Microbiol. August 2007 vol. 73 no. 16 5261-5267. 
}
\examples{
### Use the default classifier
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="rRDP"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## use rdp for classification (this needs package rRDPData) 
pred <- predict(rdp(), seq)
pred
  
attr(pred, "confidence")  

### Train a custom RDP classifier on new data
trainingSequences <- readDNAStringSet(
    system.file("examples/trainingSequences.fasta", package="rRDP"))

customRDP <- trainRDP(trainingSequences)
customRDP

testSequences <- readDNAStringSet(
    system.file("examples/testSequences.fasta", package="rRDP"))
predict(customRDP, testSequences)

## clean up
removeRDP(customRDP)
}
\keyword{model}
