\name{accuracy}
\alias{confusionTable}
\alias{accuracy}
\title{Calculate Classification Accuracy}
\description{
Calculate the classification accuracy at a given phylogenetic level.
}
\usage{
accuracy(actual, predicted, rank)
confusionTable(actual, predicted, rank)
}
\arguments{
        \item{actual}{data.frame with the actual classification hierarchy.}
        \item{predicted}{data.frame with the predicted classification hierarchy.}
        \item{rank}{rank at which the accuracy should be evaluated.}
}
%\details{
%}
\value{
The accuracy or a confusion table.
}
%\references{
%}
\examples{
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="rRDP"))

### decode the actual classification
actual <- decode_Greengenes(names(seq))
  
### use RDP to predict the classification
pred <- predict(rdp(), seq)  

### calculate accuracy
confusionTable(actual, pred, "genus")  
accuracy(actual, pred, "genus")
}
\keyword{model}