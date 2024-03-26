#######################################################################
# rRDP - Interfaces RDP
# Copyright (C) 2014 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


#' Calculate Classification Accuracy
#'
#' Calculate the classification accuracy at a given phylogenetic level.
#'
#' @name accuracy
#' @aliases confusionTable accuracy
#' @param actual data.frame with the actual classification hierarchy.
#' @param predicted data.frame with the predicted classification hierarchy.
#' @param rank rank at which the accuracy should be evaluated.
#' @return The accuracy or a confusion table.
#' @keywords model
#' @examples
#'
#' seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
#'     package = "rRDP"
#' ))
#'
#' ### decode the actual classification
#' actual <- decode_Greengenes(names(seq))
#'
#' ### use RDP to predict the classification
#' pred <- predict(rdp(), seq)
#'
#' ### calculate accuracy
#' confusionTable(actual, pred, "genus")
#' accuracy(actual, pred, "genus")
#' @export
accuracy <- function(actual, predicted, rank) {
    ct <- confusionTable(actual, predicted, rank)
    sum(diag(ct)) / sum(ct)
}

#' @rdname accuracy
#' @export
confusionTable <- function(actual, predicted, rank) {
    ranka <- colnames(actual)[grep(rank, colnames(actual), ignore.case = TRUE)]
    rankp <- colnames(predicted)[grep(rank, colnames(predicted), ignore.case = TRUE)]
    actual <- factor(actual[, ranka])
    predicted <- factor(predicted[, rankp], exclude = NULL)

    commonLevels <- sort(unique(c(levels(actual), levels(predicted))), na.last = TRUE)
    actual <- factor(actual, levels = commonLevels, exclude = NULL)
    predicted <- factor(predicted, levels = commonLevels, exclude = NULL)

    table(actual, predicted, dnn = list("actual", "predicted"))
}
