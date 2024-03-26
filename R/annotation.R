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


#' Decoding and Encoding Phylogenetic Classification Annotations
#'
#' Functions to represent, decode and encode phylogenetic classification
#' annotations used in FASTA files by RDP and the Greengenes project.
#'
#' @name classification
#' @aliases GenClass16S decode_RDP encode_RDP decode_Greengenes
#' encode_Greengenes
#' @param Kingdom Name of the kingdom to which the organism belongs.
#' @param Phylum Name of the phylum to which the organism belongs.
#' @param Class Name of the class to which the organism belongs.
#' @param Order Name of the order to which the organism belongs.
#' @param Family Name of the family to which the organism belongs.
#' @param Genus Name of the genus to which the organism belongs.
#' @param Species Name of the species to which the organism belongs.
#' @param Otu Name of the otu to which the organism belongs.
#' @param Org_name Name of the organism.
#' @param Id ID of the sequence.
#' @param annotation Annotation from a FASTA file containing the classification
#' information.
#' @param classification A \code{data.frame} created with \code{GenClass16S()}
#' with the classification information.
#' @return \code{GenClass16S()} and \code{decodeX()} return a
#' \code{data.frame}. \code{encodeX()} returns a string with the corresponding
#' annotation.
#' @keywords model
#' @examples
#'
#' seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
#'     package = "rRDP"
#' ))
#'
#' ### the FASTA annotation is read as names. This data has a Greengenes format
#' ### annotation
#' names(seq)
#'
#' classification <- decode_Greengenes(names(seq))
#' classification
#'
#' ### look at the Genus of all sequences
#' classification[, "Genus"]
#'
#' ### to train the RDP classifier, the annotations need to be in RDP format
#' annotation <- encode_RDP(classification)
#' names(seq) <- annotation
#' seq
#'
#' ### now we can train the classifier
#' customRDP <- trainRDP(seq)
#' customRDP
#'
#' ## clean up
#' removeRDP(customRDP)
#' @export
decode_Greengenes <- function(annotation) {
    fields <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "otu_")
    ## add out and org_name

    # remove leading ">"
    annotation <- sub(">", "", annotation)
    # split at "k__"
    tmp <- strsplit(annotation, " *k__")
    org_name <- gsub("'", "", sapply(tmp, "[", 1))
    ss <- strsplit(org_name, " ")
    id <- sapply(ss, "[", 1)
    org_name <- sapply(ss, FUN = function(s) paste(s[-1], collapse = "_"))

    tmp <- sapply(tmp, "[", 2)
    tmp <- paste("k__", tmp, sep = "")
    tmp <- strsplit(tmp, "(?= .__)|(?= otu_)", perl = TRUE)


    cl <- sapply(fields, FUN = function(f) {
        val <- lapply(tmp, FUN = function(i) grep(f, i, value = TRUE))
        sapply(val, FUN = function(s) if (length(s) == 1) sub("(^.__)|(otu_)", "", s) else "unknown")
    }, simplify = FALSE)

    # remove trailing ;
    cl <- lapply(cl, FUN = function(s) sub(";$", "", s))


    return(GenClass16S(
        cl[[1]], cl[[2]], cl[[3]], cl[[4]], cl[[5]],
        cl[[6]], cl[[7]], cl[[8]], org_name, id
    ))
}

## classification hierarchy for 16S
#' @rdname classification
#' @export
GenClass16S <- function(Kingdom = NA, Phylum = NA, Class = NA, Order = NA,
                        Family = NA, Genus = NA, Species = NA, Otu = NA, Org_name = NA, Id = NA) {
    ### prevent recycling
    params <- as.list(environment())
    l <- max(sapply(params, length))
    params <- lapply(params, "length<-", l)
    return(as.data.frame(do.call(cbind, params)))
}

#' @rdname classification
#' @export
encode_Greengenes <- function(classification) {
    return(paste(
        ">", classification[, "Id"], " ", classification[, "Org_name"],
        " k__", classification[, "Kingdom"], ";",
        " p__", classification[, "Phylum"], ";",
        " c__", classification[, "Class"], ";",
        " o__", classification[, "Order"], ";",
        " f__", classification[, "Family"], ";",
        " g__", classification[, "Genus"], ";",
        " s__", classification[, "Species"], ";",
        " otu_", classification[, "Otu"],
        sep = ""
    ))
}

#' @rdname classification
#' @export
decode_RDP <- function(annotation) {
    ann <- strsplit(annotation, ";")

    ret <- matrix(NA_character_, ncol = length(GenClass16S()), nrow = length(ann))
    for (i in seq_along(ann)) ret[i, seq_along((ann[[i]]))] <- ann[[i]]

    ret[, 1] <- sub(" Root$", "", ret[, 1])

    GenClass16S(
        Kingdom = as.character(ret[, 2]),
        Phylum = as.character(ret[, 3]),
        Class = as.character(ret[, 4]),
        Order = as.character(ret[, 5]),
        Family = as.character(ret[, 6]),
        Genus = as.character(ret[, 7]),
        Species = NA_character_,
        Otu = NA_character_,
        Org_name = NA_character_,
        Id = as.character(ret[, 1])
    )
}

#' @rdname classification
#' @export
encode_RDP <- function(classification) {
    h <- classification[, seq_len(6)]
    h <- apply(h, MARGIN = 1, FUN = function(x) paste(x, collapse = ";"))
    h <- gsub(";unknown", "", h)
    h <- gsub(" \\(class\\)", "", h)

    paste(classification[, "Id"], " ", "Root;", h, sep = "")
}
