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

.isRDP <- function(dir) {
  is.null(dir) ||
    file.exists(file.path(dir, "wordConditionalProbIndexArr.txt"))
}

.defaultRDPPath <- function() {
  root <- tools::R_user_dir("rRDP", "data")
  file.path(root, "data", "classifier", "16srrna")
}

.defaultRDP <- function() {
  root <- tools::R_user_dir("rRDP", "data")
  dir <- .defaultRDPPath()
  if (!.isRDP(dir)) {
    answer <- utils::askYesNo(
      "The default RDP classifier files need to be downloaded. This may take a while. Download now?",
      default = FALSE
    )
    if (!isTRUE(answer)) {
      stop("The default RDP classifier was not downloaded. Call rdp() again when you are ready to download it.")
    }
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    archive <- tempfile(fileext = ".tgz")
    on.exit(unlink(archive))
    url <- paste0("https://master.dl.sourceforge.net/project/",
                  "rdp-classifier/rdp-classifier/data.tgz")
    cat("Downloading the default classifier data\n")
    utils::download.file(url, archive, mode = "wb", quiet = FALSE)
    files <- utils::untar(archive, list = TRUE)
    files <- files[startsWith(files, "data/classifier/16srrna/")]
    if (!length(files)) stop("The RDP archive does not contain the 16S classifier.")
    utils::untar(archive, files = files, exdir = root)
  }
  normalizePath(dir)
}

#' Ribosomal Database Project (RDP) Classifier for 16S rRNA
#'
#' RDP is a naive Bayes classifier for biological sequences using 8-mers as features.
#' Use the RDP classifier (Wang et al, 2007) to classify 16S rRNA sequences.
#' This package contains
#' currently RDP version 2.14 released in August 2023.
#'
#' ## Train a new classifier
#' `trainRDP()` creates a new classifier for the data in `x` and
#' stores the classifier information in `dir`. The data in `x` needs
#' to have annotations in the following format:
#'
#' `"<ID> <Kingdom>;<Phylum>;<Class>;<Order>;<Family>;<Genus>"`
#'
#' The trained classifier can be used as long as the directory is available.
#' A created classifier can be removed with `removeRDP()`. This will
#' remove the directory which stores the classifier information permanently.
#'
#' ## Use a pretrained classifier
#' Calling `rdp()` without a directory downloads the default 16S classifier
#' from the RDP project and caches it in the user data directory. The first
#' download may take a while, and you will be asked before it starts. A custom classifier can be loaded by passing its
#' directory to `rdp(dir)`.
#'
#' ## Classify sequences
#' `rdp()` loads a classifier model from a specified directory. `predict()` can then
#' be used to create a prediction.
#'
#' @name rdp
#' @aliases rdp RDP predict predict.RDPClassifier print.RDPClassifier trainRDP removeRDP
#' @param dir directory where the classifier information is stored.
#' @param object a RDPClassifier object.
#' @param newdata new data to be classified as a [Biostrings::DNAStringSet].
#' @param confidence numeric; minimum confidence level for classification.
#' Results with lower confidence are replaced by NAs. Set to 0 to disable.
#' @param rdp_args additional RDP arguments for classification (e.g.,
#' `"-minWords 5"` to set the minimum number of words for each bootstrap trial.).
#' See RDP documentation.
#' @param x an object of class [Biostrings::DNAStringSet] with the 16S rRNA sequences for
#' training.
#' @param rank Taxonomic rank at which the classification is learned.
#' @param verbose logical; print additional information.
#' @param ... additional arguments (currently unused).
#' @return `rdp()` and `trainRDP()` return a `RDPClassifier` object.
#'
#' `predict()` returns a data.frame containing the classification results
#' for each sequence (rows). The data.frame has an attribute called
#' `"confidence"` with a matrix containing the confidence values.
#' @references Hahsler M, Nagar A (2020). "rRDP: Interface to the RDP
#' Classifier." R Package, Bioconductor.
#' \doi{10.18129/B9.bioc.rRDP}.
#'
#' RDP classifier software:
#' \url{https://sourceforge.net/projects/rdp-classifier/}
#'
#' Qiong Wang, George M. Garrity, James M. Tiedje and James R. Cole. Naive
#' Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New
#' Bacterial Taxonomy, Appl. Environ. Microbiol. August 2007 vol. 73 no. 16
#' 5261-5267. \doi{10.1128/AEM.00062-07}
#'
#' Qiong W. and Cole J.R. Updated RDP taxonomy and RDP Classifier for more
#' accurate taxonomic classification, Microbial Ecology,
#' Announcement, 4 March 2024. \doi{https://doi.org/10.1128/mra.01063-23}
#'
#' @keywords model
#' @examples
#' ## Train a RDP classifier on new data
#' trainingSequences <- readDNAStringSet(
#'     system.file("examples/trainingSequences.fasta", package = "rRDP")
#' )
#'
#' customRDP <- trainRDP(trainingSequences, dir = "test_classifier")
#' customRDP
#'
#' testSequences <- readDNAStringSet(
#'     system.file("examples/testSequences.fasta", package = "rRDP")
#' )
#' predict(customRDP, testSequences)
#'
#' ## clean up (don't do this if you want to use the classifier later on again)
#' removeRDP(customRDP)
#' 
#' ## Download a pretrained model manually
#' \dontrun{
#' # the download may take a while
#' 
#' if (!file.exists("data")) {
#'   if (!file.exists("data.tgz"))
#'     download.file(paste0("https://master.dl.sourceforge.net/project/",
#'       "rdp-classifier/rdp-classifier/data.tgz"), "data.tgz", mode = "wb")
#'   untar("data.tgz")
#' }
#'   
#' classifier <- rdp("data/classifier/16srrna")
#' classifier
#' 
#' seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
#'     package = "rRDP"
#' ))
#' 
#' # decode the actual classification
#' actual <- decode_Greengenes(names(seq))
#'
#' # use RDP to predict the classification
#' pred <- predict(classifier, seq)
#'
#' # calculate accuracy
#' confusionTable(actual, pred, "genus")
#' accuracy(actual, pred, "genus")
#' 
#' # Don't forget to delete the models if you do not need them anymore
#' }
#' @export
rdp <- function(dir = NULL) {
  if (is.null(dir)) dir <- .defaultRDP()
  dir <- normalizePath(dir)
  
  if (!.isRDP(dir)) {
    stop("Not a RDP classifier directory!")
  }
  structure(list(dir = dir), class = "RDPClassifier")
}

#' @export
print.RDPClassifier <- function(x, ...) {
  loc <- x$dir
  if (is.null(loc)) {
    loc <- "Default RDP classifier"
  }
  cat("RDPClassifier\nLocation:", loc, "\n")
  if (!.isRDP(x$dir)) {
    cat("The RDPClassifier is not valid!\n")
  }
}

#' @rdname rdp
#' @export
predict.RDPClassifier <- function(object,
                                  newdata,
                                  confidence = .8,
                                  rdp_args = "",
                                  verbose = FALSE,
                                  ...) {
  classifier <- object$dir
  x <- newdata
  
  # get temp files and change working directory
  wd <- tempdir()
  dir <- getwd()
  temp_file <- tempfile(tmpdir = wd)
  #  temp_file <- "newdata"
  on.exit({
    file.remove(Sys.glob(paste(temp_file, "*", sep = "")))
    setwd(dir)
  })
  
  # setwd(wd)
  infile <- paste(temp_file, ".fasta", sep = "")
  outfile <- paste(temp_file, "_tax_assignments.txt", sep = "")
  
  ## property?
  fp <- file.path(classifier, "rRNAClassifier.properties")
  if (.Platform$OS.type == "windows") {
    fp <- utils::shortPathName(fp)
  }
  if (!is.null(classifier)) {
    property <- c("-t", fp)
  } else {
    property <- NULL
  }
  
  writeXStringSet(x, infile, append = FALSE)
  
  cmd <- c(
    "classify",
    property,
    "-o",
    outfile,
    "-c",
    confidence,
    "-format",
    "fixrank",
    unlist(strsplit(rdp_args, split = "\\s+")),
    infile
  )
  
  if (verbose) {
    cat("Calling RDP with arguments: ", paste(cmd), "\n")
  }
  
  rdp_obj <- .jnew("edu.msu.cme.rdp.classifier.cli.ClassifierMain")
  s <- .jcall(rdp_obj, returnSig = "V", method = "main", cmd)
  
  ## run the garbage collector or Java on Windows will not close the files.
  remove(rdp_obj)
  rJava::.jgc()
  
  ## read and parse rdp output
  cl_tab <- read.table(outfile, sep = "\t")
  
  if (verbose) {
    cat("Number of hits returned br RDP ", nrow(cl_tab), "\n")
  }
  
  seq_names <- cl_tab[, 1] ## sequence names are in first column
  
  i <-
    seq(3, ncol(cl_tab), by = 3) ## start at col. 3 and use 3 columns for each tax. level
  
  ## get classification
  cl <- cl_tab[, i]
  dimnames(cl) <- list(seq_names, as.matrix(cl_tab[1, i + 1])[1, ])
  
  ## get confidence
  conf <- as.matrix(cl_tab[, i + 2])
  dimnames(conf) <- list(seq_names, as.matrix(cl_tab[1, i + 1])[1, ])
  
  ### don't show assignments with too low confidence
  if (confidence > 0) {
    cl[conf < confidence] <- NA
  }
  
  attr(cl, "confidence") <- conf
  cl
}

#' @rdname rdp
#' @export
trainRDP <-
  function(x,
           dir,
           rank = "genus",
           verbose = FALSE) {
    if (file.exists(dir)) {
      stop(
        "Classifier directory",
        sQuote(dir),
        "already exists! Choose a different directory, use removeRDP(), or remove the directory manually."
      )
    }
    
    taxonomyNames <-
      c("Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species")
    rankNumber <- pmatch(tolower(rank), tolower(taxonomyNames))
    dir.create(dir)
    l <- strsplit(names(x), "Root;")
    # annot is the hierarchy starting from kingdom down to genus (or desired rank)
    annot <- sapply(
      l,
      FUN = function(x) {
        x[2]
      }
    )
    # Make names of x so that it only has hierarchy info upto the desired rank
    names(x) <- sapply(
      strsplit(names(x), ";"),
      FUN = function(y) {
        paste(y[seq_len(rankNumber + 1)], collapse = ";")
      }
    )
    # get the indices of sequences that are to be removed because they have incomplete hierarchy information
    removeIdx <- as.integer(sapply(
      annot,
      FUN = function(y) {
        if (length(unlist(strsplit(y, ";"))) < rankNumber ||
            grepl(";;", y) ||
            grepl("unknown", y)) {
          1
        } else {
          0
        }
      }
    ))
    
    if (any(removeIdx == 1)) {
      removeIdx <- which(removeIdx == 1)
      # get greengenes ids to be removed
      idsRemoved <- as.character(sapply(
        names(x),
        FUN = function(y) {
          as.character(unlist(strsplit(y, " "))[1])
        }
      ))[removeIdx]
      x <- x[-removeIdx]
      annot <- annot[-removeIdx]
      cat(
        "Warning! The following sequences did not contain complete hierarchy information and have been ignored:",
        idsRemoved,
        "\n"
      )
    }
    if (length(x) <= 0) {
      stop("No sequences with complete information found")
    }
    # writeXStringSet(x,file.path(dir,"train.fasta"))
    h <- matrix(ncol = rankNumber, nrow = 0)
    # colnames(h) <-c("Kingdom","Phylum","Class","Order","Family","Genus")
    colnames(h) <- taxonomyNames[seq_len(rankNumber)]
    for (i in seq_along(annot)) {
      h <- rbind(h, unlist(strsplit(annot[i], ";"))[seq_len(rankNumber)])
    }
    m <- matrix(ncol = 5, nrow = 0)
    # first row of the file
    f <- "0*Root*-1*0*rootrank"
    m <- rbind(m, unlist(strsplit(f, split = "\\*")))
    taxNames <- colnames(h)
    badHierarchy <- vector()
    for (i in seq_len(nrow(h))) {
      for (j in seq_len(ncol(h))) {
        taxId <- nrow(m)
        taxonName <- h[i, j]
        if (j == 1) {
          parentTaxId <- 0
        } else {
          parentTaxId <- previousTaxId
        }
        depth <- j
        rank <- colnames(h)[j]
        
        # search if already there
        if (length(which(m[, 2] == taxonName &
                         m[, 5] == rank)) == 0) {
          str <- paste(taxId, taxonName, parentTaxId, depth, rank, sep = "*")
          m <- rbind(m, unlist(strsplit(str, split = "\\*")))
          previousTaxId <- taxId
        } else if (length(which(m[, 2] == taxonName &
                                m[, 5] == rank)) > 0) {
          row <- which(m[, 2] == taxonName & m[, 5] == rank)
          # error check -> is parent tax id same for both or not?
          # if not, remove the sequence
          if (parentTaxId != m[row, 3]) {
            # x<-x[-i]
            badHierarchy <- c(badHierarchy, i)
          }
          previousTaxId <- m[row, 1]
        }
        # end seach
      }
    }
    
    out <- apply(
      m,
      MARGIN = 1,
      FUN = function(x) {
        paste(x, collapse = "*")
      }
    )
    write(out, file = file.path(dir, "train.txt"))
    # create parsed training files
    if (length(badHierarchy) > 0) {
      warning(
        "Following sequences had bad sequence hierarchy information, so removing: ",
        names(x)[badHierarchy],
        "\n"
      )
      x <- x[-badHierarchy]
    }
    writeXStringSet(x, file.path(dir, "train.fasta"))
    
    cmd <- c(
      "train",
      "-o",
      dir,
      "-t",
      file.path(dir, "train.txt"),
      "-s",
      file.path(dir, "train.fasta")
    )
    
    if (verbose) {
      cat("Calling RDP with arguments: ", paste(cmd), "\n")
    }
    
    rdp_obj <- .jnew("edu.msu.cme.rdp.classifier.cli.ClassifierMain")
    s <- .jcall(rdp_obj, returnSig = "V", method = "main", cmd)
    
    ## run the garbage collector or Java on Windows will not close the files.
    remove(rdp_obj)
    rJava::.jgc()
    
    file.copy(system.file("java/rRNAClassifier.properties", package = "rRDP"),
              dir)
    
    rdp(dir)
  }

#' @rdname rdp
#' @export
removeRDP <- function(object) {
  ### first check if it looks like a RDP directory!
  if (!dir.exists(object$dir)) {
    stop("The given RDPClassifier/directory does not exits!")
  }
  
  if (!.isRDP(object$dir)) {
    stop("The given RDPClassifier/directory is not a valid RDP classifier!")
  }
  
  ### don't remove the default data
  
  if (identical(normalizePath(object$dir), normalizePath(.defaultRDPPath(), mustWork = FALSE))) {
    stop("The default classifier is cached for reuse and cannot be removed with removeRDP().")
  }
  status <- unlink(object$dir, recursive = TRUE)
  if (status != 0) {
    stop(
      "Removing the RDPClassifier/directory failed! ",
      "Please remove the following directory manually: ",
      object$dir
    )
  }
  
}
