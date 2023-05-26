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

.get_rdp <- function() {
  p <- system.file("java/rdp_classifier.jar", package="rRDP")
  
  if (.Platform$OS.type == "windows")
    p <- utils::shortPathName(p)
  
  return(p)
}

  
  
.isRDP <- function(dir) is.null(dir) ||
  file.exists(file.path(dir, "wordConditionalProbIndexArr.txt"))

### NULL is the default classifier in package rRDPData
rdp <- function(dir = NULL) {
  if(is.null(dir)) {
    dir <- system.file("extdata","16srrna", package="rRDPData")
    if(dir=="") stop("Install package 'rRDPData' for the trained 16S rRNA classifier.")
  }
  dir <- normalizePath(dir)

  if(!.isRDP(dir)) stop("Not a RDP classifier directory!")
  structure(list(dir = dir), class="RDPClassifier")
}

print.RDPClassifier <- function(x, ...) {
  loc <- x$dir
  if(is.null(loc)) loc <- "Default RDP classifier"
  cat("RDPClassifier\nLocation:", loc, "\n")
}

# see: https://stackoverflow.com/questions/43528059/where-does-rjava-takes-location-of-jdk-from
.javaExecutable <- function() {
    java <- ""
  
    ### Linux (set with R CMD reconfigure Java)
    if (file.exists(file.path(R.home(component="etc"), "Makeconf")))
      java <- strsplit(grep("^JAVA ", 
        readLines(file.path(R.home(component="etc"), "Makeconf")), value = TRUE), 
        split = "\\s*=\\s*")[[1L]][2L]

    ### Windows - fall back to env
    if (.Platform$OS.type == "windows")
      java <- utils::shortPathName(file.path(Sys.getenv("JAVA_HOME"), "bin", "java.exe"))
    
    if (java == "")
      stop("Install the Java JDK and run 'R CMD javareconf' or set the environment variable 'JAVA_HOME'.")

    return(java)
}
    
predict.RDPClassifier <- function(object, newdata,
  confidence=.8, rdp_args="", java_args="-Xmx1g", ...) {

  classifier <- object$dir
  x <- newdata

  # get temp files and change working directory
  wd <- tempdir()
  dir <- getwd()
  temp_file <- tempfile(tmpdir = wd)
#  temp_file <- "newdata"
  on.exit({
    file.remove(Sys.glob(paste(temp_file, "*", sep="")))
    setwd(dir)
  })

  #setwd(wd)
  infile <- paste(temp_file, ".fasta", sep="")
  outfile <- paste(temp_file, "_tax_assignments.txt", sep="")
  
  ## property?
  fp <- file.path(classifier,"rRNAClassifier.properties")
  if (.Platform$OS.type == "windows")
    fp <- utils::shortPathName(fp)
  if(!is.null(classifier)) property <- paste("-t", fp)
  else property <- ""
  
  writeXStringSet(x, infile, append=FALSE)
  
  exe <- paste(.javaExecutable(), java_args, "-jar", .get_rdp(), "classify",
               property, "-o", outfile, "-c", confidence, "-format fixrank", rdp_args, infile)
  if (system(exe, ignore.stdout=TRUE, ignore.stderr=TRUE, intern=FALSE))
    stop("Error executing: ", exe)

  ## read and parse rdp output
  cl_tab <- read.table(outfile, sep="\t")

  seq_names <- cl_tab[,1] ## sequence names are in first column

  i <- seq(3, ncol(cl_tab), by=3) ## start at col. 3 and use 3 columns for each tax. level

  ## get classification
  cl <- cl_tab[,i]
  dimnames(cl) <- list(seq_names, as.matrix(cl_tab[1,i+1])[1,])

  ## get confidence
  conf <- as.matrix(cl_tab[,i+2])
  dimnames(conf) <- list(seq_names, as.matrix(cl_tab[1,i+1])[1,])

  ### don't show assignments with too low confidence
  if(confidence>0) cl[conf < confidence] <- NA

  attr(cl, "confidence") <- conf
  cl
}

trainRDP <- function(x, dir="classifier", rank="genus", java_args="-Xmx1g")
{
  if (file.exists(dir)) stop("Classifier directory already exists! Choose a different directory or use removeRDP().")

  taxonomyNames <-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rankNumber <- pmatch(tolower(rank),tolower(taxonomyNames))
  dir.create(dir)
  l<-strsplit(names(x),"Root;")
  # annot is the hierarchy starting from kingdom down to genus (or desired rank)
  annot<-sapply(l,FUN=function(x) x[2])
  # Make names of x so that it only has hierarchy info upto the desired rank
  names(x)<-sapply(strsplit(names(x),";"), FUN=function(y) {
    paste(y[1:(rankNumber+1)],collapse=";")
  })
  #get the indices of sequences that are to be removed because they have incomplete hierarchy information
  removeIdx <- as.integer(sapply(annot, FUN=function(y) {
    if(length(unlist(strsplit(y,";")))<rankNumber || grepl(";;",y)
      || grepl("unknown",y)) 1 else 0
  }))

  if(any(removeIdx==1)) {
    removeIdx <- which(removeIdx==1)
    #get greengenes ids to be removed
    idsRemoved <- as.character(sapply(names(x),
      FUN=function(y) as.character(unlist(strsplit(y," "))[1])))[removeIdx]
    x<-x[-removeIdx]
    annot<-annot[-removeIdx]
    cat("Warning! The following sequences did not contain complete hierarchy information and have been ignored:", idsRemoved,"\n")
  }
  if (length(x) <=0) stop("No sequences with complete information found")
  #writeXStringSet(x,file.path(dir,"train.fasta"))
  h<-matrix(ncol=rankNumber,nrow=0)
  #colnames(h) <-c("Kingdom","Phylum","Class","Order","Family","Genus")
  colnames(h) <- taxonomyNames[1:rankNumber]
  for(i in 1:length(annot)) {
    h<-rbind(h,unlist(strsplit(annot[i],";"))[1:rankNumber])
  }
  m<-matrix(ncol=5,nrow=0)
  #first row of the file
  f<-"0*Root*-1*0*rootrank"
  m<-rbind(m,unlist(strsplit(f,split="\\*")))
  taxNames <- colnames(h)
  badHierarchy<-vector()
  for(i in 1:nrow(h)) {
    for(j in 1:ncol(h))	{
      taxId <- nrow(m)
      taxonName <- h[i,j]
      if (j==1)
        parentTaxId <-0
      else
        parentTaxId <- previousTaxId
      depth <- j
      rank <- colnames(h)[j]

      #search if already there
      if (length(which(m[,2]==taxonName & m[,5]==rank)) == 0) {
        str <- paste(taxId,taxonName,parentTaxId,depth,rank,sep="*")
        m<-rbind(m,unlist(strsplit(str,split="\\*")))
        previousTaxId <- taxId
      }
      else if (length(which(m[,2]==taxonName & m[,5]==rank)) > 0) {
        row <- which(m[,2]==taxonName & m[,5]==rank)
        #error check -> is parent tax id same for both or not?
        #if not, remove the sequence
        if (parentTaxId!=m[row,3]) {
          #x<-x[-i]
          badHierarchy<-c(badHierarchy,i)
        }
        previousTaxId <- m[row,1]
      }
      #end seach
    }
  }

  out<-apply(m, MARGIN=1, FUN=function(x) paste(x,collapse="*"))
  write(out, file=file.path(dir, "train.txt"))
  #create parsed training files
  if(length(badHierarchy)>0) {
    warning("Following sequences had bad sequence hierarchy information, so removing: ",
      names(x)[badHierarchy],"\n")
    x<-x[-badHierarchy]
  }
  writeXStringSet(x,file.path(dir,"train.fasta"))

  if(system(paste(.javaExecutable(), java_args, "-jar", .get_rdp(), "train",
    "-o", dir, "-t", file.path(dir,"train.txt"),
    "-s", file.path(dir, "train.fasta")),
    ignore.stderr=TRUE, ignore.stdout=TRUE, intern=FALSE)){
    unlink(dir, recursive=TRUE)
    stop("Creating the classifier failed!")
  }

  file.copy(system.file("java/rRNAClassifier.properties",
    package="rRDP"), dir)

  rdp(dir)
}

removeRDP <- function(object) {
  ### first check if it looks like a RDP directory!
  if(!.isRDP(object$dir)) stop("The given RDPClassifier/directory does not look valid! Please remove the directory manually!")

  ### don't remove the default data
  if(object$dir != normalizePath(system.file("extdata","16srrna", package="rRDPData")))
    unlink(object$dir, recursive=TRUE)
}
