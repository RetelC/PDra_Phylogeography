#################################################################
### Functions created during phylogeographic analysis of      ###
### Pterosagitta draco along a basin-scale Atlantic transect, ###
### related to processing of IonTorrent NGS output            ###
### Cas Retel MSc                                             ###
### casretel@gmail.com                                        ###
#################################################################

count.clstr <- function(file){
  # count.clstr reads in a .clstr file and returns a data frame with 
  # column 1 the cluster names present and 
  # column 2 the number of reads mapping to every cluster, 
  # ordered by column 2
  # file=path to .clstr file, after Galaxy quality filtering and clustering
  # Reference: Retel 2015 - Deep mitochondrial divergence does not imply 
  #            reproductive isolation
  lines <- readLines(file)
  ind <- which(substr(lines, 1L, 1L) == ">")
  nclust <- length(ind)
  nseq <- length(lines)-nclust
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }
  start <- ind + 1
  end <- c((ind - 1)[-1], length(lines))
  nperclust <- end-(start-1)
  return(data.frame(read=lines[ind], nperclust=nperclust, 
                    prop.reads=nperclust/sum(nperclust))[
                      order(nperclust, decreasing=T), ])
}

select.clstr <- function(file, percentage=0.3){
  # select.clstr reads in a .clstr file and returns a data frame containing 
  # all clusters to which more than (percentage) of the reads map to. 
  # file = path to .clstr file, after Galaxy quality filtering and clustering
  # percentage = proportion of total reads that should map to a cluster, 
  #              in order to select the concerning cluster
  # Reference: Retel 2015 - Deep mitochondrial divergence does not imply 
  #            reproductive isolation
  lines <- readLines(file)
  ind <- which(substr(lines, 1L, 1L) == ">")
  nclust <- length(ind)
  nseq <- length(lines)-nclust
  
  cutoff <- nseq*percentage
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }
  
  start <- ind + 1
  end <- c((ind - 1)[-1], length(lines))
  nperclust <- end-(start-1)
  ind.cutoff <- which(nperclust >= cutoff)
  return(data.frame(read=lines[ind], nperclust=nperclust, 
                    prop.reads=nperclust/sum(nperclust))[ind.cutoff, ])
}

select2.clstr <- function(file, nclust.ret=2){
  # select2.clstr reads in a .clstr file and returns a data frame containing 
  # the nclust.ret clusters to which the most reads map
  # file = path to .clstr file, after Galaxy quality filtering and clustering
  # nclust.ret = number of clusters to return
  # Reference: Retel 2015 - Deep mitochondrial divergence does not imply 
  #            reproductive isolation  
  lines <- readLines(file)
  ind <- which(substr(lines, 1L, 1L) == ">")
  nclust <- length(ind)
  nseq <- length(lines)-nclust
  
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }
  
  start <- ind + 1
  end <- c((ind - 1)[-1], length(lines))
  nperclust <- end-(start-1)
  ind.ret <- order(nperclust, decreasing=T)[1:nclust.ret]
  return(data.frame(read=lines[ind], nperclust=nperclust, 
                    prop.reads=nperclust/sum(nperclust))[ind.ret, ])
  # return(data.frame(read=lines[ind][ind.ret], nperclust=nperclust[ind.ret]))
}

readselect.fasta <- function(file, index=data.frame(read=">Cluster_0", nperclust=1)){
  # readselect.fasta was meant to be used in conjunction with either 
  # select.clstr() or selec2.clstr() - functions. 
  # Reads in sequences of the selected clusters, which can then be 
  # written to a .fasta file using seqinr::write.fasta() or a similar function
  # file = path to .fasta file containing cluster sequences 
  #        (named Cluster_0, Cluster_1, etc)
  # index = data.frame, directly compatible with select.clstr
  # ! code is adapted read.fasta() from seqinr package !
  lines <- readLines(file)
  ind <- which(substr(lines, 1L, 1L) == ">")
  nseq <- length(ind)
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }
  if(nrow(index)==0){
    out <- "N"; names(out) <- "no_majority"
    return(out)
  }
  start <- ind + 1
  end <- c((ind - 1)[-1], length(lines))
  
  # seqinr used lapply: I think using sapply doesn't screw it up
  sequences <- sapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], 
                                                       collapse = ""))
  nomseq <- sapply(seq_len(nseq), function(i) {
    firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })
  names(sequences) <- nomseq
  
  clusters <- as.character(index[, "read"])
  clusters <- substr(clusters, start=2, stop=nchar(clusters))
  index.selected <- match(clusters, table=nomseq)
  return(sequences[index.selected])
}
