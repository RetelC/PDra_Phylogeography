################################################################
### Functions created during phylogeographic analysis of     ###
### Pterosagitta draco along a basin-scale Atlantic transect ###
### Cas Retel MSc                                            ###
### casretel@gmail.com                                       ###
### Meant for personal use                                   ###
################################################################

# With the various R packages all using different formats
# to store sequences, it is worthwhile to archive how to extract
# them (and their names) from every object;
# ape::DNAbin: as.character(seqs); labels(seqs)
#     t(sapply(seqs, function(x) x))
# seqinr::seqFastadna: getSequence(seqs); labels(seqs)

pasteZeroes <- function(x){
  # Returns a character vector with zeroes pasted at the start
  # of elements of x, such that every element is the same length
  # Created to make ordering of 1, 2, 3, 10 correct; 
  # as.factor() would order them as c(1, 10, 2, 3), 
  # creating inconsistencies when using ggplot()
  x <- as.character(x)
  nmax <- max(sapply(x, nchar))
  while(any(sapply(x, nchar) != nmax)){
    x <- paste("0", x, sep="")
    x[sapply(x, nchar)>nmax] <- 
      substr(x[sapply(x, nchar)>nmax], start=2, stop=nmax+1)
  }
  return(x)
}


# pDist <- function(x1, x2){
#   # PDist() returns uncorrected divergence of two 
#   # nucleotide character vectors: nr of differences per site
#   # It is assumed that sequences are aligned, and 
#   # "N" values are excluded from calculations
#   informative <- (x1!="N" & x2!="N")
#   return(sum(x1[informative]!=x2[informative])/sum(informative))
# }
# Obsolete: use ape::dist.dna()

# JCDist <- function(x1, x2){
#   # JCDist() returns Jukes-Cantor corrected divergence of two 
#   # nucleotide character vectors
#   # It is assumed that sequences are aligned, and 
#   # "N" values are excluded from calculations
#   informative <- (x1!="N" & x2!="N")
#   pi <- sum(x1[informative]!=x2[informative])/sum(informative)
#   return (-3*log(1-(4*pi/3))/4)
# }
# Obsolete: use ape::dist.dna()

# 
# TNDist <- function(x1, x2){
#   # TNDist() returns Tajima-Nei corrected divergence of two 
#   # nucleotide character vectors. Similar to J-C, T-N recognized
#   # that unequal nucleotide frequencies will result in an 
#   # overestimation by this calculation
#   # It is assumed that sequences are aligned, and 
#   # "N" values are excluded from calculations
#   nucs <- c("A", "T", "G", "C")
#   informative <- (x1!="N" & x2!="N")
#   x1 <- x1[informative]
#   x2 <- x2[informative]
#   n <- sum(informative)
#   
#   pi <- sum(x1!=x2)/n
#   x <- sapply(nucs, function(a1)
#     sapply(nucs, function(a2)
#       sum((x1%in%a1 & x2%in%a1) + 
#             (x1%in%a2 & x2%in%a2))/n))
#   q <- sapply(nucs, function(x) sum( c(x1, x2)==x)/(2*n))
#   
#   h <- sum(unlist(lapply(1:3, function(i)
#     sapply(i:4, function(j) x[i, j]^2 / (q[i]*q[j]) ))))/2
#   b <- (1-sum(q^2) + (pi^2)/h)/2
#   return(-1*b*log(1-(pi/b)))
# }
# Obsolete: use ape::dist.dna()

pDist2 <- function(x1, x2){
  # PDist2() returns uncorrected divergence of two 
  # nucleotide character vectors: nr of differences per site
  # Different from PDist(), it allows not only "ACTGN-", but also
  # "DHKMRSTWY"
  # It is assumed that sequences are aligned, and 
  # "N" values are excluded from calculations
  require("seqinr")
  informative <- which(x1!="N" & x2!="N" & x1!="-" & x2!="-")
  gaps <- rbind((x1=="-"), (x2=="-"))
  unequal <- sum(!sapply(informative, function(i) 
    any(amb(x1[i])%in%amb(x2[i])))) + sum(colSums(gaps)==1)
  return(unequal/(length(informative)+sum(colSums(gaps)>0)))
}
# Obsolete: use ape::dist.dna()


distClades <- function(seq, cladeseqs, cladefact){
  # distClades was created to assign new sequences clades, 
  # based on an available sequence set that already have clades
  # seq = matrix of (unassigned sequence) alignment
  # cladeseqs = matrix of sequences that already were assigned a clade
  # cladefact = factor giving clades corresponding to cladeseqs
  # returns a data frame with for every new sequence (per row)
  # the average uncorrected divergence to all sequences per clade (columns)
  if(!is.matrix(seq)) seq <- as.matrix(seq)
  if(!is.matrix(cladeseqs)) cladeseqs <- as.matrix(cladeseqs)
  cladefact <- as.factor(cladefact)
  n <- nrow(seq)
  nc <- length(levels(cladefact))
  out <- matrix(0, nrow=n, ncol=nc)
  
  for(i in 1:n){
    for(cla in 1:nc){
      claind <- which(cladefact==levels(cladefact)[cla])
      out[i, cla] <- mean(sapply(claind, function(j)
        dist.dna(seq[i, ], cladeseqs[j, ], model="raw")))
    }
  }
  return(out)
}

revComp <- function(x){
  # Returns the reverse complement of a sequence, 
  # of either string or character vector format
  library(seqinr)
  if(length(x) == 1){
    return(c2s(comp(rev(s2c(x)), force=F)))
  }else{
    return(comp(rev(x), force=F))
  }
}

seqsToLocus <- function(seqs, mname="Marker1"){
  # from an aligned sequence matrix, removes invariable positions
  # and returns a matrix of ordinal locus values
  if(class(seqs) == "DNAbin") seqs <- as.character(seqs)
  out <- data.frame(V1=rep(0, nrow(seqs)))
  hnum <- 1
  while(any(out$V1==0)){
    uncl <- which(out==0)
    out$V1[uncl[1]] <- hnum
    for(i in uncl[-1]){
      if(all.equal(seqs[uncl[1], ], seqs[uncl[i], ])==TRUE){
        out$V1[i] <- hnum
      }
    }
    hnum <- hnum+1
  }
  colnames(out) <- mname
  return(out)
}

seqsToDataframe <- function(seqs, ...){
  # from a sequence matrix (with individuals as rows), 
  # returns a data frame consisting of columns: 
  # ...-arguments as factors, followed by all polymorphic positions
  # ! only works for my own IonTorrent data, can't handle ambiguous
  # symbols etc. !
  if(class(seqs) == "DNAbin") seqs <- as.character(seqs)
  factors <- list(...)
  if(any(sapply(factors, length)!=nrow(seqs))){
    stop("Sequence and factor lengths do not match")
  }
  colnames(seqs) <- paste("pos", 1:ncol(seqs), sep="")
  loci <- seqs[, apply(seqs, 2, function(x) length(unique(x))!=1)]
  df.out <- data.frame(subject=rownames(loci))
  if(length(factors)>0){
    for(i in 1:length(factors)){
      df.out <- cbind(df.out, factors[[i]])
    }
  }
  colnames(df.out) <- c("subject", names(factors))
  df.out <- cbind(df.out, loci)
  return(df.out)
}

lengthdet <- function(...){
  # convenience function
  fa <- list(...)
  return(length(fa))
}

haploToNum <- function(x, diploid=F){
  # from a sequence matrix in seqsToDataframe-format, replaces 
  # nucleotide symbols with numeric values 1:4, to be compatible 
  # to hierfstat::test.xxx-functions. 
  # Compatible with seqsToDataframe(), and to this function only. 
  # if diploid=T, haploToNum expects a sequence matrix with
  # two sequences per individual, positioned after each other
  # Hence, rows 1, 2 are one individual, as are rows 205, 206. 
  # If homozygous, diplotype of two identical haplotypes is expected
  if(diploid & (nrow(x)%%2)) stop("Odd number of sequences")
  
  x <- as.matrix(x[, grep("pos", colnames(x))])
  symbollist <- x %>% as.character %>% unique
  out <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
  
  for(i in 1:length(symbollist)){
    out[x==symbollist[i]] <- i
  }
  if(diploid){
    out <- sapply(2*(1:(nrow(x)/2)), function(i)
      paste(out[i-1, ], out[i, ], sep="")) %>% t
  }
  out <- apply(out, 2, as.numeric)
  return(out)
}

nameToIndiv <- function(x, appendHetState=F){
  # Convenience function: From a vector of sequence labels of the form
  # Pdra_AMT22_<stationnr>_<indivnr>_<heterozygositytag>, 
  # returns "<stationnr>_<indivnr>"
  if(!appendHetState){
    x %>% strsplit(split="AMT22_") %>%
      (function(x) sapply(x, "[[", 2)) %>%
      substr(start=1, stop=5)
  }else{
    x %>% strsplit(split="AMT22_") %>%
      (function(x) sapply(x, "[[", 2)) %>%
      substr(start=1, stop=7)
  }
}
nameToStation <- function(x){
  # Convenience function: From a vector of sequence labels of the form
  # Pdra_AMT22_<stationnr>_<indivnr>_<heterozygositytag>, 
  # returns "<stationnr>"
  x %>% strsplit(split="AMT22_") %>%
    (function(x) sapply(x, "[[", 2)) %>%
    substr(start=1, stop=2) %>% as.factor
}
nameToBiome <- function(x){
  # Convenience function: From a vector of sequence labels of the form
  # Pdra_AMT22_<stationnr>_<indivnr>_<heterozygositytag>, 
  # returns biome
  # ! Only works for sequences sampled at thirteen stations used in 
  # P. draco research !
  out <- x %>% strsplit(split="AMT22_") %>%
    (function(x) sapply(x, "[[", 2)) %>%
    substr(start=1, stop=2) %>% 
    (function(x) as.factor(1 + (x>22) + (x>32) + (x>48) + (x>61)))
  levels(out) <- c("N temp", "N gyre", "Equat", "S gyre", "S temp")
  return(out)
}

  
