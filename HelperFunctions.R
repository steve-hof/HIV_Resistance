#### helper functions

# function to check that the mutations have been entered correctly
# if the letters are entered as lower case, they are converted to upper case

# convert insertions or deletions to # or ~
convert.muts <- function(muts.in){
  muts.in5 <- which(nchar(muts.in) == 5)
  muts.in6 <- which(nchar(muts.in) == 6)
  for(mi in muts.in5){
    postmp <- substr(muts.in[mi],1,2)
    if(substr(muts.in[mi],3,5) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],3,5) == "del") muts.in[mi] <- paste0(postmp,"~")   
  }
  for(mi in muts.in6){
    postmp <- substr(muts.in[mi],1,3)
    if(substr(muts.in[mi],4,6) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],4,6) == "del") muts.in[mi] <- paste0(postmp,"~")   
  }
  return(muts.in)
}

check.muts <- function(muts.in){
  # all entries should be nchar=3
  if(any(nchar(muts.in) < 3) | any(nchar(muts.in) > 4))
    stop('All entries in argument "muts.in" should be between 3 and 4 characters long.')
  
  muts.in3 <- muts.in[nchar(muts.in) == 3]
  muts.in4 <- muts.in[nchar(muts.in) == 4]
  
  # all should have numbers first two or three characters and a letter for the last 
  if(!all(toupper(substr(muts.in3,3,3))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last 
         character.')
  if(any(is.na(as.numeric(substr(muts.in3,1,2)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
  
  if(!all(toupper(substr(muts.in4,4,4))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last 
         character.')
  if(any(is.na(as.numeric(substr(muts.in4,1,3)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
}



# function to create the design matrix X with the input mutations/positions
buildX <- function(dat, mut, ps){
  X <- matrix(NA, nrow=nrow(dat), ncol=length(mut))

  # loop through all positions
  for(p in unique(ps)){
    p1 <- substr(dat[,p],1,1)  # first mutation at this position
    p2 <- substr(dat[,p],2,2)
    for(ind in which(ps==p)){
      X[,ind] <-  as.numeric(p1==as.character(mut[ind]) | 
                               p2==as.character(mut[ind]))  
    }
  }  
  colnames(X) <- paste0(ps,mut)  
  return(X)
}


#########################################
### help function, extract mutation name
# called by "mut.names()" function
extractMut <- function(mut.set)
{
  require(tidyr)	
  X <- separate(mut.set, V2, paste0("V", 2:10), sep=" ")
  mut.set.name <- c()
  for (i in 1:nrow(X)){
    for (j in 1:(ncol(X)-1)){
      if(nchar(as.character(X[i,1]))<2) {X[i,1] <- paste0(0,X[i,1])}
      if (!(is.na(X[i,j+1]) | X[i,j+1]=="")){
        mut.set.name <- append(mut.set.name, paste0(X[i,1],X[i,j+1]))
      }
    }
  }
  return(mut.set.name)
}

#### Extract mutation patterns
mut.name <- function(Mut)
{
  
  if (Mut=="TSM.NRTI"){
    MutPt <- read.table('https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/NP_TSM/NRTI',
                        header=F, sep="\t", comment.char="@",
                        stringsAsFactors=FALSE)
  }
  if (Mut=="Comp.NRTI"){
    MutPt <- read.table('https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/COMPLETE/NRTI',
                        header=F, sep="\t", comment.char="@",
                        stringsAsFactors=FALSE)
  }
  if (Mut=="Exp.NRTI") {
    MutPt <- read.table('https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/MUTATIONLISTS/IASUSA/NRTI',
                        header=F, sep="\t", comment.char="@",
                        stringsAsFactors=FALSE)
  }
  Pt_name <- extractMut(MutPt)
  return(Pt_name)
}


##################
# Load data      #
##################

# function to load data from HIV Drug Resistance Database https://hivdb.stanford.edu/pages/genopheno.dataset.html

# requires R package 'tidyr'


### Inputs:

# dataset specifies which dataset the function is to be run on
# must be either "PR", "NRTI", or "NNRTI"

# muts.in is a character vector of the mutations to be included as predictor variables 
# each entry in the vector is 3 characters long
  # a letter representing the mutation in the first character
  # followed by two numbers representing the position

# drug is a string indicating the drug to be used for the dependent variable
# must be equal to one of the following: "FPV", "ATV", "IDV", "LPV", "NFV", 
# "RTV", "SQV", "TPV", "DRV"

# min.muts is the minimum number of sequences that a mutation must appear in.
# If a mutation appears in too few sequences, it is removed from the model.



### output :
# the data matrix for analysis


load.data <- function(dataset="NRTI", drug="LPV", min.muts=10, muts.in=mut.comp)
{
  
  require(tidyr)
  
  # check that arguments are entered correctly
  if(!is.character(muts.in)){
    stop('The argument "muts.in" must be a character vector.')
  }
  if(!is.character(drug) | length(drug)!=1){
    stop('The argument "drug" must be a character vector of length 1.')
  }
  
  if(!is.character(dataset) | !dataset%in%c("PI","NRTI","NNRTI")){
    stop('The argument "dataset" must be a character vector equal to "PI", 
         "NRTI", or "NNRTI".')
  }
  
  if(drug=="3TC") drug <- "X3TC"
  
  muts.in.conv <- convert.muts(muts.in)
  check.muts(muts.in.conv)
  
  # get the amino acids and positions for the mutations to be included in the model
  mut <- ifelse(nchar(muts.in)==3,toupper(substr(muts.in,3,3)),
                toupper(substr(muts.in,4,4)))
  ps <- suppressWarnings(ifelse(nchar(muts.in)==3,as.numeric(substr(muts.in,1,2)),
                                as.numeric(substr(muts.in,1,3))))  
  
  ## automatically read in the data using url
  if(dataset=="PI"){ 
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
                      header=TRUE, sep="\t",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,10:108]
  }
  if(dataset=="NRTI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)    
    dat[dat=="."] <- NA
    posu <- dat[,8:247]  
  }
  if(dataset=="NNRTI"){ 
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)       
    dat[dat=="."] <- NA
    posu <- dat[,6:245]  
  }
  
  # construct design matrix for OLS
  X <- buildX(posu, mut, ps)
  
  # construct dependent variable
  drugcol <- which(names(dat)==drug)    
  Y <- as.numeric(dat[,drugcol])  # absolute measure
  Ylog10 <- log10(Y)
  df.log <- data.frame(Y=Ylog10, X=X)
  
  # remove all rows with missing values
  rem.rows <- unique(which(is.na(df.log),arr.ind=TRUE)[,1])
  df.log.cc <- df.log[-rem.rows,]  # complete case
  
  # remove mutations that are rare
  rare.muts <- which(colSums(df.log.cc[,-1])<min.muts)
  if(length(rare.muts)>0){
    message(paste0(muts.in[rare.muts],
                   " excluded from the model because it appears in fewer than ",
                   min.muts," sequences.\n"))
    df.log.cc <- df.log.cc[,-(rare.muts+1)]  
  }
  
  return(df.log.cc)
}
