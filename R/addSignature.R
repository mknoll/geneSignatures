#' @title Add signature
#' @export
addSignatureGroup <- function(obj, group, meta=NULL) {
    ## Custom supplied file?
    pck <- F
    if (is.null(meta)) {
	meta <- system.file("extdata", "meta.csv", package="geneSignatures")
	pck <- T
    }
    fl <- read.csv(meta)
    if (colnames(fl)[1] == "X") {
	fl <-fl[,-1]
    }
    fl <- fl[which(fl$GROUP == group),,drop=F]
    for (i in 1:length(fl[,1])) {
	if (pck) {
	    tmp <- read.csv(system.file("extdata", fl$FILE[i], package="geneSignatures"))
	} else {
	    tmp <- read.csv(fl$FILE[i])
	}
	obj <- addSignatureDF(obj, tmp, ID=as.character(fl$FILE[i]))
    }
    return(obj)
}

#' @title Add signature providid via data.frame
#'
#' @description Adds a gene signature, provided as data.frame
#' with two columns.
#'
#' @param obj geneSign instance
#' @param df data.frame with two columns (first column: Gene,
#' second columns: Direction of alterations: "up" or "down").
#' columnnames: Gene, Arbitrary_Name
#' @param ID Arbitrary ID of signature
#'
#' @export
#'
#' @return Updated geneSign object
addSignatureDF <- function(obj, df, ID) {
    print(ID)
    #TODO: sanity checks
    #check if genes from df are expressed in cells in obj@expr
    dforig<-df
    df <- df[which(df$Gene %in% rownames(obj@expr)),]
    if(isTRUE(all.equal(dforig$Gene[order(dforig$Gene)],df$Gene[order(df$Gene)]))==TRUE){
      print("all genes expressed")
    }else{
      rest<-outersect(dforig$Gene,df$Gene)
      print("genes not expressed in cells, therefore excluded:")
      print(rest)
    }

    lst <- obj@signatures

    ### Weight
    weight <- rep(1, length(df[,1]))
    if ("Weight" %in% colnames(df)) {
	## check if sum is 1
    }
    lst[[length(lst)+1]] <- list(ID=ID,
				 data=df,
				 WEIGHT=weight)
    obj@signatures <- lst

    return(obj)
}



#' @title Load Signature Metadata
#' @description Data frame with filenames and group assignments
#' to use for analysis
#' @export
loadSignatureMetadata <- function(obj, df=NULL) {
	return(obj)
}

