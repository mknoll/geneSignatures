#' @title Calculate fraction of similarly altered genes
#'
#' @description For each gene, over- or underexpression is
#' determined by comparing each value to the row-wise values
#' calculated via fun (default: median; median of all samples).
#' The number of similiarly altered genes compared to the gene
#' signature is calculated. Fractions (for each signature)
#' and binary classifications, based on the fractCut and CI
#' parameters are returned.
#'
#' @param obj singGene instance
#' @param fun Function to determine up/downregulation based
#' on all values (all samples) for each gene
#' @param fractCut Fraction cutoffs to classify sample as
#' not signature associated or signature associated. If
#' CI is set to 0, samples with fractions below the first
#' elements of the fractCut vector are set to 0, above the
#' second element of fractCut, set to 1. Other samples are
#' classified as NA.
#' @param CI Calculate Wilson confidence interval for fractions
#' and compare upper 95% limits for the lower cutoff (first
#' element of fractCut) and lower 95% limits for upper cutoffs
#' (second element of fractCut)
#'
#' @import Hmisc
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export
#'
#' @return Updated singGene object
fract <- function(obj, fun=median, fractCut=c(0.5, 0.5), CI=0) {
    # TODO: check if necessarydata is available
    res <- obj@fraction

    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)

    out <- foreach(i=1:length(obj@signatures)) %dopar% {
	tmp <- obj@expr[which(rownames(obj@expr) %in% obj@signatures[[i]]$data$Gene),
			, drop = F]
	if (length(tmp[,1]) >= 1) {
	    sg <- obj@signatures[[i]]$data
	    sg <- sg[match(rownames(tmp), sg$Gene), ,drop=F]
	    cut <- apply(tmp, 1, fun)
	    alt <- apply(tmp, 2, function(x) ifelse(x > cut, "up",
						    "down"))
	    if (class(alt) == "character") {
		alt <- data.frame(t(alt))
		rownames(alt) <- names(cut)
	    }
	    mtch <- apply(alt, 2, function(x) length(which(x == sg[,2])))
	    nMx <- length(unique(sg$Gene))

	    #Wilson confidence interval
	    #int <- binconf(mtch, nMx, method="Wilson")

	    if (CI == 1) {
		int <- binconf(mtch, nMx, method="Wilson")
		grp <- ifelse(int[,"Lower"] > fractCut[2], "UP",
			      ifelse(int[,"Upper"] <= fractCut[1], "DOWN", NA))
	    } else {
		int <- data.frame(FREQ=mtch/nMx)
		grp <- ifelse(int[,1] > fractCut[2], "UP",
			      ifelse(int[,1] <= fractCut[1], "DOWN", NA))
	    }

	    list(ID=obj@signatures[[i]]$ID,
		 #ret <- list(ID=obj@signatures[[i]]$ID,
		 mtch=mtch/nMx,
		 CI=int,
		 grp=grp)
	} else {
	    list(ID=obj@signatures[[i]]$ID,
		 #ret <- list(ID=obj@signatures[[i]]$ID,
		 mtch=NA,
		 CI=NA,
		 grp=NA)
	}
    }
    res <- c(res, out)

    obj@fraction <- res

    doParallel::stopImplicitCluster()

    return(obj)
}


#' @title Plot calculated fractions of singature
#' overlaps
#'
#' @description Plonts continuous or categorical
#' classifications of samples based on signatures
#' and returns the respective matrices
#'
#' @param obj singGene object
#' @param type plot continous fractions (cont)
#' or classified sample matrix (cat)
#'
#' @export
#'
#' @return list with elements GRP (classified) and
#' N (continuous fraction)
plotFract <- function(obj, type="cont") {
    tmp <- obj@fraction
    ids <- lapply(tmp, function(x) x$ID)
    print(ids)

    ## numbers of matches
    mt <- lapply(tmp, function(x) x$mtch)
    mt <- do.call(rbind, mt)
    rownames(mt) <- ids
    if (type == "cont") {
	tryCatch({
	    pheatmap(mt[complete.cases(mt),])
	},error=function(e) { print(e) })
    }

    ## Cat
    m <- lapply(tmp, function(x) x$grp)
    m <- do.call(rbind, m)
    rownames(m) <- ids
    m <- data.matrix(m)
    m[which(m == "UP")] <- 1
    m[which(m == "DOWN")] <- 0

    m <- data.matrix(apply(m, 2, function(x) as.numeric(as.character(x))))
    rownames(m) <- ids
    if (type != "cont") {
	tryCatch({
	    pheatmap(m[complete.cases(m),])
	}, error=function(e) { print(e) })
    }

    return(list(GRP=m, N=mt))
}


