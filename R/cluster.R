#' @title Cluster data
#' 
#' @import pheatmap
#' @export
cluster <- function(obj, nClust=2, silent=T, hcl_method="all", dist_method="all") {
	if (hcl_method == "all") {
	    hcl_method <- c("ward.D", "single", "complete", "average", "mcquitty","median", "centroid", "ward.D2")
	}
	if (dist_method == "all") {
	    dist_method <- c("euclidean", "maximum", "manhattan", "canberra","binary", "minkowski")
	}


	res <- obj@cluster
	
	for (i in 1:length(obj@signatures)) {
		tmp <- obj@expr[which(rownames(obj@expr) %in% obj@signatures[[i]]$data$Gene),]
		print(dim(tmp))
		#tmp <- tmp[complete.cases(tmp*0),,drop=F]
		if (obj@clusterZ == 1) {
			tmp <- t(scale(t(tmp)))
		}	
		tmp[which(is.infinite(tmp))] <- 0
		tmp[which(is.na(tmp))] <- 0
		
		for (mh in hcl_method) {
		    for (md in dist_method) {
			pm <- pheatmap(tmp, silent=silent, clustering_method=mh, 
				       clustering_distance_cols=md)
			pm.clust <- cutree(pm$tree_col, nClust)
		    
			res[[length(res)+1]] <- list(ID=obj@signatures[[i]]$ID,
						     nClust=nClust, 
						     clust=pm.clust,
						     hclust_method=mh,
						     dist_method=md)
		    }
		}
	}	

	obj@cluster <- res

	
	return(obj)
}

#' @title Plot cluster assignment
#' @import pheatmap
#' @export 
plotCluster <- function(obj, cm="ward.D2") {
	tmp <- obj@cluster 
	dat <- lapply(tmp, function(x) x$clust)
	ids <- lapply(tmp, function(x) x$ID)


	dat <- do.call(rbind, dat)
	rownames(dat) <- paste("RN", 1:length(dat[,1]))

	anR <- data.frame(SIGN=unlist(ids))
	rownames(anR) <- rownames(dat)

	pheatmap(dat, clustering_method=cm, annotation_row=anR, cluster_rows=F)

	return(list(data=dat))
}
