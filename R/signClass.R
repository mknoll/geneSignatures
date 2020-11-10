# S4 Helper classes to allow NULL values
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("dfOrNULL", c("data.frame", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))


#' A S4 class representing an analysis instance 
#'
#' @name signClass-class
#' @exportClass signClass
incClass <- setClass("signClass",
                     slots=c(dateCreated="POSIXct",
			    expr="dfOrNULL",
			    signatures="listOrNULL",
			    cluster="listOrNULL",
			    clusterZ="numeric",
			    fraction="listOrNULL",
			    text="characterOrNULL"
			    ))


#' @title Constructor incClass
#' @description Constructor for new analysis
setMethod("initialize", "signClass",
          function(.Object,
		   expr=data.frame,
                   ...) {
	      	.Object@dateCreated=Sys.time()
		.Object@expr=expr
		.Object@clusterZ=1

		.Object
	  })

