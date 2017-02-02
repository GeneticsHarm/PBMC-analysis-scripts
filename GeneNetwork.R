##
## Object for the creation of co-expression network graph
##
## Author: Harm Brugge
##
## Parameters:
##    data: Gene expression matrix where rows are samples and columns are genes
##    start.gene: The starting point of the network
##    min.correlation: Cut-off value for an edge between genes
##
##

GeneNetwork <- setClass(
  "GeneNetwork",
  
  slots = c(
    data = "matrix",
    start.gene = "character",
    network = "matrix"
  ),
  
  validity = function(object) {
    if(ncol(object@data) != nrow(object@data)) {
      return("Matrix should be symetrical")
    }
    return(TRUE)
  })

setGeneric("GeneNetwork", function(data, start.gene, min.correlation) standardGeneric("GeneNetwork")) 

setMethod("GeneNetwork", signature(data="matrix", start.gene="character", min.correlation="numeric"), 
          definition = function(data, start.gene, min.correlation) {
            
            object <- new("GeneNetwork", data=data, start.gene=start.gene)
            
            #is.correlated <- (object@data[start.gene,] > min.correlation | object@data[start.gene,] < -min.correlation)
            is.correlated <- order(abs(object@data[start.gene,]), decreasing = T)[1:5]
            
            
            correlated.genes <- object@data[start.gene, is.correlated]
            init.length <- length(correlated.genes)
            
            if (init.length <= 1) {
                stop("No correlated genes found")
            }
            
            object@network <- matrix(0L, nrow = init.length, ncol = init.length, dimnames = list(names(correlated.genes), names(correlated.genes)))
            object@network[start.gene,] <- 1
            object@network[,start.gene] <- 1
            
            for (gene in names(correlated.genes)) {
              object <- addGene(object, gene, min.correlation)
            }
            
            return(object)
          }
) 

setGeneric("addGene",  function(object, gene.name, min.correlation) standardGeneric("addGene"))
  
setMethod(f="addGene", signature("GeneNetwork", "character", "numeric"),
          definition = function(object, gene.name, min.correlation) {
            
              ## 
              #is.correlated <- (object@data[gene.name,] > min.correlation | object@data[gene.name,] < -min.correlation)
              
              # take top 5 correlated genes
              is.correlated <- order(abs(object@data[gene.name,]), decreasing = T)[1:5]
              correlation.vector <- object@data[gene.name, is.correlated]
            
              if(!gene.name %in% rownames(object@network)) {
                object@network <- rbind(object@network, 0L)
                object@network <- cbind(object@network, 0L)
                matrix.size <- dim(object@network)[1]
                rownames(object@network)[matrix.size] <- gene.name
                colnames(object@network)[matrix.size] <- gene.name
              }
              
              for (gene in names(correlation.vector)) {
                if (!gene %in% rownames(object@network)) {
                  object@network <- rbind(object@network, 0L)
                  row.count <- dim(object@network)[1]
                  rownames(object@network)[row.count] <- gene
                  object@network[row.count, gene.name] <- 1
                  
                  object@network <- cbind(object@network, 0L)
                  colnames(object@network)[row.count] <- gene
                  object@network[gene.name, row.count] <- 1
                } else {
                  object@network[gene.name, gene] <- 1
                  object@network[gene, gene.name] <- 1
                }
              }
              
              return(object)
            })



