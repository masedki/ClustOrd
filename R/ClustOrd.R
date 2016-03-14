clustOrd <- function(x, val, path, iterMH=200, nbSmall=25, iterSmall=5, nbKeep=10, iterKeep=50, tolKeep=0.01, init=NULL){  
  g <- val[1]
  number <- val[2]
  reference <- BuildS4object(x, g, iterMH, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)
  reference <- OptimizeBIC(reference)
  gc(verbose = FALSE)
  cat("K", g, " Chain",number, " BestICL", round(reference@criteria@ICL), " Model", t(reference@model@omega),"\n")
  if (missing(path)==FALSE)
    save(reference, file = paste(path,"/nbclasses",g,"/number",number,".rda", sep=""))
  return(reference)
}

Cleanmodel <- function(vec){
  passage <- matrix(NA, 2,length(vec))
  passage[1,] <- vec
  passage[2,1] <- 1
  locmax <- 1
  for (u in 2: length(vec)){
    if (passage[1,u] %in% passage[1,1:(u-1)]){
      passage[2,u] <- passage[2, which(passage[1,1:(u-1)] == passage[1,u])[1]]
    }else{
      locmax <- locmax + 1
      passage[2,u] <- locmax
    }
  }
  return(passage)
}

CleanOutput <- function(reference){
  reference@model@omega <- as.numeric(reference@model@omega) + 1 
  g <- reference@model@g
  reference@param@pi <- as.numeric(reference@param@pi)
  reference@param@epsilon <-   reference@param@epsilon[, 1:max(reference@model@omega)]
  rownames(reference@param@epsilon) <- paste("Class", 1:g, sep="-")
  colnames(reference@param@epsilon) <- paste("Block", 1:max(reference@model@omega), sep="-")
  length(reference@param@alpha) <- max(reference@model@omega)
  length(reference@param@beta) <- max(reference@model@omega)
  names(reference@param@alpha) <- paste("Block", 1:max(reference@model@omega), sep="-")
  names(reference@param@beta) <- paste("Block", 1:max(reference@model@omega), sep="-")
  for (b in 1:length(reference@param@alpha)){
    names(reference@param@alpha[[b]]) <- colnames(reference@data@data)[which(reference@model@omega==b)]
    for (j in 1:length(reference@param@alpha[[b]])){
      reference@param@alpha[[b]][[j]] <- matrix(reference@param@alpha[[b]][[j]] ,nrow = g)
      colnames(reference@param@alpha[[b]][[j]]) <- as.character(0:(reference@data@modalities-1))
      rownames(reference@param@alpha[[b]][[j]]) <- paste("Class",1:g,sep="-")
    } 
    if (sum(reference@model@omega == b)>1){
      reference@param@beta[[b]]<- matrix(reference@param@beta[[b]] ,nrow = g)
      colnames(reference@param@beta[[b]]) <- as.character(0:(reference@data@modalities-1))
      rownames(reference@param@beta[[b]]) <- paste("Class",1:g,sep="-") 
    }
  }
  
  
  #   interm <- reference
  #   passage <- matrix(NA,2, reference@data@d)
  # 
  #   passage <- Cleanmodel(reference@model@omega)
  #   reference@model@omega <- passage[2,]
  #   
  #   reference@param@epsilon <- as.matrix(interm@param@epsilon[,unique(passage[1,])])
  #   colnames(reference@param@epsilon) <- paste("Block", 1:ncol(reference@param@epsilon), sep="-")
  #   rownames(reference@param@epsilon) <- paste("Class", 1:g, sep="-")
  #   
  #   reference@param@alpha <- list()
  #   reference@param@beta <- list()
  #   for (b in 1:ncol(reference@param@epsilon)){
  #     avt <- passage[1,which(passage[2,]==b)[1]]
  #     reference@param@alpha[[b]] <- interm@param@alpha[[avt]]
  #     if (sum(passage[2,]==b)>1){
  #       reference@param@beta[[b]] <- interm@param@beta[[avt]]
  #       names(reference@param@beta)[[b]] <- paste("Block",b,sep="-")
  #     }else{
  #       reference@param@beta[[b]] <- NULL
  #     }
  #   }
  #   names(reference@param@alpha) <- paste("Block", 1:ncol(reference@param@epsilon), sep="-")
  # 
  # 
  #   reference@detailsMH@Bestmodel <- t(apply(reference@detailsMH@Bestmodel, 1, f <- function(vec){return(Cleanmodel(vec)[2,])}))
  #   reference@detailsMH@Currentmodel <- t(apply(reference@detailsMH@Currentmodel, 1, f <- function(vec){return(Cleanmodel(vec)[2,])}))
  #   reference@detailsMH@Candidatemodel <- t(apply(reference@detailsMH@Candidatemodel, 1, f <- function(vec){return(Cleanmodel(vec)[2,])}))
  #   
  #   colnames(reference@detailsMH@allbic) <- c("Best", "Current", "Candidate")
  
  return(reference)
}