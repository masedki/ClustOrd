setClass(
  Class = "ClustOrdpartitions", 
  representation = representation(zMAP="numeric",  tik="matrix"), 
  prototype = prototype(zMAP=numeric(), tik=matrix(0,0,0))
)

setClass(
  Class = "ClustOrdstrategy", 
  representation = representation(iterMH="numeric", nbSmall="numeric", iterSmall="numeric", nbKeep="numeric", iterKeep="numeric", tolKeep="numeric"), 
  prototype = prototype(iterMH=numeric(), nbSmall=numeric(), iterSmall=numeric(), nbKeep=numeric(), iterKeep=numeric(), tolKeep=numeric())
) 

## Constructeur de la classe S4 ClustOrdstrategy
ClustOrdstrategy <- function(iterMH, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep){
  if( nbKeep > nbSmall)
    nbKeep <- nbSmall
  new("ClustOrdstrategy", iterMH=iterMH, nbSmall=nbSmall, iterSmall=iterSmall, nbKeep=nbKeep, iterKeep=iterKeep, tolKeep=tolKeep)
}


setClass(
  Class = "ClustOrdmodel", 
  representation = representation(g="numeric", omega="numeric"), 
  prototype = prototype(g=numeric(), omega=numeric())
)

setClass(
  Class = "ClustOrddetailsMH", 
  representation = representation(Bestmodel="matrix", Currentmodel="matrix", Candidatemodel="matrix", allbic="matrix"), 
  prototype = prototype(Bestmodel=matrix(), Bestmodel=matrix(), Bestmodel=matrix())
)


setClass(Class = "ClustOrdcriteria", 
         representation = representation(loglikelihood="numeric", BIC="numeric", ICL="numeric", nbparam="numeric"), 
         prototype = prototype(loglikelihood=numeric(), BIC=numeric(), ICL=numeric(), nbparam=numeric())
)

setClass(
  Class = "ClustOrdparam", 
  representation = representation(pi="numeric", alpha="list", beta="list", epsilon="matrix"), 
  prototype = prototype(pi=numeric(), alpha=list(), beta=list(), epsilon=matrix())
)


setClass(
  Class = "ClustOrddata", 
  representation = representation(
    n="numeric",
    d="numeric",
    data="matrix",
    modalities="numeric"
  ), 
  prototype = prototype(
    n=numeric(),
    d=numeric(),
    data=matrix(),
    modalities=numeric()
  )
)

setClass(
  Class = "ClustOrdresults", 
  representation = representation(data="ClustOrddata", criteria="ClustOrdcriteria", partitions="ClustOrdpartitions",
                                  model="ClustOrdmodel", strategy="ClustOrdstrategy", param="ClustOrdparam", detailsMH="ClustOrddetailsMH"), 
  prototype = prototype(data=new("ClustOrddata"), criteria=new("ClustOrdcriteria"), partitions=new("ClustOrdpartitions"),
                        model=new("ClustOrdmodel"), strategy=new("ClustOrdstrategy"), param=new("ClustOrdparam"), detailsMH=new("ClustOrddetailsMH"))
)

BuildS4object <- function(x, g, iterMH, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep){
  data <- new("ClustOrddata", n=nrow(x), d=ncol(x), data=x, modalities=max(x)+1)
  strategy <- ClustOrdstrategy(iterMH, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)
  output <- new("ClustOrdresults", data=data, 
                strategy=strategy, 
                model=new("ClustOrdmodel", g=g, omega=Cleanmodel(sample(1:(ncol(x)), ncol(x), replace = T))[2,]-1),
                criteria=new("ClustOrdcriteria", BIC=-Inf, ICL=-Inf, loglikelihood=-Inf, nbparam=0))
  return(output)
}
