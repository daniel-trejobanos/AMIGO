BROptions <- setClass(
  "BROptions",
  slots = c(
    chainLength = "numeric",
    burnin = "numeric",
    seed   = "numeric",
    numThread = "numeric",
    numThreadSpawned = "numeric",
    preProcessChunks = "numeric",
    thin = "numeric",
    S= "numeric",
    numGroups = "numeric",
    mS = "numeric",
    groupFile = "character",
    title = "character",
    analysisType = "character",
    bayesType = "character",
    phenotypeFile = "character",
    bedFile = "character",
    mcmcSampleFile = "character",
    optionFile = "character",
    compress = "logical",
    dataType = "numeric"
    ),
  prototype=list(
    chainLength = 10000,
    burnin = 5000,
    seed =   as.numeric(Sys.time()),
    numThread = 1,
    numThreadSpawned = 0,
    preProcessChunks = 1,
    thin = 5,
    S = c(0.01,0.001,0.0001),
    title="brr",
    analysisType = "bayes",
    bayesType= "C",
    phenotypeFile = "",
    bedFile = "",
    mcmcSampleFile = "bayesOutput.csv",
    optionFile = "",
    numGroups =2,
    dataType =1
  )
  
) 


