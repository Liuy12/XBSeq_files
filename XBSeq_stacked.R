library('DESeq2')
library('matrixStats')
library('pracma')
library('locfit')

setClass( "XBSeqDataSet", 
          contains = "DESeqDataSet",
          slots = c( 
            fitInfo = "environment",
            dispTable = "character",
            conditions = 'factor',
            dispEst = 'list')
)


setValidity( "XBSeqDataSet", function( object ) {
  if(!is.factor(object@conditions))
    return("conditions have to be factors")
  TRUE
})


setGeneric("fitInfo", function(object, name=NULL) standardGeneric("fitInfo"))


setMethod('fitInfo', signature(object = "XBSeqDataSet"), 
          function( object, name){
            if( length( ls( object@fitInfo ) ) == 0 )
              stop( "No fits available. Call 'estimateSCV' first." )
            if( length( ls( object@fitInfo ) ) > 1 && is.null(name) )
              stop( "More than one fitInfo object available. Specify by name. (See 'ls(XB@fitInfo)' for a list.)" )
            if( length( ls( object@fitInfo ) ) == 1 && is.null(name) )
              name = ls( object@fitInfo )[ 1 ]
            object@fitInfo[[ name]]
          }
)


setMethod("conditions", signature(object="XBSeqDataSet"),
          function( object, ... ) {
            if(length(list(...))!=0)
              warning("in conditions: Ignoring second and/or further arguments.")
            conds <- object@conditions
            names( conds ) <- colnames( counts(object,1) )
            conds
          })   


setReplaceMethod("conditions", signature(object="XBSeqDataSet"),
                 function( object, value ) {
                   object@conditions <- factor( value )
                   validObject( object )
                   object
                 })


setMethod("dispTable", signature(object="XBSeqDataSet"),
          function( object ) {
            object@dispTable
          })   


setReplaceMethod("dispTable", signature(object="XBSeqDataSet"),
                 function( object, value ) {
                   object@dispTable <- value
                   validObject( object )
                   object
                 })   


setGeneric("dispEst", function(object, varname = NA) standardGeneric("dispEst"))


setGeneric("dispEst<-", function(object, varname = NA, value) standardGeneric("dispEst<-"))


setMethod("dispEst", signature(object="XBSeqDataSet"),
          function(object, varname) {
            if(!is.na(varname))
              object@dispEst[[varname]]
            else
              object@dispEst
          })


setReplaceMethod("dispEst", signature(object="XBSeqDataSet"),
                 function(object, varname,  value) {
                   if(!is.na(varname))
                     object@dispEst[[varname]] <- value
                   else
                     object@dispEst <- value
                   validObject( object )
                   object
                 })


XBSeqDataSet <- function(counts, bgcounts, conditions, sizeFactors=NULL, ...)
{
  counts <- as.matrix(counts)
  bgcounts <-  as.matrix(bgcounts)
  conditions <- as.factor(conditions)
  if(any(round(counts) != counts) | any(round(bgcounts) != bgcounts))
    stop("The input data have to be integer!")
  mode( counts ) <- "integer"
  mode( bgcounts ) <- "integer"
  if( nrow( counts )!= nrow( bgcounts ) ){
    if(is.null(rownames(counts)) | is.null(rownames(bgcounts)))
      stop("Please provide gene symbols or other unique ids as rownames")
    else if(sum(is.na(as.numeric(rownames(counts)))) != 0 | sum(is.na(as.numeric(rownames(bgcounts)))) != 0)
      stop('Please provide meaningful names as rownames rather than arabic numerals')
    MissedRecord <- which( rownames( counts) %in% setdiff( rownames( counts ),rownames( bgcounts ) ) )
    bgcounts <- insertRow( bgcounts, repmat(apply(bgcounts,2,mean),length(MissedRecord),1) ,MissedRecord)
  }
  if( is.null( sizeFactors ) ) {
    sizeFactors <- rep( NA_real_, ncol(counts) )
  }
  assays <- list(counts=counts, bgcounts=bgcounts)
  colData <- data.frame(conditions = conditions)
  rownames(colData) <- colnames(counts)
  colData <- DataFrame(colData, row.names=rownames(colData))
  se <- SummarizedExperiment(assays, colData = colData)
  stopifnot(length(conditions) == ncol(counts))
  rvft <- rep(NA_character_, length(levels(conditions)))
  XB <- DESeqDataSet(se, formula(~conditions))
  XB <- new("XBSeqDataSet",
            XB, 
            fitInfo = new.env(hash=TRUE),
            dispTable = rvft,
            conditions = conditions)
  return(XB)
}







setGeneric('estimateRealCount', function(object) standardGeneric('estimateRealCount'))


setMethod('estimateRealCount', signature(object = 'XBSeqDataSet'),
          function(object){
            signal <- assay(object, 1) - assay(object, 2)
            signal[signal < 0] <- 0
            assay(object,3) <- signal
            object
          })


setMethod('counts', signature(object = 'XBSeqDataSet'),
          function(object, slot = 3, normalized = FALSE){
            if(length(assays(object)) == 2 & slot == 3)
              stop('Only two assays exist. Call "estimateRealCount" first')
            if (!normalized) {
              return(assay(object, slot))
            } 
            else if (is.null(sizeFactors(object)) | any(is.na(sizeFactors(object)))) {
              stop("first calculate size factors, add normalizationFactors, or set normalized=FALSE")
            } else {
              return(t(t(assay(object,slot)) / sizeFactors(object)))
            }
          })

setGeneric("estimateSCV",function(object, ...) standardGeneric('estimateSCV'))


setMethod("estimateSCV", signature(object="XBSeqDataSet"),
          function( object, method = c( "pooled", "per-condition", "blind" ),
                    sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
                    fitType = c("local", "parametric"),
                    locfit_extra_args=list(), lp_extra_args=list(), ... )
          {
            t1 <- Sys.time()
            stopifnot( is( object, "XBSeqDataSet" ) )
            if(any(c(is.null(sizeFactors(object)))))
              stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
            method <- match.arg( method )
            sharingMode <- match.arg( sharingMode )
            fitType <- match.arg( fitType )
            if( length(list(...)) != 0 )
              warning( "in estimateSCV: Ignoring extra argument(s)." )
            if( sharingMode == "gene-est-only" && length(conditions(object))/length(levels(conditions(object))) <= 2 )
              warning( "in estimateSCV: sharingMode=='gene-est-only' will cause inflated numbers of false positives unless you have many replicates." )
            #Remove results from previous fits
            index <- which(! colnames(dispEst(object)) %in% paste( "disp", object@dispTable, sep="_" ))
            if(length(index)){
              dispEst(object) <- dispEst(object, index)
              object@dispTable <- character()
              object@fitInfo = new.env( hash=TRUE )
            }
            if( method == "blind" ) {
              data <- getCountParams(counts(object), sizeFactors(object))
              data_var <- getSignalVars(counts(object, 1), counts(object, 2))
              SCVf <- getSCV(data$baseMean,
                             data_var, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
              object@fitInfo[[ "blind" ]] <- list(
                perGeneSCVEsts = SCVf$SCV,
                SCVFunc = SCVf$SCVfunc,
                fittedSCVEsts = SCVf$SCVfunc( data$baseMean ),
                df = ncol(counts(object)) - 1,
                sharingMode = sharingMode )
              a <- rep( "blind", length( levels( conditions(object) ) ) )
              names(a) <- levels( conditions(object) )
              object@dispTable <- a
            } else if( method == "per-condition" ) {
              replicated <- names( which( tapply( conditions(object), conditions(object), length ) > 1 ) )
              if( length( replicated ) < 1 )
                stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, if you have crossed factors." )
              nonreplicated <- names( which( tapply( conditions(object), conditions(object), length ) == 1 ) )
              overall_basemeans <- rowMeans( counts( object, normalized=TRUE ) )
              for( cond in replicated ) {
                cols <- conditions(object)==cond
                data <- getCountParams(counts(object)[ , cols ], sizeFactors(object)[ cols ] )
                data_var <- getSignalVars(counts(object, 1)[, cols], counts(object, 2)[, cols])
                SCVf <- getSCV( data$baseMean,
                                data_var, sizeFactors(object)[cols], fitType, locfit_extra_args, lp_extra_args )
                object@fitInfo[[ cond ]] <- list(
                  perGeneSCVEsts = SCVf$SCV,
                  SCVFunc = SCVf$SCVfunc,
                  fittedSCVEsts = SCVf$SCVfunc( overall_basemeans ),     # Note that we do not use bmv$baseMean here
                  df = sum(cols) - 1,
                  sharingMode = sharingMode ) }
              
              object@dispTable <- sapply( levels(conditions(object)), function( cond )
                ifelse( cond %in% replicated, cond, "max" ) )
              
            } else if( method == "pooled" ) {
              conds <- conditions(object)
              if( !any( duplicated( conds ) ) )
                stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors." )
              data <- getCountParamsPooled( counts(object), sizeFactors(object), conds )
              baseMeans <- data$baseMean
              data_var <- getSignalVars(counts(object, 1), counts(object, 2))
              SCVf <- getSCV(data$baseMean,
                             data$baseVar, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
              df <- ncol(counts(object)) - length(unique(conds))
              object@fitInfo[[ "pooled" ]] <- list(
                perGeneSCVEsts = SCVf$SCV,
                SCVFunc = SCVf$SCVfunc,
                fittedSCVEsts = SCVf$SCVfunc( baseMeans ),
                df = df,
                sharingMode = sharingMode )
              a <- rep( "pooled", length( levels( conditions(object) ) ) )
              names(a) <- levels( conditions(object) )
              dispTable(object) = a
            } else
              stop(sprintf("Invalid method '%s'.", method))
            for( n in ls(object@fitInfo) )
              dispEst(object, paste( "disp", n, sep="_" ) ) <-
              switch( sharingMode,
                      `fit-only`      = object@fitInfo[[ n ]]$fittedSCVEsts,
                      `gene-est-only` = {
                        a <- object@fitInfo[[ n ]]$perGeneSCVEsts
                        a[ is.nan(a) ] <- 0
                        pmax( a, 1e-8 ) },
                      `maximum`       = pmax( object@fitInfo[[ n ]]$fittedSCVEsts, object@fitInfo[[ n ]]$perGeneSCVEsts, na.rm=TRUE ),
                      stop(sprintf("Invalid sharingMode '%s'.", sharingMode))
              ) ## switch
            
            if( "max" %in% object@dispTable )
              dispEst(object, "disp_max") <- do.call( pmax,
                                                      c( dispEst(object, which(colnames(dispEst(object)) %in% paste( "disp", object@dispTable, sep="_" )) ), na.rm=TRUE ) )
            validObject( object )
            return(object)
          })


XBSeqTest <- function(XB, condA, condB, pvals_only=FALSE, method = c('NP', 'MLE'))
{
  stopifnot( is( XB, "XBSeqDataSet" ) )
  if( all( is.na( dispTable(XB) ) ) )
    stop( "Call 'estimateSCV' first." )
  if( dispTable(XB)[condA] == "blind") {
    if( fitInfo( XB, "blind" )$sharingMode != "fit-only" )
      warning( 'You have used \'method="blind"\' in estimateSCV without also setting \'sharingMode="fit-only"\'. This will not yield useful results.' )
  }
  stopifnot(condA %in% levels(conditions(XB)))
  stopifnot(condB %in% levels(conditions(XB)))
  colA <- conditions(XB)==condA
  colB <- conditions(XB)==condB
  
  rawScvA <- dispEst(XB, paste( "disp", dispTable(XB)[condA], sep="_" ))
  rawScvB <- dispEst(XB, paste( "disp", dispTable(XB)[condB], sep="_" ))
  
  pval <- XBSeqTestForMatrices(
    counts(XB)[,colA],
    counts(XB)[,colB],
    counts(XB, slot=2)[,colA],
    counts(XB, slot=2)[,colB],
    sizeFactors(XB)[colA],
    sizeFactors(XB)[colB],
    rawScvA,
    rawScvB, 
    method = method)
  if( pvals_only )
    return(pval)
  else {
    data <- getCountParams( counts(XB), sizeFactors(XB)[colA|colB] )
    dataA <- getCountParams( counts(XB)[,colA], sizeFactors(XB)[colA] )
    dataB <- getCountParams( counts(XB)[,colB], sizeFactors(XB)[colB] )
    return(data.frame(
      id    = rownames( counts(XB) ),
      baseMean  = data$baseMean,
      baseMeanA = dataA$baseMean,
      baseMeanB = dataB$baseMean,
      foldChange = dataB$baseMean / dataA$baseMean,
      log2FoldChange = log2( dataB$baseMean / dataA$baseMean ),
      pval = pval,
      padj = p.adjust( pval, method="BH" ),
      stringsAsFactors = FALSE )) 
  }
}


# a wrapper function once and for all
XBSeq <- function(counts, bgcounts, conditions, method='pooled', sharingMode='maximum', fitType='local', pvals_only=FALSE, paraMethod='NP'){
  if(!is.factor(conditions))
    conditions <- as.factor(conditions)
  XB <- XBSeqDataSet(counts, bgcounts, conditions)
  XB <- estimateRealCount(XB)
  XB <- estimateSizeFactors(XB)
  XB <- estimateSCV(XB, method=method, sharingMode=sharingMode, fitType=fitType)
  Teststas <- XBSeqTest(XB, levels(conditions)[1L], levels(conditions)[2L], pvals_only=pvals_only, method = paraMethod)
  Teststas
}








getCountParams <- function( counts, sizeFactors ) {
  
  # Divides the counts by sizeFactors and calculates the estimates for
  # base means and variances for each gene
  
  data.frame(
    baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
    baseVar = rowVars( t( t(counts) / sizeFactors ) ) )
}


getCountParamsPooled <- function( counts, sizeFactors, conditions ) {
  
  basecounts <- t( t(counts) / sizeFactors )
  replicated_sample <- conditions %in% names(which(table(conditions)>1))
  df <- sum(replicated_sample) - length( unique( conditions[ replicated_sample ] ) )
  
  data.frame(
    baseMean = rowMeans( basecounts ),
    baseVar =
      rowSums(
        sapply(
          tapply(
            ( seq_len(ncol(counts)) )[ replicated_sample ],
            factor( conditions[ replicated_sample ] ),
            function(cols)
              rowSums( ( basecounts[,cols] - rowMeans(basecounts[,cols]) )^2 ) ),
          identity ) ) / df )
}

getSCV <- function( means,
                    variances, sizeFactors, fitType = c( "parametric", "local" ),
                    locfit_extra_args=list(), lp_extra_args=list(), adjustForBias=TRUE ) {
  
  fitType <- match.arg( fitType )
  
  xim <- mean( 1/sizeFactors )
  SCVAll <- ( variances - xim * means ) / means^2
  
  variances <- variances[ means > 0 ]
  SCV <- SCVAll[ means > 0 ]
  means <- means[ means > 0 ]
  
  if( adjustForBias )
    SCV <- adjustScv( SCV, length( sizeFactors ) )
  
  if( fitType == "local" ) {
    
    fit <- do.call( "locfit", c(
      list(
        variances ~ do.call( "lp", c( list( log(means) ), lp_extra_args ) ),
        family = "gamma" ),
      locfit_extra_args ) )
    
    rm( means )
    rm( variances )
    
    if( adjustForBias )
      ans <- function( q )
        adjustScv(
          pmax( ( predict_helper( fit, log(q) ) - xim * q ) / q^2, 1e-8 ),
          length(sizeFactors) )
    else
      ans <- function( q )
        pmax( ( predict_helper( fit, log(q) ) - xim * q ) / q^2, 1e-8 )
    
    # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
    # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.
    
  } else if( fitType == "parametric" ) {
    
    ans <- parametricscvFit( means, SCV )
    
  } else
    stop( "Unknown fitType." )
  
  attr( ans, "fitType" ) <- fitType
  list( SCV=SCVAll, SCVfunc=ans )
}


parametricscvFit <- function( means, disps )
{
  coefs <- c( .1, 1 )
  iter <- 0
  while(TRUE) {
    residuals <- disps / ( coefs[1] + coefs[2] / means )
    good <- which( (residuals > 1e-4) & (residuals < 15) )
    fit <- glm( disps[good] ~ I(1/means[good]),
                family=Gamma(link="identity"), start=coefs )
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if( !all( coefs > 0 ) )
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateSCV')" )
    if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
      break
    iter <- iter + 1
    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break }
  }
  
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function( q )
    coefs[1] + coefs[2] / q
  attr( ans, "coefficients" ) <- coefs
  ans
}


predict_helper <- function( fit, x )
{
  # A wrapper around predict to avoid the issue that predict.locfit cannot
  # propagate NAs and NaNs properly.
  
  res <- rep.int( NA_real_, length(x) )
  res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
  res
}


XBSeqTestForMatrices <- function( countsA, countsB, bgcountsA, bgcountsB, sizeFactorsA, sizeFactorsB,
                                  SCVA, SCVB , method = c('NP', 'MLE'))
{
  method <- match.arg(method, c('NP', 'MLE'))
  if(ncol(countsA) < 5 & method == 'MLE')
    warning('Non-parametric estimation method is recommended for experiments with replicates smaller than 5')
  
  kAs <- apply(countsA, 1, sum)
  kBs <- apply(countsB, 1, sum)
  
  mus <- rowMeans(cbind(
    t( t( countsA ) / sizeFactorsA ),
    t( t( countsB ) / sizeFactorsB)))
  
  if(method == 'NP'){
    signalmuA <- mus*sum(sizeFactorsA)
    signalmuB <- mus*sum(sizeFactorsB)
    
    signalVarsA <- pmax( mus * sum( sizeFactorsA ) + SCVA * mus^2 * sum(sizeFactorsA^2),
                         mus * sum( sizeFactorsA ) * (1+1e-8) )
    signalVarsB <- pmax( mus * sum( sizeFactorsB ) + SCVB * mus^2 * sum(sizeFactorsB^2),
                         mus * sum( sizeFactorsB ) * (1+1e-8) )
    
    sizeA <- signalmuA^2/(signalVarsA-signalmuA)
    sizeB <- signalmuB^2/(signalVarsB-signalmuB)
  }
  else{
    musA <- mus*mean(sizeFactorsA)
    musB <- mus*mean(sizeFactorsB)
    
    VarsA <- pmax(mus * mean(sizeFactorsA) + SCVA * mus^2 * mean(sizeFactorsA^2),
                  mus * mean(sizeFactorsA) * (1+1e-8))
    VarsB <- pmax(mus * mean(sizeFactorsB) + SCVB * mus^2 * mean(sizeFactorsB^2),
                  mus * mean(sizeFactorsB) * (1+1e-8))
    
    lambda <- rowMeans(cbind(
      t(t(bgcountsA) / sizeFactorsA),
      t(t(bgcountsB) / sizeFactorsB)))
    
    lambdaA <- lambda*mean(sizeFactorsA)
    lambdaB <- lambda*mean(sizeFactorsB)
    
    
    ParamsA <- sapply(1:nrow(countsA), function(i) estimation_param_PoissonNB_MLE(countsA[i,] + bgcounts[i,],
                                                                                  bgcountsA[i,],
                                                                                  musA[i]^2/(VarsA[i]-musA[i]),
                                                                                  (VarsA[i]-musA[i])/musA[i],
                                                                                  lambdaA[i]
    ))
    ParamsA <- matrix(unlist(ParamsA), ncol=3, byrow = TRUE)
    ParamsB <- sapply(1:nrow(countsB), function(i) estimation_param_PoissonNB_MLE(countsB[i,] + bgcounts[i,],
                                                                                  bgcountsB[i,],
                                                                                  musB[i]^2/(VarsB[i]-musB[i]),
                                                                                  (VarsB[i]-musB[i])/musB[i],
                                                                                  lambdaB[i]
    ))
    ParamsB <- matrix(unlist(ParamsB), ncol=3, byrow = TRUE)
    
    sizeA <- ParamsA[,1] * sum(sizeFactorsA)
    sizeB <- ParamsB[,1] * sum(sizeFactorsB)
    
    signalmuA <- ParamsA[,1] * ParamsB[,2] * sum(sizeFactorsA)
    signalmuB <- ParamsB[,1] * ParamsB[,2] * sum(sizeFactorsB)
  }
  
  sapply( seq(along=kAs), function(i) {
    
    if(kAs[i] == 0 & kBs[i] == 0)
      return( NA )
    
    # probability of all possible counts sums with the same total count:
    ks <- 0 : (kAs[i] + kBs[i])
    ps <- dnbinom(ks, mu = signalmuA[i], size = sizeA[i]) *
      dnbinom( kAs[i] + kBs[i] - ks, mu = signalmuB[i], size = sizeB[i])
    
    # probability of observed count sums:
    pobs <- dnbinom( kAs[i], mu = signalmuA[i], size = sizeA[i]) *
      dnbinom( kBs[i], mu = signalmuB[i], size = sizeB[i])
    
    if( kAs[i] * sum( sizeFactorsB ) < kBs[i] * sum( sizeFactorsA ) )
      numer <- ps[ 1 : (kAs[i]+1) ]
    else
      numer <- ps[ (kAs[i]+1) : length(ps) ]
    min( 1, 2 * sum(numer) / sum(ps) )
  })
}


insertRow <- function(existingDF, newrow, r) {
  for(i in 1:length(r))
    existingDF[seq(r[i]+1,nrow(existingDF)+1),] <- existingDF[seq(r[i],nrow(existingDF)),]
  existingDF[r[i],] <- newrow[i]
  existingDF
}


getSignalVars<-function( counts, bgcounts){
  if(!is.matrix(counts))
    counts <- as.matrix(counts)
  if(!is.matrix(bgcounts))
    bgcounts <- as.matrix(bgcounts)
  bgsizeFactors <- estimateSizeFactorsForMatrix(bgcounts)
  lambda <- rowMeans( sweep(bgcounts, 2, bgsizeFactors, '/' ) )
  observe_param <- getCountParams(counts, estimateSizeFactorsForMatrix(counts))
  observe_sf <- estimateSizeFactorsForMatrix(counts)
  tho <- c()
  for(i in 1:nrow(counts) )
  {
    if(lambda[i]==0){
      temp <- 0
      if ( length( which(counts[i,]==0) ) >= ncol(counts)/2 )
        temp <- NA
    }
    else {
      if(length( which(counts[i,]==0) ) >= ncol(counts)/2 ){
        temp <- NA
      }
      else {
        temp2 <- cor( counts[i,]*observe_sf, bgcounts[i,]*bgsizeFactors )
        if( is.na(temp2))
          temp <- 0
        else
          temp <-temp2
      }
    }
    tho <- c(tho, temp)
  }
  fullvar <- observe_param$baseVar
  signalvar <- fullvar + lambda - 2*tho*sqrt(fullvar)*sqrt(lambda)
  as.matrix(signalvar)
}


prepareScvBiasCorrectionFits <- function( maxnrepl=15, mu=100000, ngenes=10000,
                                          true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ) )
  lapply( 2:maxnrepl, function( m ) {
    est_raw_scv <- sapply( true_raw_scv, function( alpha ) {
      k <- matrix( rnbinom( ngenes*m, mu=mu, size=1/alpha ), ncol=m )
      k <- k[ rowSums(k)>0, ]
      mean( rowVars(k) / rowMeans(k)^2 ) } )
    locfit( true_raw_scv ~ lp( est_raw_scv, nn=.2 ) ) } )


adjustScv <- function( scv, nsamples ) {
  stopifnot( nsamples > 1 )
  if(!exists(("scvBiasCorrectionFits")))
    data( "scvBiasCorrectionFits" )
  if( nsamples - 1 > length( scvBiasCorrectionFits ) )
    scv
  else
    ifelse( scv > .02,
            pmax( predict_helper( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv ),
            scv )   # For scv < .02, our fit is too coarse, but no correction seems necessary anyway
}


Loglikhood <- function(counts, bgcounts){
  function(para) {
    alpha <- para[1]
    beta <- para[2]
    lambda <- para[3] 
    -sum(ddelap(counts, alpha = alpha, beta = beta, lambda = lambda, log = TRUE), dpois(bgcounts, lambda = lambda, log = TRUE))
  }
}


estimation_param_PoissonNB_MLE <- function(counts, bgcounts, alpha, beta, lambda){
  if(any(is.na(c(alpha, beta, lambda))))
    list(alpha = alpha, beta = beta, lambda = lambda)
  else{
    mle <- try(optim(c(alpha, beta, lambda), Loglikhood(counts, bgcounts),
                     method = 'L-BFGS-B', lower=c(0,0,0)), silent = TRUE)
    if(class(mle) == 'try-error')
      list(alpha = alpha, beta = beta, lambda = lambda)
    else{
      mle <- optim(c(alpha, beta, lambda), Loglikhood(counts, bgcounts),
                   method = 'L-BFGS-B', lower=c(0,0,0))
      if (mle$convergence>0 | any(is.na(mle$par)))
        list(alpha = alpha, beta = beta, lambda = lambda)
      else
        list(alpha = mle$par[1], beta = mle$par[2], lambda = mle$par[3])
    }
  }
}