require('pracma')
require('matrixStats')
require('locfit')
require('Biobase')
require('ggplot2')

setGeneric("estimateSizeFactorsXBSeq",function(object, ...) standardGeneric('estimateSizeFactorsXBSeq'))
setMethod("estimateSizeFactorsXBSeq", signature(object="XBSeqDataSet"),
          function( object, locfunc=median, ... ) {
            if( length(list(...)) != 0 )
              warning( "in estimateSizeFactors: Ignoring extra argument(s)." )
            sizeFactors(object) <- estimateSizeFactorsForMatrixXBSeq( counts(object), locfunc )
            object
          })


setGeneric("estimateSCV",function(object, ...) standardGeneric('estimateSCV'))
setMethod("estimateSCV", signature(object="XBSeqDataSet"),
          function( object, observe, background, method = c( "pooled", "pooled-CR", "per-condition", "blind" ),
                    sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
                    fitType = c( "parametric", "local" ),
                    locfit_extra_args=list(), lp_extra_args=list(),
                    modelFrame = NULL, modelFormula = count ~ condition, ... )
          {
            stopifnot( is( object, "XBSeqDataSet" ) )
            if( any( c( is.na( sizeFactors(object) ) ) ) )
              stop( "NAs found in size factors. Have you called already 'estimateSizeFactorsXBSeq'?" )
            method <- match.arg( method )
            sharingMode <- match.arg( sharingMode )
            fitType <- match.arg( fitType )
            if( length(list(...)) != 0 )
              warning( "in estimateSCV: Ignoring extra argument(s)." )
            if( object@multivariateConditions && ! method %in% c( "blind", "pooled", "pooled-CR" ) )
              stop( "You have specified multivariate conditions (i.e., passed a data frame with conditions). In this case, you cannot use method 'per-condition'." )
            if( sharingMode == "gene-est-only" && length(pData(object)$condition)/length(levels(pData(object)$condition)) <= 2 )
              warning( "in estimateDispersions: sharingMode=='gene-est-only' will cause inflated numbers of false positives unless you have many replicates." )
            #Remove results from previous fits
            fData(object) <- fData(object)[ , ! colnames(fData(object)) %in% paste( "disp", object@dispTable, sep="_" ), drop=FALSE ]
            object@dispTable <- character()
            object@fitInfo = new.env( hash=TRUE )
            if( method == "blind" ) {
              data <- getCountParams( counts(object), sizeFactors(object) )
              data_var <- getsignalVars(observe, background)
              SCVf <- getSCV( data$baseMean,
                              data_var, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
              object@fitInfo[[ "blind" ]] <- list(
                perGeneSCVEsts = SCVf$SCV,
                SCVFunc = SCVf$SCVfunc,
                fittedSCVEsts = SCVf$SCVfunc( data$baseMean ),
                df = ncol(counts(object)) - 1,
                sharingMode = sharingMode )
              if( object@multivariateConditions )
                dispTable(object) <- c( "_all" = "blind" )
              else {
                a <- rep( "blind", length( levels( conditions(object) ) ) )
                names(a) <- levels( conditions(object) )
                object@dispTable <- a }

            } else if( method == "per-condition" ) {
              replicated <- names( which( tapply( conditions(object), conditions(object), length ) > 1 ) )
              if( length( replicated ) < 1 )
                stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors." )
              nonreplicated <- names( which( tapply( conditions(object), conditions(object), length ) == 1 ) )
              overall_basemeans <- rowMeans( counts( object, normalized=TRUE ) )
              for( cond in replicated ) {
                cols <- conditions(object)==cond
                data <- getCountParams( counts(object)[ , cols ], sizeFactors(object)[ cols ] )
                data_var <- getsignalVars(observe[,cols], background[,cols])
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

            } else if( method == "pooled" || method == "pooled-CR" ) {

              if( method == "pooled" ) {

                if( object@multivariateConditions ) {
                  if( is.null( modelFrame ) )
                    modelFrame <- pData(object)[ , colnames(pData(object)) != "sizeFactor" ]
                  conds <- modelToConditionFactor( modelFrame ) }
                else
                  conds <- conditions(object)
                if( !any( duplicated( conds ) ) )
                  stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors." )
                data <- getCountParamsPooled( counts(object), sizeFactors(object), conds )
                baseMeans <- data$baseMean
                data_var <- getsignalVars(observe, background)
                SCVf <- getSCV( data$baseMean,
                                data$baseVar, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
                df <- ncol(counts(object)) - length(unique(conds))

              } else {  # method == "pooled-CR"
                if( is.null( modelFrame ) )
                  modelFrame <- pData(object)[ , colnames(pData(object)) != "sizeFactor", drop=FALSE ]
                baseMeans <- rowMeans( counts( object, normalized=TRUE ) )
                SCVf <- getSCVCoxReid( counts(object), modelFormula, modelFrame,
                                sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
                df <- NA
              }
              object@fitInfo[[ "pooled" ]] <- list(
                perGeneSCVEsts = SCVf$SCV,
                SCVFunc = SCVf$SCVfunc,
                fittedSCVEsts = SCVf$SCVfunc( baseMeans ),
                df = df,
                sharingMode = sharingMode )

              dt = if( object@multivariateConditions ) c( "_all" = "pooled" ) else character(0)

              if("condition" %in% colnames(pData(object))) {
                a <- rep( "pooled", length( levels( conditions(object) ) ) )
                names(a) <- levels( conditions(object) )
                dt = c(dt, a)
              }
              dispTable(object) = dt

            } else
              stop(sprintf("Invalid method '%s'.", method))

            for( n in ls(object@fitInfo) )
              fData(object)[[ paste( "disp", n, sep="_" ) ]] <-
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
              fData(object)[["disp_max"]] <- do.call( pmax,
                                                      c( fData(object)[ , colnames(fData(object)) %in% paste( "disp", object@dispTable, sep="_" ), drop=FALSE ], na.rm=TRUE ) )

            validObject( object )
            object
          })


XBSeqTest <- function( XB, condA, condB, pvals_only=FALSE )
{
  stopifnot( is( XB, "XBSeqDataSet" ) )
  if( all( is.na( dispTable(XB) ) ) )
    stop( "Call 'estimateSCV' first." )
  if( dispTable(XB)[condA] == "blind") {
    if( fitInfo( XB, "blind" )$sharingMode != "fit-only" )
      warning( 'You have used \'method="blind"\' in estimateSCV without also setting \'sharingMode="fit-only"\'. This will not yield useful results.' )
  }
  stopifnot( condA %in% levels(conditions(XB)) )
  stopifnot( condB %in% levels(conditions(XB)) )
  colA <- conditions(XB)==condA
  colB <- conditions(XB)==condB

  rawScvA <- fData(XB)[ , paste( "disp", dispTable(XB)[condA], sep="_" ) ]
  rawScvB <- fData(XB)[ , paste( "disp", dispTable(XB)[condB], sep="_" ) ]

  pval <- XBSeqTestForMatrices(
    counts(XB)[,colA],
    counts(XB)[,colB],
    sizeFactors(XB)[colA],
    sizeFactors(XB)[colB],
    rawScvA,
    rawScvB )

  if( pvals_only )
    pval
  else {
    data <- getCountParams( counts(XB), sizeFactors(XB)[colA|colB] )
    dataA <- getCountParams( counts(XB)[,colA], sizeFactors(XB)[colA] )
    dataB <- getCountParams( counts(XB)[,colB], sizeFactors(XB)[colB] )
    data.frame(
      id    = rownames( counts(XB) ),
      baseMean  = data$baseMean,
      baseMeanA = dataA$baseMean,
      baseMeanB = dataB$baseMean,
      foldChange = dataB$baseMean / dataA$baseMean,
      log2FoldChange = log2( dataB$baseMean / dataA$baseMean ),
      pval = pval,
      padj = p.adjust( pval, method="BH" ),
      stringsAsFactors = FALSE ) }
}


# a wrapper function once and for all
XBSeq <- function( observe, background, conditions, method='pooled', sharingMode='maximum', fitType='parametric', pvals_only=FALSE ){
  Signal <- estimateRealcount(observe, background)
  XB <- newXBSeqDataSet(Signal,conditions)
  XB <- estimateSizeFactorsXBSeq( XB )
  XB <-estimateSCV( XB, observe, background, method=method, sharingMode=sharingMode, fitType=fitType )
  Teststas <- XBSeqTest( XB, levels(conditions)[1L], levels(conditions)[2L], pvals_only=pvals_only )
  Teststas
}

estimateSizeFactorsForMatrixXBSeq <- function( counts, locfunc = median )
{
  loggeomeans <- rowMeans( log(counts) )
  apply( counts, 2, function(cnts)
    exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
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


modelToConditionFactor <- function( modelMatrix ) {
  if(!is.matrix(modelMatrix))
    modelMatrix<-as.matrix(modelMatrix)
  nr = nrow(modelMatrix)
  if(nr<2)
    stop("nrow(modelMatrix) must be >=2.")

  mmconds <- 1:nr
  for( i in 2:nr )
    for( j in 1:(i-1) )
      if( all( modelMatrix[i,] == modelMatrix[j,] ) ) {
        mmconds[i] = mmconds[j]
        break
      }
  factor( as.integer( factor( mmconds ) ) )
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


getSCVCoxReid <- function( counts, modelFormula, modelFrame,
                           sizeFactors, fitType = c( "parametric", "local" ),
                           locfit_extra_args=list(), lp_extra_args=list(), initialGuess=.1 )
{
  if( as.character(modelFormula[2]) == "count" )
    modelFormula <- modelFormula[-2]   # the '[-2]' removes the lhs, i.e., the response
  mm <- model.matrix( modelFormula, modelFrame )
  disps <- apply( counts, 1, function( y ) {
    fit <- try(
      glm.fit( mm, y, family=MASS::negative.binomial( initialGuess ), offset=log(sizeFactors) ),
      silent=TRUE )
    if( inherits( fit, "try-error" ) )
      NA
    else {
      if( df.residual(fit) == 0 )
        stop( "No residual degrees of freedom. Most likely the design is lacking sufficient replication." )
      exp(
        optimize(
          function(logalpha)
            LogLikelihood( exp(logalpha), mm, y, fitted.values(fit) ),
          log( c( 1e-11, 1e5 ) ),
          maximum=TRUE
        )$maximum ) } } )

  means <- colMeans( t(counts) / sizeFactors )
  xim <- mean( 1/sizeFactors )

  if( fitType == "local" ) {

    fit <- do.call( "locfit", c(
      list(
        disps[means>0] ~ do.call( "lp", c( list( log(means[means>0]) ), lp_extra_args ) ),
        family = "gamma" ),
      locfit_extra_args ) )

    rm( means )

    ans <- function( q )
      pmax( ( predict_helper( fit, log(q) ) - xim * q ) / q^2, 1e-8 )

    # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
    # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.

  } else if( fitType == "parametric" )

    ans <- parametricDispersionFit( means, disps )

  else
    stop( "Unkknown fitType." )

  attr( ans, "fitType" ) <- fitType
  list( disps=disps, dispFunc=ans )

}


fitInfo <- function( XB, name=NULL )
{
  stopifnot( is( XB, "XBSeqDataSet" ) )
  if( length( ls( XB@fitInfo ) ) == 0 )
    stop( "No fits available. Call 'estimateDispersions' first." )
  if( length( ls( XB@fitInfo ) ) > 1 && is.null(name) )
    stop( "More than one fitInfo object available. Specify by name. (See 'ls(XB@fitInfo)' for a list.)" )
  if( length( ls( XB@fitInfo ) ) == 1 && is.null(name) )
    name = ls( XB@fitInfo )[ 1 ]
  XB@fitInfo[[ name]]
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
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
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


LogLikelihood <- function( disp, mm, y, muhat )
{
  # calculate the log likelihood:
  if(length(disp) != length(y)){
    disp <- rep(disp, length(y))
  }

  ll <- sum( sapply( seq(along=y), function(i)
    dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) ) )

  # transform the residuals, i.e., y - muhat, to the linear
  # predictor scale by multiplying them with the derivative
  # of the link function, i.e., by 1/muhat, and add this to the
  # linear predictors, log(muhat), to get the predictors that
  # are used in the IWLS regression
  z <- log(muhat) + ( y - muhat ) / muhat

  # the variance function of the NB is as follows
  v0 <- muhat + disp * muhat^2

  # transform the variance vector to linear predictor scale by
  # multiplying with the squared derivative of the link to
  # get the (reciprocal) weights for the IWLS
  w <- 1 / ( ( 1 / muhat )^2 * v0 )

  # All we need from the IRLS run is the QR decomposition of
  # its matrix
  qrres <- qr( mm*sqrt(w) )

  # from it, we extract we leverages and calculate the Cox-Reid
  # term:
  cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )

  # return the profile log likelihood:
  ll - cr
}


predict_helper <- function( fit, x )
{
  # A wrapper around predict to avoid the issue that predict.locfit cannot
  # propagate NAs and NaNs properly.

  res <- rep.int( NA_real_, length(x) )
  res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
  res
}


XBSeqTestForMatrices <- function( countsA, countsB, sizeFactorsA, sizeFactorsB,
                                  SCVA, SCVB )
{
  kAs <- rowSums( countsA )
  kBs <- rowSums( countsB )

  mus <- rowMeans( cbind(
    t( t( countsA ) / sizeFactorsA ),
    t( t( countsB ) / sizeFactorsB ) ) )

  signalmuA <- mus*sum(sizeFactorsA)
  signalmuB <- mus*sum(sizeFactorsB)

  signalVarsA <- pmax( mus * sum( sizeFactorsA ) + SCVA * mus^2 * sum(sizeFactorsA^2),
                     mus * sum( sizeFactorsA ) * (1+1e-8) )
  signalVarsB <- pmax( mus * sum( sizeFactorsB ) + SCVB * mus^2 * sum(sizeFactorsB^2),
                     mus * sum( sizeFactorsB ) * (1+1e-8) )

  sapply( seq(along=kAs), function(i) {

    if( kAs[i] == 0 & kBs[i] == 0 )
      return( NA )

    # probability of all possible counts sums with the same total count:
    ks <- 0 : ( kAs[i] + kBs[i] )
    ps <- dnbinom(                   ks, mu = signalmuA[i], size = signalmuA[i]^2/(signalVarsA[i]-signalmuA[i]) ) *
      dnbinom( kAs[i] + kBs[i] - ks, mu = signalmuB[i], size = signalmuB[i]^2/(signalVarsB[i]-signalmuB[i]) )

    # probability of observed count sums:
    pobs <- dnbinom( kAs[i], mu = signalmuA[i], size = signalmuA[i]^2/(signalVarsA[i]-signalmuA[i]) ) *
      dnbinom( kBs[i], mu = signalmuB[i], size = signalmuB[i]^2/(signalVarsB[i]-signalmuB[i]) )

    #stopifnot( na.omit(pobs == ps[ kAs[i]+1 ]) )
    if( kAs[i] * sum( sizeFactorsB ) < kBs[i] * sum( sizeFactorsA ) )
      numer <- ps[ 1 : (kAs[i]+1) ]
    else
      numer <- ps[ (kAs[i]+1) : length(ps) ]
    min( 1, 2 * sum(numer) / sum(ps) )
  } )
}

makeExampleXBSeqDataSet <- function( )
{
  ngenes <- 10000
  q0 <- rexp( ngenes, rate=1/250 )
  is_DE <- runif( ngenes ) < .3
  lfc <- rnorm( ngenes, sd=2 )
  q0A <- ifelse( is_DE, q0 * 2^(  lfc/2 ), q0 )
  q0B <- ifelse( is_DE, q0 * 2^( -lfc/2 ), q0 )
  true_sf <- c( 1., 1.3, .7, .9, 1.6 )
  conds <- c( "A", "A", "B", "B", "B" )
  m <- t( sapply( seq_len(ngenes), function(i)
    sapply( 1:5, function( j )
      rnbinom( 1, mu = true_sf[j] * ifelse( conds[j]=="A", q0A[i], q0B[i] ),
               size = 1/.2 ) ) ) )
  colnames(m) <- c( "A1", "A2", "B1", "B2", "B3" )
  rownames(m) <- paste( "gene", seq_len(ngenes),
                        ifelse( is_DE, "T", "F" ), sep="_" )
  newXBSeqDataSet( m, conds )
}

newXBSeqDataSetFromHTSeqCount <- function( sampleTable, directory="." )
{
  l <- lapply( as.character( sampleTable[,2] ), function(fn)
    read.table( file.path( directory, fn ) ) )
  if( ! all( sapply( l, function(a) all( a$V1 == l[[1]]$V1 ) ) ) )
    stop( "Gene IDs (first column) differ between files." )
  tbl <- sapply( l, function(a) a$V2 )
  rownames(tbl) <- l[[1]]$V1
  colnames(tbl) <- sampleTable[,1]
  specialRows <- rownames(tbl) %in% c( "no_feature", "ambiguous",
                                       "too_low_aQual", "not_aligned", "alignment_not_unique" )
  tbl <- tbl[ !specialRows, ]
  if( ncol(sampleTable) == 3 )
    newXBSeqDataSet( tbl, sampleTable[,3] )
  else
    newXBSeqDataSet( tbl, sampleTable[,-(1:2)] )
}


insertRow <- function(existingDF, newrow, r) {
  for(i in 1:length(r))
    existingDF[seq(r[i]+1,nrow(existingDF)+1),] <- existingDF[seq(r[i],nrow(existingDF)),]
  existingDF[r[i],] <- newrow[i]
  existingDF
}


estimateRealcount<- function( observe, background ){
  if( nrow( as.matrix(observe) )!= nrow( as.matrix(background) ) ){
    MissedRecord <- which( rownames( observe) %in% setdiff( rownames( observe ),rownames( background ) ) )
    background <- insertRow( background, repmat(apply(background,2,mean),length(MissedRecord),1) ,MissedRecord)
  }
  signal <- observe - background
  signal[signal<0] <- 0
  signal
}


getsignalVars<-function( counts, bgcounts){
  bgsizeFactors <- estimateSizeFactorsForMatrixXBSeq(bgcounts)
  lambda <- rowMeans( sweep(bgcounts, 2, bgsizeFactors, '/' ) )
  observe_param <- getCountParams(counts, estimateSizeFactorsForMatrixXBSeq(counts))
  observe_sf <- estimateSizeFactorsForMatrixXBSeq(counts)
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
  signalvar
}

# Note: The following function is never called; it is here only for
# documentation purposes, as it has been used to produce the data object
# scvBiasCorrectionFits, which is stored in the file
# inst/scvBiasCorrectionFits.rda, gets loadewd by the line after this
# function and is used by the function adjustScvForBias
#
# To do: the correct place for this type of documentation is a vignette
#

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
  if(exists(("scvBiasCorrectionFits.rda")))
    load( "scvBiasCorrectionFits.rda" )
  else{
    scvBiasCorrectionFits<-prepareScvBiasCorrectionFits(maxnrepl=15, mu=100000, ngenes=10000,
                                 true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ))
    save(scvBiasCorrectionFits,file='scvBiasCorrectionFits.rda')
  }
  if( nsamples - 1 > length( scvBiasCorrectionFits ) )
    scv
  else
    ifelse( scv > .02,
            pmax( predict_helper( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv ),
            scv )   # For scv < .02, our fit is too coarse, but no correction seems necessary anyway
}


MAplot <- function(stats, ylim, padj=T, pcuff=0.1, lfccuff=1, linecol='red3',
                   xlab='mean of normalized counts', ylab=expression(log[2]~fold~change))
{
  if(!(is.data.frame(stats) && all(c("baseMean", "log2FoldChange") %in% colnames(stats))))
    stop("'stats' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  if(padj)
    col = ifelse(stats$padj>=pcuff, "gray32", "red")
  else
    col = ifelse(stats$pval>=pcuff, "gray32", "red")
  col = col[ stats$baseMean != 0 ]
  y = stats$log2FoldChange[ stats$baseMean != 0 ]
  stats = subset(stats, baseMean!=0)
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(y[is.finite(y)]), probs=0.99) * 1.1
  shape = ifelse(y<ylim[1], 6, ifelse(y>ylim[2], 2, 16) )
  stats$log2FoldChange = pmax( ylim[1], pmin(ylim[2], y) )

  ggplot() + geom_point( data=stats,aes( x=baseMean, y=log2FoldChange ), color=col, shape=shape ) +
    ylim(ylim) + geom_hline(yintercept=0,colour=linecol,size=1)  + scale_x_log10() +
    labs( x=xlab, y=ylab )
}


plotSCVEsts = function( XB, name=NULL, ymin, linecol='red3',
                        xlab = "mean of normalized counts", ylab = "SCV", ... )
{
  px = rowMeans( counts( XB, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]
  py = fitInfo(XB, name=name)$perGeneSCVEsts[sel]
  if(missing(ymin))
    ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
  py = pmax(py, ymin)
  shape = ifelse(py<ymin, 6, 16)
  sel1 = complete.cases(px)&complete.cases(py)
  px = px[sel1]
  py = py[sel1]
  shape = shape[sel1]
  fitd = data.frame(px=px,py=py)
  fx = 10^seq( -.5, 5, length.out=100 )
  fy = fitInfo(XB, name=name)$SCVFunc(fx)
  fitl = data.frame(fx=fx, fy=fy)
  ggplot() + geom_point( data=fitd, aes( x=px, y=py), shape=shape) +
    geom_line( data=fitl, aes ( x=fx, y=fy), shape=shape, col=linecol, size=1.5, alpha=0.6) + scale_x_log10() +
    scale_y_log10() + labs(x=xlab, y=ylab)
}

setClass( "XBSeqDataSet",
   contains = "eSet",
   representation = representation(
      fitInfo = "environment",
      dispTable = "character",
      multivariateConditions = "logical" ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), XBSeqDataSet = "1.1.0" ) ) )
)

setValidity( "XBSeqDataSet", function( object ) {
  if( length( object@multivariateConditions ) != 1 )
    return( "multivariateConditions is not scalar." )
  if( ! "sizeFactor"  %in% names(pData(object)) )
    return( "phenoData does not contain a 'sizeFactor' columns.")
  if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
    return( "The 'sizeFactor' column in phenoData is not numeric." )
  if( ! object@multivariateConditions ) {
    if( ! "condition"  %in% names(pData(object)) )
      return( "phenoData does not contain a 'condition' columns." )
    if( ! is( pData(object)$`condition`, "factor" ) )
      return( "The 'condition' column in phenoData is not a factor." )
  }
  if( !is.integer( counts(object) ) )
    return( "the count data is not in integer mode" )
  if( any( counts(object) < 0 ) )
    return( "the count data contains negative values" )
  TRUE
} )

setMethod("counts", signature(object="XBSeqDataSet"),
          function(object, normalized=FALSE ) {
            if(!normalized) {
              assayData(object)[["counts"]]
            } else {
              if(any(is.na( sizeFactors(object)))) {
                stop( "Please first calculate size factors or set normalized=FALSE")
              } else {
                t(t( assayData(object)[["counts"]] ) / sizeFactors(object) )
              }
            }
          })

setReplaceMethod("counts", signature(object="XBSeqDataSet", value="matrix"),
                 function( object, value ) {
                   assayData(object)[[ "counts" ]] <- value
                   validObject(object)
                   object
                 })

setMethod("sizeFactors", signature(object="XBSeqDataSet"),
          function(object) {
            sf <- pData(object)$sizeFactor
            names( sf ) <- colnames( counts(object) )
            sf
          })

setReplaceMethod("sizeFactors", signature(object="XBSeqDataSet", value="numeric"),
                 function( object, value ) {
                   pData(object)$sizeFactor <- value
                   validObject( object )
                   object
                 })

setMethod("conditions", signature(object="XBSeqDataSet"),
          function( object, ... ) {
            if(length(list(...))!=0)
              warning("in conditions: Ignoring second and/or further arguments.")
            if( object@multivariateConditions && !("condition" %in% colnames(pData(object))))
              stop( "Could not find 'condition' column in pData." )
            conds <- pData(object)$`condition`
            names( conds ) <- colnames( counts(object) )
            conds
          })

setReplaceMethod("conditions", signature(object="XBSeqDataSet"),
                 function( object, value ) {
                   if( object@multivariateConditions )
                     stop( "The 'conditions<-' accessor is only for simple single-factor conditions, but you have specified multivariate conditions. Access them via 'pData<-'." )
                   pData(object)$`condition` <- factor( value )
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



newXBSeqDataSet <- function( countData, conditions, sizeFactors=NULL,
      phenoData = NULL, featureData = NULL )
{
   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) )
      stop( "The input data has to be integer!" )
   mode( countData ) <- "integer"
   if( is.null( sizeFactors ) ) {
      sizeFactors <- rep( NA_real_, ncol(countData) )
   } else
      warning( "The 'sizeFactor' argument is deprecated. Use 'estimateSizeFactors'." )
   if( is.null( phenoData ) )
      phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   if( is.null( featureData ) )
      featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )

   phenoData$`sizeFactor` <- sizeFactors
   varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <-
      "size factor (relative estimate of sequencing depth)"

   if( is( conditions, "matrix" ) )
      conditions <- as.data.frame( conditions )

   if( is( conditions, "data.frame" ) || is( conditions, "AnnotatedDataFrame" ) ) {
      stopifnot( nrow( conditions ) == ncol( countData ) )
      conditions <- as( conditions, "AnnotatedDataFrame" )
      dimLabels( conditions ) <- dimLabels( phenoData )
      rownames( pData(conditions) ) <- rownames( pData(phenoData) )
         # TODO: What if the rownames were set?
      phenoData <- combine( phenoData, conditions )
      multivariateConditions <- TRUE
      rvft <- c( `_all` = NA_character_ )
   } else {
      conditions <- factor( conditions )
      stopifnot( length( conditions ) == ncol( countData ) )
      phenoData$`condition` <- factor( conditions )
      varMetadata( phenoData )[ "condition", "labelDescription" ] <-
         "experimental condition, treatment or phenotype"
      multivariateConditions <- FALSE
      rvft <- rep( NA_character_, length(levels(conditions)) )
   }

   XB <- new( "XBSeqDataSet",
      assayData = assayDataNew( "environment", counts=countData ),
      phenoData = phenoData,
      featureData = featureData,
      multivariateConditions = multivariateConditions,
      fitInfo = new.env( hash=TRUE ),
      dispTable = rvft )

   XB
}

