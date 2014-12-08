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
  kAs <- apply( countsA,1,sum )
  kBs <- apply( countsB ,1,sum)

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
  if(file.exists(("scvBiasCorrectionFits.rda")))
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