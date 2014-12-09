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
                    fitType = c("local", "parametric"),
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
XBSeq <- function( observe, background, conditions, method='pooled', sharingMode='maximum', fitType='local', pvals_only=FALSE ){
  Signal <- estimateRealcount(observe, background)
  XB <- newXBSeqDataSet(Signal,conditions)
  XB <- estimateSizeFactorsXBSeq( XB )
  XB <-estimateSCV( XB, observe, background, method=method, sharingMode=sharingMode, fitType=fitType )
  Teststas <- XBSeqTest( XB, levels(conditions)[1L], levels(conditions)[2L], pvals_only=pvals_only )
  Teststas
}