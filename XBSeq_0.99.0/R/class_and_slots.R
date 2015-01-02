require('pracma')
require('matrixStats')
require('locfit')
require('Biobase')
require('ggplot2')
require('methods') 

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

