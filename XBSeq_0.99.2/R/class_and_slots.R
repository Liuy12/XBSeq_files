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
} )


setMethod("conditions", signature(object="XBSeqDataSet"),
          function( object, ... ) {
            if(length(list(...))!=0)
              warning("in conditions: Ignoring second and/or further arguments.")
            conds <- object@conditions
            names( conds ) <- colnames( counts(object) )
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


XBSeqDataSet <- function( countData, conditions, sizeFactors=NULL, ...)
{
  countData <- as.matrix( countData )
  conditions <- as.factor(conditions)
  if( any( round( countData ) != countData ) )
    stop( "The input data has to be integer!" )
  mode( countData ) <- "integer"
  if( is.null( sizeFactors ) ) {
    sizeFactors <- rep( NA_real_, ncol(countData) )
  } 
  colData <- data.frame(condition = factor( conditions ))
  rownames(colData) <- colnames(countData)
  stopifnot( length( conditions ) == ncol( countData ) )
  rvft <- rep( NA_character_, length(levels(conditions)) )
  XB <- DESeqDataSetFromMatrix(countData, colData, ~condition )
  XB <- new( "XBSeqDataSet",
             XB, 
             fitInfo = new.env( hash=TRUE ),
             dispTable = rvft,
             conditions = conditions)
  return(XB)
}
