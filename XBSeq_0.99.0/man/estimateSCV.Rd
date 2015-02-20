\name{estimateSCV}
\docType{methods}
\alias{estimateSCV}
\alias{estimateSCV,XBSeqDataSet-method}
\title{Estimate squared coefficient of variation for each gene}
\description{
   A similar method is applied to estimate the SCV for each gene based on the method used in DESeq
}
\usage{
\S4method{estimateSCV}{XBSeqDataSet}( object, observe, background,
   method = c( "pooled", "pooled-CR", "per-condition", "blind" ),
   sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
   fitType = c("local","parametric"),
   locfit_extra_args=list(), lp_extra_args=list(),
   modelFrame = NULL, modelFormula = count ~ condition, ... )
}
\arguments{
  \item{object}{
      a \code{XBSeqDataSet} with size factors.
  }
  \item{observe}{
  The observed read count from exonic regions
  }
  \item{background}{
  The background noise from non-exonic regions
  }
  \item{method}{
     There are three ways how the empirical dispersion can be computed:

     \itemize{
     \item \code{pooled} - Use the samples from all conditions with
     replicates to estimate a single pooled empirical dispersion value,
     called "pooled", and assign it to all samples.

     \item \code{pooled-CR} - Fit models according to \code{modelFormula} and
     estimate the dispersion by maximizing a Cox-Reid adjusted profile
     likelihood (CR-APL). This method
     is much slower than \code{method=="pooled"} but works also with crossed
     factors (as may occur, e.g., in designs with paired samples).
     Usually, you will need to specify the model formula, which should be the
     same as the one used later in the call to \code{nbinomFitGLMs} for fitting the full model.

     Note: The method of using CR-APL maximization for this application has been
     developed by McCarthy, Chen and Smyth [Nucl. Acid Res., 2012 and been
     first implemented in edgeR (in 2010). DESeq optimizes the expression
     for the CR-APL given in McCarthy et al.'s paper, but does not use the
     weigthed maximum likelihood scheme proposed there.

     \item \code{per-condition} - For each condition with replicates, compute
     a gene's empirical dispersion value by considering the data from samples for this
     condition.  For samples of unreplicated conditions, the maximum
     of empirical dispersion values from the other conditions is used.  If \code{object} has a
     multivariate design (i.e., if a data frame was passed instead of a
     factor for the \code{condition} argument in
     \code{\link{newXBSeqDataSet}}), this method is not available. (Note: This method was called ``normal'' in
        previous versions.)

     \item \code{blind} - Ignore the sample labels and compute a
     gene's empirical dispersion value as if all samples were replicates of a
     single condition. This can be done even if there are no biological
     replicates. This method can lead to loss of power; see
     the vignette for details. The single estimated dispersion condition
     is called "blind" and used for all samples.

     } }

   \item{sharingMode}{
     After the empirical dispersion values have been computed for each
     gene, a dispersion-mean relationship is fitted for sharing
     information across genes in order to reduce variability of the
     dispersion estimates. After that, for each gene, we have two values: the
     empirical value (derived only from this gene's data), and the
     fitted value (i.e., the dispersion value typical for genes with an
     average expression similar to those of this gene). The
     \code{sharingMode} argument specifies which of these two values
     will be written to the \code{featureData}'s \code{disp_} columns
     and hence will be used by the functions \code{\link{XBSeqTest}}

     \itemize{
     \item \code{fit-only} - use only the fitted value, i.e., the
     empirical value is used only as input to the fitting, and then
     ignored. Use this only with very \emph{few} replicates, and when
     you are not too concerned about false positives from dispersion outliers, i.e. genes
     with an unusually high variability.

     \item \code{maximum} - take the maximum of the two values. This is
     the conservative or prudent choice, recommended once you have at
     least three or four replicates and maybe even with only two replicates.

     \item \code{gene-est-only} - No fitting or sharing, use only the
     empirical value. This method is preferable when the number of
     replicates is large and the empirical dispersion values are
     sufficiently reliable. If the number of replicates is small, this
     option may lead to many cases where the dispersion
     of a gene is accidentally underestimated and a false positive arises in
     the subsequent testing.
     }
  }

  \item{fitType}{

     \itemize{
     \item \code{parametric} - Fit a dispersion-mean relation of the
     form \code{dispersion = asymptDisp + extraPois / mean} via a robust
     gamma-family GLM. The coefficients \code{asymptDisp} and \code{extraPois}
     are given in the attribute \code{coefficients} of the \code{dispFunc}
     in the \code{fitInfo} (see below).
     }

     \itemize{
     \item \code{local} - Use the locfit package to fit a dispersion-mean
     relation, as described in the DESeq paper.
     }

  }

  \item{locfit_extra_args, lp_extra_args}{
      (only for \code{fitType=local})
      Options to be passed to the \code{locfit} and to the \code{lp}
      function of the locfit package. Use this to adjust the local
      fitting. For example, you may pass a value for \code{nn} different
      from the default (0.7) if the fit seems too smooth or too rough by
      setting \code{lp_extra_agrs=list(nn=0.9)}. As another example, you
      can set \code{locfit_extra_args=list(maxk=200)} if you get the
      error that locfit ran out of nodes. See the documentation of the
      \code{locfit} package for details. In most cases, you will not
      need to provide these parameters, as the defaults seem to work
      quite well.}

   \item{modelFrame}{
      By default, the information in \code{conditions(object)} or \code{pData(object)} is used to
      determine which samples are replicates (see \code{\link{newXBSeqDataSet}}).
      For \code{method="pooled"}, a data frame can be passed
      here, and all rows that are identical in this data frame are considered
      to indicate replicate samples in \code{object}. For \code{method="pooled-CR"},
      the data frame is used in the fits. For the other methods, this argument is
      ignored.
   }

   \item{modelFormula}{
      For \code{method="pooled-CR"}, this is the formual used for the dispersion fits.
      For all other methods, this argument is ignored.
   }

   \item{...}{extra arguments are ignored}

}
\details{
   The details regarding which option to choose can be found in the DESeq help page. Generally
   speaking, if you have less number of replicates (<=3), set \code{method="pooled"}. Otherwise,
   try \code{method="per-condition"}. We revised the code to estimate the variance of the true
   signal by using variance sum law rather than calculate the variance directly.
}
\value{
   The \code{XBSeqDataSet} cds, with the slots \code{fitInfo} and
   \code{featureData} updated.
}
\author{
   Yuanhang Liu
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\examples{
  set.seed(1990)
   conditions <- factor(c('C1','C1','C2','C2'))
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   Signal <- estimateRealcount(observe_signal, background_noise)
   XB <- newXBSeqDataSet(Signal,conditions)
   XB <- estimateSizeFactorsXBSeq( XB )
   XB <- estimateSCV( XB, observe_signal, background_noise,fitType='local')
   str( fitInfo( XB ) )
}
