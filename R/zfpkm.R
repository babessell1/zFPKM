#AUTHOR: Ron Ammar <ron.ammar@bms.com>, John R Thompson <john.thompson@bms.com>
#COMPANY: Bristol-Myers Squibb Co.
#DATE: 2015-07-08
#OBJECTIVE: Perform the zFPKM transform on RNA-seq FPKM data.
#SUMMARY:
#  Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is based on
#  the publication by Hart et al., 2013 (Pubmed ID 24215113).
#CHANGELOG:
# JRT 8Sep2015:
#   Modified zFPKMTransformDF to optionally allow plotting and saving the
#   plot to a PNG file and set an output path.
#   Modified PlotGaussianFitDF to move legend to the top and optionally remove facet
#   titles on each individual plot (to allow more data to be plotted and still visually
#   interpretable).
#    Modfied to use double colon references to external libraries
# JRT 4Mar2016:
#   Added "floor" parameter to set lower limit on x axis log2FPKM value
# RA 10May2017:
#   Updated plotting options. Minor refactor. Using checkmate.
# RA 1Jun2017:
#   Separating zFPKM calc and plotting as per Bioconductor review suggestion.
# RA 10Jul2017:
#   Style changes for Bioconductor submission.
# BAB 10Mar2022:
#   Rolling average to filter insignificant local maxima and implementation of fitting
#   Gaussian to highest local maxima (with respect to X) rather than global maxima


#' zFPKM Transformation
#'
#' Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is
#' based on the publication by Hart et al., 2013 (Pubmed ID 24215113). Reference recommends
#' using zFPKM > -3 to select expressed genes.  Validated with encode open/closed promoter chromatin structure epigenetic
#' data on six of the ENCODE cell lines.  Works well for gene level data using FPKM or TPM. Does not appear to calibrate well
#' for transcript level data.
#'
#' @author Ron Ammar, \email{ron.ammar@bms.com}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/24215113}
#' @keywords zFPKM
#'
#' @param fpkmDF A SummarizedExperiment or data frame containing raw FPKM (or TPM)
#'  values. Each row corresponds to a gene/transcript and each column corresponds
#'  to a sample.
#'  NOTE: these are NOT log_2 transformed. Also, the rownames are gene/transcript
#'  names and NOT included as a separate column
#' @param assayName When input is a SummarizedExperiment, names the specific
#'  assay. Typically one of "fpkm" or "tpm" [default = "fpkm"]
#'
#' @return zFPKM data frame
#'
#' @examples
#' library(dplyr)
#' gse94802 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_normalized_FPKM.csv.gz"
#' temp <- tempfile()
#' download.file(gse94802, temp)
#' fpkm <- read.csv(gzfile(temp), row.names=1)
#' MyFPKMdf <- select(fpkm, -MGI_Symbol)
#'
#' zfpkm <- zFPKM(MyFPKMdf)
#'
#' @import checkmate dplyr ggplot2 tidyr SummarizedExperiment
#'
#' @export
zFPKM<- function(fpkmDF, min_thresh, assayName="fpkm") {

  assert(checkDataFrame(fpkmDF), checkClass(fpkmDF, "SummarizedExperiment"),
         combine="or")

  return(zFPKMTransform(fpkmDF, min_thresh, assayName)[[2]])
}


removeNanInfRows <- function(fpkm) {
  # Remove FPKM rows containing all NaN values. These are most likely a result
  # of effective lengths = 0 when calculating FPKM.
  return(fpkm[which(!apply(fpkm, 1, function(r) all(is.nan(r) | is.infinite(r)))), ])
}


#' zFPKM Transformation
#'
#' Perform the zFPKM transform on RNA-seq FPKM data. This algorithm is
#' based on the publication by Hart et al., 2013 (Pubmed ID 24215113). Reference recommends
#' using zFPKM > -3 to select expressed genes.  Validated with encode open/closed promoter chromatin structure epigenetic
#' data on six of the ENCODE cell lines.  Works well for gene level data using FPKM or TPM. Does not appear to calibrate well
#' for transcript level data.
#'
#' @author Ron Ammar, \email{ron.ammar@bms.com}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/24215113}
#' @keywords zFPKM
#'
#' @param fpkmDF A SummarizedExperiment or data frame containing raw FPKM (or TPM)
#'  values. Each row corresponds to a gene/transcript and each column corresponds
#'  to a sample.
#'  NOTE: these are NOT log_2 transformed. Also, the rownames are gene/transcript
#'  names and NOT included as a separate column
#' @param assayName When input is a SummarizedExperiment, names the specific
#'  assay. Typically one of "fpkm" or "tpm" [default = "fpkm"]
#' @param FacetTitles use to label each facet with the sample name [default = FALSE]
#' @param PlotXfloor Lower limit for X axis (log2FPKM units) [default = -20] set to NULL to disable
#'
#' @return Displays plots of zFPKM distributions
#'
#' @examples
#' library(dplyr)
#' gse94802 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_normalized_FPKM.csv.gz"
#' temp <- tempfile()
#' download.file(gse94802, temp)
#' fpkm <- read.csv(gzfile(temp), row.names=1)
#' MyFPKMdf <- select(fpkm, -MGI_Symbol)
#'
#' zFPKMPlot(MyFPKMdf)
#'
#' @import checkmate dplyr ggplot2 tidyr SummarizedExperiment zoo
#'
#' @export
zFPKMPlot <- function(fpkmDF, min_thresh, assayName="fpkm", FacetTitles=FALSE, PlotXfloor=-20) {

  assert(checkDataFrame(fpkmDF), checkClass(fpkmDF, "SummarizedExperiment"),
         combine="or")
  PlotGaussianFitDF(zFPKMTransform(fpkmDF, min_thresh, assayName)[[1]], FacetTitles, PlotXfloor)
}


zFPKMTransform <- function(fpkmDF, min_thresh, assayName) {
  # Helper function for zFPKM output and plotting. Do not call directly.
  #
  # Args:
  #   Defined by zFPKM() and zFPKMPlot()
  #
  # Returns:
  #   Internal objects used for zFPKM calculations and plotting.

  if (is(fpkmDF, "SummarizedExperiment"))  {
    fpkmDF <- assay(fpkmDF, assayName)
  }
  fpkmDF <- removeNanInfRows(fpkmDF)
  zFPKMDF <- data.frame(row.names=row.names(fpkmDF))
  outputs <- list()
  for (c in colnames(fpkmDF)) {
    output <- zFPKMCalc(fpkmDF[, c], min_thresh)
    zFPKMDF[, c] <- output[["z"]]
    outputs[[c]] <- output
  }

  return(list(outputs, zFPKMDF))
}


zFPKMCalc <- function(fpkm, min_thresh) {
  # Performs the zFPKM transform on RNA-seq FPKM data. This involves fitting a
  # Gaussian distribution based on the right side of the FPKM distribution, as
  # described by Hart et al., 2013 (Pubmed ID 24215113). The zFPKM transformed
  # FPKM values represent normalized FPKM data, and, according to Hart et al.,
  # a zFPKM value >= -3 can be considered expressed (has an active promoter).
  #
  # Args:
  #   fpkm: a vector of raw FPKM values. NOTE: these are NOT log_2 transformed.
  #
  # Returns:
  #   The zFPKM transformed vector of the input FPKM data

  if (!is.numeric(fpkm)) {
    stop("argument 'fpkm' must be numeric")
  }

  # log_2 transform the FPKM values
  fpkmLog2_filt <- log(fpkm[fpkm>min_thresh], base=2)
  fpkmLog2 <- log(fpkm, base=2)
  #print(fpkmLog2[1:100])

  # Compute kernel density estimate
  d <- density(fpkmLog2_filt)


  # calculate rolling average
  perc <- as.integer(0.1*length(d[["y"]]) + 1) # 10% roll avg interval

  d[["roll_y"]] <- zoo::rollmean(d[["y"]], perc)

  # from https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
  find_maxima <- function (x, m = 1){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
  }

  local_maxes <- find_maxima(d[["roll_y"]][floor(0.3*length(d[["roll_y"]])):0.7*length(d[["roll_y"]])])
  fit_max <- max(local_maxes) + as.integer(perc/2)
  

  # Set the maximum point in the density as the mean for the fitted Gaussian
  mu <- d[["x"]][fit_max] # get max with respect to x) local maxima of rolling
  max_y <- d[["y"]][fit_max]
  cnt <- 0
  
  while  ( (max_y < 0.1*max(d[["y"]])) && (cnt < 20) && FALSE) { # while selected local max y is less than 20% of actual maximum 
	cnt <- cnt + 1
    perc <- as.integer((0.2-(cnt*0.01))*length(d[["y"]]) + 1) # rm 1 percent from roll avg interval per iteration

    #d[["roll_y"]] <- filter(data.frame(d[["y"]]), f_2perc, sides=2)
    d[["roll_y"]] <- zoo::rollmean(d[["y"]], perc)

	local_maxes <- local_maxes[floor(0.3*length(local_maxes)):0.8*length(local_maxes)]
    fit_max <- max(local_maxes) + as.integer(perc/2)
    # Set the maximum point in the density as the mean for the fitted Gaussian
    #mu <- d[["x"]][which.max(d[["y"]])]
    mu <- d[["x"]][fit_max] # get max with respect to x) local maxima of rolling
    max_y <- d[["y"]][fit_max]
	
	#if ( which.max(local_maxes) < which(local_maxes == max_y)
	
  }

  if ( FALSE && (max_y < 0.1*max(d[["y"]])) ) {  
    mu <- d[["x"]][which.max(d[["y"]])]
    max_y <- max(d[["y"]]) # if doesnt work use regular zFPKM calculation
	# TODO: FAILURE MESSAGE AND OUTPUT TO LIST?
  }

  if ( TRUE ) {
	#local_maxes <- local_maxes[floor(0.3*length(local_maxes)):length(local_maxes)]
	
	local_maxes <- sort(local_maxes, decreasing=TRUE)
	while ( length(local_maxes) > 1 && d[["x"]][local_maxes[2]] > d[["x"]][local_maxes[1]] ) {	
		local_maxes <- local_maxes[-1]
	}
	fit_max <- local_maxes[1]
	mu <- d[["x"]][fit_max] # get max with respect to x) local maxima of rolling
	max_y <- d[["y"]][fit_max]
  }

  # Determine the standard deviation
  U <- mean(fpkmLog2[fpkmLog2 > mu])
  stdev <- (U - mu) * sqrt(pi / 2)

  # Compute zFPKM transform
  zFPKM <- (fpkmLog2 - mu) / stdev

  result <- ZFPKMResult(zFPKM, d, mu, stdev, max_y)

  return(result)
}


ZFPKMResult <- function(zfpkmVector, density, mu, stdev, max_y) {
  # S3 class to store zFPKM vector and related metrics to be used in plotting

  zfpkmRes <- list(
    z = zfpkmVector,
    d = density,
    m = mu,
    s = stdev,
    max_y = max_y
  )

  class(zfpkmRes) <- append(class(zfpkmRes), "zFPKM")
  return(zfpkmRes)
}


PlotGaussianFitDF <- function(results, FacetTitles=TRUE, PlotXfloor) {
  # Plot a grid of the log_2(FPKM) density and the fitted Gaussian (scale both
  # to the same density scale so that the curves overlap).

  megaDF <- data.frame()

  for (name in names(results)) {
    result <- results[[name]]
    d <- result[["d"]]
    mu <- result[["m"]]
    stdev <- result[["s"]]
    max_y <- result[["max_y"]]

    # Get max of each density and then compute the factor to multiply the
    # fitted Gaussian. NOTE: This is for plotting purposes only, not for zFPKM
    # transform.
    fitted <- dnorm(d[["x"]], mean=mu, sd=stdev)

    #maxFPKM <- max(d[["y"]])
    maxFPKM <- max_y
    maxFitted <- max(fitted)

    scaleFitted <- fitted * (maxFPKM / maxFitted)

    df <- data.frame(sample_name=name, log2fpkm=d[["x"]], fpkm_density=d[["y"]],
                     fitted_density_scaled=scaleFitted)


    megaDF <- megaDF %>% dplyr::bind_rows(df)
  }

  megaDFG <- megaDF %>% tidyr::gather(source, density, -c(log2fpkm, sample_name))
  labels <- unique(megaDFG$sample_name)

  maxX = max(megaDFG[["log2fpkm"]])
  maxY = max(d[["y"]])

  p <- ggplot2::ggplot(megaDFG, ggplot2::aes(x=log2fpkm, y=density, color=source)) +
    #ggplot2::facet_wrap(~ sample_name) +
    ggplot2::facet_wrap(vars(sample_name)) +
    ggplot2::geom_line(alpha=0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="log2(FPKM)", y="[scaled] density")  +
    ggplot2::theme(legend.position="top") +
    ggplot2::xlim(PlotXfloor, maxX)


  print(p)
}
