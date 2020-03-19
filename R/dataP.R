#' Daily NADP precipitation data
#'
#' Daily NADP data of precipitation depth in inches^3 for all sites from 1979 to 2017.
#'
#' @docType data
#'
#' @usage data(preDaily)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references National Atmospheric Deposition Program (NRSP-3). 2019. NADP Program Office, Wisconsin State Laboratory of Hygiene,
#' 465 Henry Mall, Madison, WI 53706.
#' (\href{http://nadp.slh.wisc.edu/data/NTN/ntnAllsites.aspx}{NADP/NTN})
#'
# @source \href{http://nadp.slh.wisc.edu/data/NTN/ntnAllsites.aspx}{NADP/NTN}
#'
#' @examples
#' data(preDaily)
#'  amount   <- attr(preDaily, "amount")
#' siteID <- preDaily$siteID

"preDaily"
