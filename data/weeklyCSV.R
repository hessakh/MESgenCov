#' Weekly NADP data
#'
#' Weekly NADP data of 10 deposits measured in mg/L for all sites from 1978 to January 2019.
#'
#' @docType data
#'
#' @usage data(weeklyCSV)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references National Atmospheric Deposition Program (NRSP-3). 2019. NADP Program Office, Wisconsin State Laboratory of Hygiene,
#' 465 Henry Mall, Madison, WI 53706.
#' (\href{http://nadp.slh.wisc.edu/data/NTN/ntnAllsites.aspx}{NADP/NTN})
#'
#' @source \href{http://nadp.slh.wisc.edu/data/NTN/ntnAllsites.aspx}{NADP/NTN}
#'
#' @examples
#' data(weeklyCSV)
#'  SO4   <- attr(weeklyCSV, "SO4")
#' siteID <- weeklyCSV$siteID

"weeklyCSV"
