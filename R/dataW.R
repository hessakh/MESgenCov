#' Weekly NADP data
#'
#' Weekly NADP data of 10 deposits measured in mg/L for all sites from 1978 to January 2019.
#'
#' @docType data
#'
#' @usage data(weeklyConc)
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
#' data(weeklyConc)
#'  SO4   <- attr(weeklyConc, "SO4")
#' siteID <- weeklyConc$siteID

"weeklyConc"