#' NADP monitor GPS location for all sites
#'
#' Location data for all sites currently or formerly active in the NTN network.
#'
#' @docType data
#'
#' @usage data(NADPgeo)
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
#' data(NADPgeo)
#'  Latitude   <- attr(NADPgeo, "Latitude (DD)")
#' siteID <- NADPgeo$`Site ID`

"NADPgeo"

