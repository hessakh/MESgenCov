#' Instance of dataframe of residuals for 50 sites
#' Sites included:1] "WV18" "AK03" "CA75" "PR20" "ND11" "TX03" "ME00" "FL11" "MT05" "CO00" "MO03"
#' "WA14" "MI99" "AL10" "MT00" "NY20" "OK17" "AZ03" "CA76" "NE99" "IN20" "ID11"
#' "TX22" "SC06" "LA12" "UT01" "MN16" "MA01" "TN00" "AR03" "OK29" "NE15" "IA08"
#' "MI09" "FL03" "WY00" "NC03" "AZ99" "IL63" "MT98" "NJ99" "WA24" "NV05" "TX16"
#'  "MN27" "WY99" "NY10" "WY08" "MS30" "OR97"

#' For dates: 01/01/86 00:00 - 12/31/94 00:00
#' For compound: SO4
#' Missing values filled in use fitted LambertW distribution
#' Site list produced by maxDistSites("01/01/86 00:00","12/31/94 00:00",50,200,"SO4",1)
#' 
#' @docType data
#'
#' @usage data(dfRes50)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#'
#' @examples
#' 
#' data(dfRes50)
#' dfRes <- dfRes50
#' dfRes

"dfRes50"

