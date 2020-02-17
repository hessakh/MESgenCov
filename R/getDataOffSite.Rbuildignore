#' @keywords internal

### Downloads weekly data and precipitation data using site links if they aren't already in the working folder
getDataOffSite <- function(){
  #----------downloads weekly data-----------
  if(!exists("weeklyCSV")){
    if(!file.exists("NTN-All-w.csv")){
      url      <- "http://nadp.slh.wisc.edu/datalib/ntn/weekly/NTN-All-w.csv"
      folder   <- getwd()
      downloadCSV("NTN-All-w.csv",url,folder)}
    else {
      weeklyCSV <- read.csv("NTN-All-w.csv", stringsAsFactors = F)
    }
      }else{data("weeklyCSV")}
  #----------------------------------------------

  #----------downloads precipitation data and combines it into one file-------------
  ###make sure you're in package
  if(!exists("preDailyCSV")){
    sites     <- unique(toupper(weeklyCSV$siteID))
    if(file.exists("dailydataReq.csv")){
      preCSV0           <- read.csv("dailydataReq.csv", stringAsFactors = F)
      colnames(preCSV0) <- c("siteID","labno","starttime","endtime","amount","type")
      preCSV0$amount    <- as.numeric(preCSV0$amount)
    }else if(file.exists("NTN-All-p-mod.csv")){
      preCSV0           <- read.csv("NTN-All-p-mod.csv")
      sitesp            <- unique(preCSV0$siteID)
      sites             <- setdiff(sites,sitesp)
    }else{preCSV0 <- NULL}

    #put precipitation data in pdata file
    if(!file.exists("pdata")){
    	dir.create("pdata",showWarnings = FALSE)
    }
    setwd("pdata")
    folder <- getwd()
    burli  <- "http://nadp.slh.wisc.edu/datalib/ntn/daily/NTN-"
     #download site precipitation data
    	for(i in seq_along(sites)){
    		tryCatch(
    			{sitei <- sites[i]
    			urlpi  <- paste0(burli,sitei,"-d.csv")
    			namepi <- paste0("NTN-",sitei,"-d.csv")
    			if(!file.exists(namepi)){
    			  downloadCSV(namepi,urlpi,folder)
    			}
    			},
    			warning = function(war) paste0("Missing site: ",sites[i], ", try url manually: ", urlpi))
    	}

    #combine data ##maybe not all i
    allp <- combinePrecipitation(sites,preCSV0,"NTN-All-p.csv")
    preDailyCSV <- allp
    setwd("..")
  }else{data("preDailyCSV")}
#-------------------------------------------------------------------------------
}
