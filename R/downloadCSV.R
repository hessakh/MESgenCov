#' @keywords internal
### code from Automated Data Collection with R by C. Rubba, D. Nyhuis, S. Munzert

downloadCSV <- function(filename,fileurl, folder){
	#dir.create(folder,showWarnings = FALSE) #creates new folder
	if(!file.exists(paste0(folder, "/", filename))){
		download.file(fileurl, destfile = paste0(folder, "/", filename))
		Sys.sleep(1) #slows executions by a second
	}
}

