download.NDVI <- function(URL) {

  # Wrap function in an if/else loop that checks if the URL is valid
  if (length(URL) == 1 & is.character(URL) & substr(URL,1,4)=="http") {

    # Read in data
    modat <- read.csv(URL)

    # Convert Digital Numbers (more efficient for data storage) to NDVI
    modat$NDVI <- modat$value*0.0001

    # Convert calendar_date to class "Date"
    modat$calendar_date <- as.Date(as.character(modat$calendar_date))

    # Return the data
    return(modat)

  } else {

    # If the URL is not valid return...
    print(paste("download.NDVI: Input URL not provided correctly",URL))
  }
}