
# You should not need to call these directly                    


list2df <- function(X){
    data.frame(t(sapply(X,function(x) unlist(x))))
}

getDataFromWS <- function(cmd,Q){
    svc <- paste('http://httr-dev.epa.gov/api/httr/v1/',cmd,sep="")
    
    out <- tryCatch(
        {httr::GET(
                     url=svc,
                     query=Q
                    ) -> resp
        httr::content(resp)
        },
        error=function(cond) {
            message(paste("Query could not be exectuted: check params"))
            message("Here's the error message:")
            message(cond)
            return(NA)
        }
        )
    out
}
                        
