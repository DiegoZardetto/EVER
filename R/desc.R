`desc` <-
function (deskott, descfun = NULL, ...)
#######################################################################
#  Descrive sinteticamente un oggetto di classe kott.design.          #
#  NOTA: Attraverso 'descfun' e' possibile specificare quale funzione #
#        utilizzare per descrivere i dati del campione.               #
#        La funzione 'descfun' deve poter agire su oggetti di classe  #
#        dataframe (possibili esempi: head, str, summary, ...).       #
#  NOTA: L'argomento '...' serve a passare parametri opzionali alla   #
#        funzione 'descfun'.                                          #
#######################################################################
{
    if (!inherits(deskott, "kott.design")) 
        stop("Object ", substitute(deskott), " must be of class kott.design")
    if (inherits(deskott, "kott.cal.design")) {
        cal <- "Calibrated"
    }
    else {
        cal <- NULL    
    }
    cat(cal,"DAGJK replicated survey data\n")
    cat(" >",attr(deskott, "nrg"),"random groups\n")
    cat(" > ")
    if (attr(deskott, "strata") != FALSE) {
        cat("Stratified ")
        Strata <- names(model.frame(attr(deskott, "strata"), deskott[1, ]))
        nStrata <- length(unique(deskott[, Strata])) 
    }
    ids.charvect <- names(model.frame(attr(deskott, "ids"), deskott[1, ]))
    nstadi <- length(ids.charvect)
    id.PSU <- ids.charvect[1]
    nPSU <- length(unique(deskott[, id.PSU]))
    nFSU <- nrow(deskott)
    if (nPSU == nFSU) {
        cat("Independent Unit Sampling design")
        if (attr(deskott, "strata") != FALSE) 
            cat(paste("\n   - [", nStrata , "] strata", sep = ""))
        cat("\n")
    }
    else {
        cat(nstadi, "- Stage Cluster Sampling design")
        if (attr(deskott, "strata") != FALSE) 
            cat(paste("\n   - [", nStrata , "] strata", sep = ""))
        nclus <- sapply(ids.charvect, function(var) length(unique(deskott[, 
            var])))
        cat(paste("\n   - [", paste(nclus, collapse = ","), 
            "] clusters", sep = ""))
        if (!identical(attr(deskott, "self.rep.str"), FALSE)) {
            cat(paste("\n   - [", attr(deskott, "nvar.PSU"), "] variance PSUs", sep = ""))
        }
        cat("\n")
    }
    cat("\nCall:\n")
    print(attr(deskott, "call"))
    cat("\n")
    if (!is.null(descfun)) {
        cat("\n*****************************************\n")
        cat("               survey data                   ")
        cat("\n*****************************************\n")
        descfun(deskott, ...)
    }
}