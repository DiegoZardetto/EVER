`kott.addvars` <-
function(deskott, ...)
##########################################################
#  Aggiorna i dati contenuti in un oggetto kott.design.  #
#  Gli argomenti ... devono essere del tipo tag = expr   #
#  (un tag puo' essere o un identificatore o una stringa #
#  di caratteri).                                        #
#  La nuova variabile creata avra' nome 'tag' e valori   #
#  ottenuti valutando 'expr' su deskott.                 #
#  NOTA: Eventuali espressioni sprovviste di tag in ...  #
#        vengono ignorate e non hanno, dunque, nessun    #
#        effetto sul valore di ritorno di kott.addvars.  #
##########################################################
{
    if (!inherits(deskott, "kott.design")) 
        stop("Object ", substitute(deskott), " must be of class kott.design")
    dots <- substitute(list(...))[-1]
    newnames <- names(dots)
    if (length(dots) < 1) 
        return(deskott)
    if (is.null(newnames)) {
        # Solo espressioni prive di tag in (...): esco.
        warning("Untagged input expressions have been dropped")
        return(deskott)
    }
    if (any(newnames %in% names(deskott))) 
        stop("Cannot modify pre-existing ", 
            substitute(deskott)," variables")
    if (any(newnames == "")) {
        # Anche espressioni prive di tag in (...): le rimuovo.
        warning("Untagged input expressions have been dropped")
    }
    full.newnames <- newnames[newnames != ""]
    full.dots <- dots[newnames != ""]
    attr.tocopy <- attributes(deskott)
    for (j in seq(along = full.dots)) {
        deskott[, full.newnames[j]] <- eval(full.dots[[j]], deskott, 
            parent.frame())
    }
    attributes(deskott) <- attr.tocopy
    names(deskott) <- c(attr.tocopy[["names"]], full.newnames)
    deskott
}