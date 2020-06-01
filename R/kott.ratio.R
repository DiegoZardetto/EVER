`kott.ratio` <-
function (deskott, num, den, by = NULL, vartype = c("se", "cv", "cvpct", "var"), 
    conf.int = FALSE, conf.lev = 0.95)
#######################################################################################
#  Calcola (su oggetti di classe kott.design) la stima del rapporto fra totali        #
#  (o medie) di variabili numeriche ed i corrispondenti errori standard               #
#  ed intervalli di confidenza, nelle sottopopolazioni definite dai livelli delle     #
#  variabili di 'by'.                                                                 #
#  NOTA: La formula da passare per 'by' deve essere del tipo by = ~var1 : ... : varn  #
#        (ogni operatore nella formula verra' comunque interpretato come ":").        #
#  NOTA: Gli intervalli di confidenza sono calcolati usando la distribuzione t di     #
#        Student con nrg-1 gradi di liberta'.                                         #
#  NOTA: La funzione chiama kottby.user() passando al parametro 'user.estimator' la   #
#        funzione ratio().                                                            #
#######################################################################################
{
    if (!inherits(num, "formula")) 
        stop("Numerator variables must be supplied as a formula")
    num.char <- names(model.frame(num, deskott[1, ]))
    na.fail(deskott, num.char)
    num.typetest <- sapply(num.char, function(y) is.numeric(deskott[, 
        y]))
    if (!all(num.typetest)) 
        stop("Numerator variables must be numeric")
    if (!inherits(den, "formula")) 
        stop("Denominator variables must be supplied as a formula")
    den.char <- names(model.frame(den, deskott[1, ]))
    na.fail(deskott, den.char)
    den.typetest <- sapply(den.char, function(y) is.numeric(deskott[, 
        y]))
    if (!all(den.typetest)) 
        stop("Denominator variables must be numeric")
    #if (length(num.char)!=length(den.char))
    #   warning("Different number of variables for numerator and denominator: recycling rule")
    if (missing(vartype)) 
        vartype <- "se"
    vartype <- match.arg(vartype, several.ok = TRUE)
    vartype <- unique(vartype)
    ratio <- function(d, w, num, den) {
    #############################################
    #  Stimatore rapporto fra totali (o medie)  #
    #  di variabili quantitative.               #
    #############################################
        r <- function(num, den) {
            numest <- sum(d[, w] * d[, num])
            denest <- sum(d[, w] * d[, den])
            if (isTRUE(all.equal(denest,0))) {
                warning("Denominator total estimate vanishes in (sub)population: ",
                        "return NaN for ratio estimate and/or SE")
                return(NaN)
            }
            numest/denest
            }
        out <- mapply(r, num, den)
        names(out) <- paste(num, den, sep = "/")
        out
    }
    out <- kottby.user(deskott = deskott, by = by, user.estimator = ratio, 
        vartype = vartype, conf.int = conf.int, conf.lev = conf.lev, df = attr(deskott, "nrg") - 1, 
        num = num.char, den = den.char)
    # Formatta la lista out: un dataframe per ogni livello delle variabili di by (se specificate)
    if (is.null(by)) {
        res <- as.data.frame(out, row.names = names(out[[1]]))
        names(res) <- names(out)
        res
    }
    else {
        res <- lapply(out, as.data.frame, row.names = names(out[[1]][[1]]))
        for(i in seq_along(res)) {colnames(res[[i]]) <- names(out[[1]])}
        res
    }
}