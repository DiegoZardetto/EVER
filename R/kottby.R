`kottby` <-
function (deskott, y, by = NULL, estimator = c("total", "mean"),
    vartype = c("se", "cv", "cvpct", "var"), conf.int = FALSE, conf.lev = 0.95)
#######################################################################################
#  Calcola (su oggetti di classe kott.design) le stime dei totali o delle medie       #
#  di piu' variabili ed i corrispondenti errori standard ed intervalli di confidenza  #
#  nelle sottopopolazioni definite dai livelli delle variabili di 'by'.               #
#  NOTA: A seconda che la variabile di stima y sia di tipo a) numeric, b) factor la   #
#        funzione calcola:                                                            #
#        se estimator="total" -> a) la stima del totale di y                          #
#                                b) la stima delle frequenze assolute di y            #
#        se estimator="mean"  -> a) la stima dela media di y                          #
#                                b) la stima delle frequenze relative di y            #
#  NOTA: La formula da passare per  'y' deve essere del tipo  y = ~var1 + ... + varn  #
#        (ogni operatore nella formula verra' comunque interpretato come "+")         #
#  NOTA: La formula da passare per 'by' deve essere del tipo by = ~var1 : ... : varn  #
#        (ogni operatore nella formula verra' comunque interpretato come ":")         #
#  NOTA: Gli intervalli di confidenza sono calcolati usando la distribuzione t di     #
#        Student con nrg-1 gradi di liberta'.                                         #
#  NOTA: Il valore di ritorno della funzione puo' essere un dataframe o una lista e   #
#        la sua struttura dipende dalla natura dell'input.                            #
#######################################################################################
{
    if (!inherits(deskott, "kott.design")) 
        stop("Object ", substitute(deskott), " must be of class kott.design")
    few.obs(deskott)
    if (!inherits(y, "formula")) 
        stop("Variables of interest must be supplied as a formula")
    y.charvect <- names(model.frame(y, deskott[1, ]))
    na.fail(deskott, y.charvect)
    typetest <- sapply(y.charvect, function(y) is.factor(deskott[, 
        y]) || is.numeric(deskott[, y]))
    if (!all(typetest)) 
        stop("Variables of interest must be numeric or factor")
    if (!identical(by, NULL)) {
        if (!inherits(by, "formula")) 
            stop("'by' variables must be supplied as a formula")
        by.charvect <- names(model.frame(by, deskott[1, ]))
        na.fail(deskott, by.charvect)
        typetest <- sapply(by.charvect, function(y) is.factor(deskott[, 
            y]))
        if (!all(typetest)) 
            stop("'by' variables must be factor")
        few.obs(deskott, by.charvect)
    }
    estimator <- match.arg(estimator)
    if (missing(vartype)) 
        vartype <- "se"
    vartype <- match.arg(vartype, several.ok = TRUE)
    vartype <- unique(vartype)
    vartype.pos <- pmatch(vartype, eval(formals(sys.function())$vartype))
    if (any(is.na(vartype.pos))) 
        stop("Unavailable vartype")
    variabilities <- function(se, cv, cvpct, var, which.one){
        var.mat <- cbind(se, cv, cvpct, var)[, which.one, drop = FALSE]
		colnames(var.mat) <- c("SE", "CV", "CV%", "Var")[which.one]
        var.mat
    }
    if (!is.logical(conf.int)) 
        stop("Parameter 'conf.int' must be logical")
    if (!is.numeric(conf.lev) || conf.lev < 0 || conf.lev > 1) 
        stop("conf.lev must be between 0 and 1")
    kottby1 <- function(deskott, y, by = NULL, estimator, vartype.pos, 
        conf.int, conf.lev) {
    #########################################################
    #  Calcola (su oggetti di classe kott.design) la stima  #
    #  del totale (o della media) di una sola variabile ed  #
    #  il relativo errore standard (ed intervallo di        #
    #  confidenza), nelle sottopopolazioni definite dai     #
    #  livelli delle variabili di 'by'.                     #
    #  NOTA: 'y' e 'by' devono essere vettori character.    #
    #  NOTA: Gli intervalli di confidenza sono calcolati    #
    #        usando la distribuzione t di Student con       #
    #        nrg-1 gradi di liberta'.                       #
    #########################################################
        if (is.null(by)) 
            return(kottestim1(deskott, y, estimator, vartype.pos, conf.int, 
                conf.lev))
        dfby <- deskott[, by]
        yvect <- deskott[, y]
        if (is.numeric(yvect)) {
            out <- sapply(split(deskott, dfby, drop = TRUE), function(des) kottestim1(des, 
                y, estimator, vartype.pos, conf.int, conf.lev))
            return(as.data.frame(out))
        }
        if (is.factor(yvect)) {
            out <- lapply(split(deskott, dfby, drop = TRUE), function(des) kottestim1(des, 
                y, estimator, vartype.pos, conf.int, conf.lev))
            return(out)
        }
    }
    kottestim1 <- function(deskott, y, estimator, vartype.pos, conf.int, conf.lev) {
    #############################################################################
    #  Calcola (su oggetti di classe kott.design) la stima del totale o della   #
    #  media di una sola variabile ed il corrispondente errore standard (ed     #
    #  intervallo di confidenza).                                               #
    #  NOTA: A seconda che la variabile di stima y sia di tipo a) numeric,      #
    #        b) factor la funzione calcola:                                     #
    #        se estimator="total" -> a) la stima del totale di y                #
    #                                b) la stima delle frequenze assolute di y  #
    #        se estimator="mean"  -> a) la stima dela media di y                #
    #                                b) la stima delle frequenze relative di y  #
    #  NOTA: 'y' deve essere di tipo character.                                 #
    #  NOTA: Gli intervalli di confidenza sono calcolati usando la              #
    #        distribuzione t di Student con nrg-1 gradi di liberta'.            #
    #############################################################################
        estim1 <- function(data, y, w, estim) {
            total1 <- function(data, y, w) {
                yvect <- data[, y]
                wvect <- data[, w]
                if (is.numeric(yvect)) 
                  ty <- sum(yvect * wvect)
                if (is.factor(yvect)) {
                  yvect <- factor(yvect)             # rimuove gli (eventuali) empty levels di yvect
                  ty <- tapply(wvect, yvect, sum)
                }
                ty
            }
            mean1 <- function(data, y, w) {
                yvect <- data[, y]
                wvect <- data[, w]
                wsum <- sum(wvect)
                if (is.numeric(yvect)) 
                  ty <- sum(yvect * wvect)
                if (is.factor(yvect)) {
                  yvect <- factor(yvect)             # rimuove gli (eventuali) empty levels di yvect
                  ty <- tapply(wvect, yvect, sum)
                }
                my <- ty/wsum
                my
            }
            switch(estim, total = total1, mean = mean1)
        }
        nrg <- attr(deskott, "nrg")
        w <- attr(deskott, "weights")
        w.char <- names(model.frame(w, deskott[1, ]))
        yvect <- deskott[, y]
        if (is.factor(yvect)) {
            yvect <- factor(yvect)                   # rimuove gli (eventuali) empty levels di yvect
            full.levname <- paste(y, levels(yvect), sep = ".")
        }
        est.fun <- estim1(deskott, y, w.char, estimator)
        e <- est.fun(deskott, y, w.char)
        er <- sapply(1:nrg, function(r) est.fun(deskott, y, paste(w.char, 
            r, sep = "")))
        if (length(e) == 1) {
            # sempre T se yvect e' numeric, se yvect e' factor T solo se ha un unico livello non empty
            var <- ((nrg - 1)/nrg) * sum((er - e)^2)
            se <- sqrt(var)
            cv <- se/e
            cvpct <- 100*cv
            vars <- rbind(variabilities(se = se, cv = cv, cvpct = cvpct, var = var, which.one = vartype.pos))
            if (!identical(conf.int, FALSE)) {
                l.conf <- confidence(estim = e, se = se, df = (nrg - 
                  1), alpha = conf.lev)[1]
                u.conf <- confidence(estim = e, se = se, df = (nrg - 
                  1), alpha = conf.lev)[2]
                out <- cbind(e, vars, l.conf, u.conf)
                l.conf.tag <- paste("l.conf(", round(100*conf.lev,1), "%)", sep="")
                u.conf.tag <- paste("u.conf(", round(100*conf.lev,1), "%)", sep="")
                dimnames(out) <- list(ifelse(!(is.factor(yvect)), 
                  y, full.levname), c(estimator, colnames(vars), l.conf.tag, 
                  u.conf.tag))
            }
            else {
                out <- cbind(e, vars)
                dimnames(out) <- list(ifelse(!(is.factor(yvect)), 
                  y, full.levname), c(estimator, colnames(vars)))
            }
            return(as.data.frame(out))
        }
        else {
            ecol <- cbind(e)
            emat <- matrix(ecol, nrow(er), ncol(er))
            var <- cbind(((nrg - 1)/nrg) * rowSums((er - emat)^2))
            se <- sqrt(var)
            cv <- se/ecol
            cvpct <- 100*cv
            vars <- variabilities(se = se, cv = cv, cvpct = cvpct, var = var, which.one = vartype.pos)
            if (!identical(conf.int, FALSE)) {
                l.conf <- confidence(estim = ecol, se = se, df = (nrg - 
                  1), alpha = conf.lev)[, 1]
                u.conf <- confidence(estim = ecol, se = se, df = (nrg - 
                  1), alpha = conf.lev)[, 2]
                out <- cbind(ecol, vars, l.conf, u.conf)
                l.conf.tag <- paste("l.conf(", round(100*conf.lev,1), "%)", sep="")
                u.conf.tag <- paste("u.conf(", round(100*conf.lev,1), "%)", sep="")
                colnames(out) <- c(estimator, colnames(vars), l.conf.tag, 
                  u.conf.tag)
            }
            else {
                out <- cbind(ecol, vars)
                colnames(out) <- c(estimator, colnames(vars))
            }
            rownames(out) <- full.levname
            return(as.data.frame(out))
        }
    }
    if (identical(by, NULL)) {
        out <- lapply(y.charvect, function(y) kottby1(deskott, 
            y, by, estimator, vartype.pos, conf.int, conf.lev))
    }
    else {
        out <- lapply(y.charvect, function(y) kottby1(deskott, 
            y, by.charvect, estimator, vartype.pos, conf.int, conf.lev))
    }
    names(out) <- y.charvect
    if (length(out) == 1) 
        out <- out[[1]]
    out
}