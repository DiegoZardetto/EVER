`kottby.user` <-
function (deskott, by = NULL, user.estimator, na.replace = NULL, 
    vartype = c("se", "cv", "cvpct", "var"), conf.int = FALSE, conf.lev = 0.95, 
    df = attr(deskott, "nrg") - 1, ...)
#############################################################################
#  Calcola, su oggetti di classe kott.design, la stima, l'errore standard   #
#  e l'intervallo di confidenza di uno stimatore arbitrario (anche privo    #
#  di una rappresentazione analitica), nelle sottopopolazioni definite dai  #
#  livelli delle variabili di 'by'.                                         #
#  NOTA: La formula da passare per 'by' deve essere del tipo:               #
#        by=~var1:...:varn                                                  #
#        (ogni operatore in formula verra' comunque interpretato come ":")  #
#  NOTA: Il valore attuale del parametro formale 'user.estimator' deve      #
#        essere di classe function.                                         #
#        La segnatura della funzione (supponiamo si chiami user.estfun)     #
#        deve essere necessariamente la seguente:                           #
#                                                                           #
#        user.estfun=function(data,wname,...){body}                         #
#                                                                           #
#        Il body di user.estfun() deve contenere tutte le istruzioni che    #
#        consentirebbero di calcolare la stima campionaria dello stimatore  #
#        desiderato a partire dai dati campionari contenuti nel dataframe   #
#        'data' usando i pesi contenuti nella sua colonna di nome 'wname'.  #
#  NOTA: La funzione kottby.user() deve poter invocare la funzione          #
#        user.estfun() sul dataframe dei dati di origine dell'oggetto       #
#        'deskott' (accessibile mediante attr(deskott,"data")).             #
#        Conseguentemente l'utente deve assicurarsi, nel costruire la       #
#        funzione user.estfun, che le istruzioni del body facciano          #
#        riferimento a variabili effettivamente contenute in tale           #
#        dataframe.                                                         #
#  NOTA: L'argomento '...' in kottby.user() consente di passare alla        #
#        funzione user.estfun() parametri ulteriori rispetto a 'data' e     #
#        'wname' (ad esempio valori "esterni" che si desidera non siano     #
#        ricalcolati per ogni random sample in deskott).                    #
#  NOTA: Nel caso in cui la funzione user.estfun() faccia uso di quantita'  #
#        che si riferiscono all'intero oggetto deskott ANCHE QUANDO la      #
#        stima di user.estfun() DEBBA ESSERE CALCOLATA IN SOTTOPOPOLAZIONI  #
#        e' possibile referenziare deskott mediante la funzione 'global()'  #
#        si veda la funzione poverty() come esempio.                        #
#  NOTA: Non e' necessario che il valore di ritorno della funzione passata  #
#        mediante il parametro 'user.estimator' sia un singolo valore       #
#        atomico (puo' essere un vettore, una matrice, un array, ...).      #
#        In ogni caso esso verra' tacitamente convertito in array dalla     #
#        funzione kottby.user().                                            #
#  NOTA: Gli intervalli di confidenza per default sono calcolati usando la  #
#        distribuzione t di Student con nrg-1 gradi di liberta'.            #
#############################################################################
{
    if (!inherits(deskott, "kott.design")) 
        stop("Object ", substitute(deskott), " must be of class kott.design")
    few.obs(deskott)
    if (!inherits(user.estimator, "function")) 
        stop("Object ", substitute(user.estimator), " is not a function")
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
    if (missing(vartype)) 
        vartype <- "se"
    vartype <- match.arg(vartype, several.ok = TRUE)
    vartype <- unique(vartype)
    vartype.pos <- pmatch(vartype, eval(formals(sys.function())$vartype))
    if (any(is.na(vartype.pos))) 
        stop("Unavailable vartype")
    variabilities <- function(se, cv, cvpct, var, which.one){
        var.list <- list(se, cv, cvpct, var)[which.one, drop = FALSE]
		names(var.list) <- c("SE", "CV", "CV%", "Var")[which.one]
        var.list
    }
    if (!is.logical(conf.int)) 
        stop("Parameter 'conf.int' must be logical")
    if (!is.numeric(conf.lev) || conf.lev < 0 || conf.lev > 1) 
        stop("conf.lev must be between 0 and 1")
    if (missing(df)) {
        if (!identical(conf.int, FALSE)) 
            warning("Unspecified degrees of freedom: used default value (df = ", df, ")")
    }
    else {
        if (!identical(conf.int, FALSE) && (!is.numeric(df) || df < 1)) 
            stop("Need at least 1 degree of freedom")
    }
    kottestim <- function(deskott, user.estimator, na.replace = NULL, 
        vartype.pos, conf.int, conf.lev, df, ...) {
    #############################################################################
    #  Calcola, su oggetti di classe kott.design, la stima, l'errore standard   #
    #  e l'intervallo di confidenza di uno stimatore arbitrario (anche privo    #
    #  di una rappresentazione analitica).                                      #
    #  NOTA: Gli intervalli di confidenza sono calcolati usando la              #
    #        distribuzione t di Student con 'df' gradi di liberta'.             #
    #############################################################################
        nrg <- attr(deskott, "nrg")
        w <- attr(deskott, "weights")
        w.char <- names(model.frame(w, deskott[1, ]))
        e <- user.estimator(deskott, w.char, ...)
        e <- as.array(e)
        er <- lapply(1:nrg, function(r) user.estimator(deskott, 
            paste(w.char, r, sep = ""), ...))
        er <- array(unlist(er), dim = c(dim(e), nrg))
        ecycle <- array(rep(e, nrg), dim = c(dim(e), nrg))
        diff2 <- (er - ecycle)^2
        var <- ((nrg - 1)/nrg) * apply(diff2, 1:length(dim(e)), 
            sum)
        var <- as.array(var)
        dimnames(var) <- dimnames(e)
        se <- sqrt(var)
        dimnames(se) <- dimnames(e)
        cv <- se/e
        dimnames(cv) <- dimnames(e)
        cvpct <- 100*cv
        dimnames(cvpct) <- dimnames(e)
        vars <- variabilities(se = se, cv = cv, cvpct = cvpct, var = var, which.one = vartype.pos)
        if (!identical(conf.int, FALSE)) {
            l.conf <- array(mapply(FUN = confidence, e, se, MoreArgs = list(df = df, 
                alpha = conf.lev))[1, ], dim = dim(e))
            dimnames(l.conf) <- dimnames(e)
            l.conf.tag <- paste("l.conf(", round(100*conf.lev,1), "%)", sep="")
            u.conf <- array(mapply(FUN = confidence, e, se, MoreArgs = list(df = df, 
                alpha = conf.lev))[2, ], dim = dim(e))
            dimnames(u.conf) <- dimnames(e)
            u.conf.tag <- paste("u.conf(", round(100*conf.lev,1), "%)", sep="")
        }
        if (!is.null(na.replace)) {
            e[is.na(e)] <- na.replace
            var[is.na(var)] <- na.replace
            se[is.na(se)] <- na.replace
            cv[is.na(cv)] <- na.replace
            cvpct[is.na(cvpct)] <- na.replace
            if (!identical(conf.int, FALSE)) {
                l.conf[is.na(l.conf)] <- na.replace
                u.conf[is.na(u.conf)] <- na.replace
            }
        }
        if (!identical(conf.int, FALSE)) {
            estimate <- list(estimate = e)
			l.conf <- list(l.conf); names(l.conf) <- l.conf.tag
            u.conf <- list(u.conf); names(u.conf) <- u.conf.tag
            c(estimate, vars, l.conf, u.conf)
        }
        else {
            estimate <- list(estimate = e)
            c(estimate, vars)
        }
    }
    if (identical(by, NULL)) {
        kottestim(deskott, user.estimator, na.replace, vartype.pos, 
            conf.int, conf.lev, df, ...)
    }
    else {
        dfby <- deskott[, by.charvect, drop = FALSE]
        lapply(split(deskott, dfby, drop = TRUE), function(des) kottestim(des, 
            user.estimator, na.replace, vartype.pos, conf.int, conf.lev, df, ...))
    }
}