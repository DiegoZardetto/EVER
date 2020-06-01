`kott.regcoef` <- 
function (deskott, model, by = NULL, vartype = c("se", "cv", "cvpct", "var"), 
    conf.int = FALSE, conf.lev = 0.95) 
############################################################################################
#  Calcola (su oggetti di classe kott.design) la stima dei coefficienti di regressione di  #
#  un modello lineare ed i corrispondenti errori standard ed intervalli di confidenza,     #
#  nelle sottopopolazioni definite dai livelli delle variabili di 'by'.                    #
#  NOTA: La formula 'model' deve contenere un termine di risposta. Tale termine deve       #
#        essere una singola variabile di tipo numeric.                                     #
#  NOTA: La formula da passare per 'by' deve essere del tipo by = ~var1 : ... : varn       #
#        (ogni operatore nella formula verra' comunque interpretato come ":")              #
#  NOTA: Gli intervalli di confidenza sono calcolati usando la distribuzione t di Student  #
#        con nrg-rank(QR) gradi di liberta'.                                               #
#  NOTA: La funzione chiama kottby.user() passando al parametro 'user.estimator' la        #
#        funzione beta().                                                                  #
############################################################################################
{
    if (!inherits(model, "formula")) 
        stop("Linear model must be supplied as a formula")
    model.vars <- all.vars(model)
    if (!all(model.vars %in% names(deskott)))
        stop("Linear model references variables not available in ",substitute(deskott))
	na.fail(deskott, model.vars)
    if (attr(terms(model),"response")==0)
        stop("Must supply a response variable")
    lhs <- model[[2]]
    if (length(lhs)!=1)
        stop("Must supply only one response variable")
	lhs <- deparse(lhs)
    if (!(is.numeric(deskott[,lhs]))) 
        stop("Response variable must be numeric")	
    lhs <- as.formula(paste("~",lhs,"-1"))
    # Calcolo del numero p di parametri indipendenti nel modello:
    # p e' il rango numerico della QR decomposition della model matrix del rhs di model (sull'intero campione)
    lhs.mat <- model.matrix(lhs,model.frame(lhs,data=deskott,na.action=na.omit,drop.unused.levels=TRUE))
    rhs.mat <- model.matrix(model,model.frame(model,data=deskott,na.action=na.omit,drop.unused.levels=TRUE))
    weights <- attr(deskott, "weights")
    weights.char <- names(model.frame(weights, deskott[1, ]))
    wsqrt <- sqrt(deskott[,weights.char])
    lhs.mat <- lhs.mat*wsqrt
    rhs.mat <- rhs.mat*wsqrt
    q <- qr(rhs.mat)
    model.rank <- q$rank
    # fine p
    nrg <- attr(deskott, "nrg")
    if ( (nrg - model.rank) < 1 )
        stop("Cannot estimate ",model.rank," coefficients with ",nrg," random groups (need ",
             1 + model.rank," at least)")
    if (missing(vartype)) 
        vartype <- "se"
    vartype <- match.arg(vartype, several.ok = TRUE)
    vartype <- unique(vartype)
    beta<-function(d,w,model,lhs){
    ################################################
    #  Calcola i coefficienti beta di regressione  #
    #  per il modello lhs ~ rhs con pesi w.        #
    #  NOTA: Usa la QR-decomposition.              #
    ################################################
        lhs.mat <- model.matrix(lhs,model.frame(lhs,data=d,na.action=na.omit,drop.unused.levels=TRUE))
        rhs.mat <- model.matrix(model,model.frame(model,data=d,na.action=na.omit))
        wsqrt <- sqrt(d[,w])
        lhs.mat <- lhs.mat*wsqrt
        rhs.mat <- rhs.mat*wsqrt
        q <- qr(rhs.mat)
        qr.coef(q,lhs.mat)
    }
    out <- kottby.user(deskott = deskott, by = by, user.estimator = beta, vartype = vartype, 
                conf.int = conf.int, conf.lev = conf.lev, df = attr(deskott, "nrg") - model.rank, 
                model = model, lhs = lhs)
    # Formatta la lista out: un dataframe per ogni livello delle variabili di by (se specificate)
    if (is.null(by)) {	
        res <- as.data.frame(out, row.names = rownames(out[[1]]))
        names(res) <- names(out)
        res
    }
    else {
        res <- lapply(out, as.data.frame, row.names = rownames(out[[1]][[1]]))
        for (i in seq_along(res)) {names(res[[i]]) <- names(out[[1]])}
        res
    }
}