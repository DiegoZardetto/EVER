`kott.quantile` <- 
function (deskott, y, probs = c(0.25, 0.5, 0.75), by = NULL, 
    vartype = c("se", "cv", "cvpct", "var"), conf.int = FALSE, conf.lev = 0.95) 
#######################################################################################
#  Calcola (su oggetti di classe kott.design) la stima dei quantili di una variabile  #
#  numerica ed i corrispondenti errori standard ed intervalli di confidenza, nelle    #
#  sottopopolazioni definite dai livelli delle variabili di 'by'.                     #
#  NOTA: La formula da passare per 'by' deve essere del tipo by = ~var1 : ... : varn  #
#        (ogni operatore nella formula verra' comunque interpretato come ":")         #
#  NOTA: Se esiste un valore osservato di y, y*, tale che FCUM(y*)=probs la funzione  #
#        ritorna y*. Altrimenti la funzione interpola linearmente fra i 2 valori      #
#        osservati con valori di FCUM piu' prossimi a probs, y*- e y*+.               #
#  NOTA: Gli intervalli di confidenza sono calcolati usando la distribuzione t di     #
#        Student con nrg-1 gradi di liberta'.                                         #
#  NOTA: La funzione chiama kottby.user() passando al parametro 'user.estimator' la   #
#        funzione w.quantiles().                                                      #
#######################################################################################
{
    if (!inherits(y, "formula")) 
        stop("Variable of interest must be supplied as a formula")
    y.char <- names(model.frame(y, deskott[1, ]))
    if (length(y.char) > 1) 
        stop("Can specify only one variable of interest")
    na.fail(deskott, y.char)
    typetest <- sapply(y.char, function(y) is.numeric(deskott[, 
        y]))
    if (!typetest) 
        stop("Variables of interest must be numeric")
    if (any(probs < 0 | probs > 1)) 
        stop("'probs' values must be between 0 and 1")
    if (missing(vartype)) 
        vartype <- "se"
    vartype <- match.arg(vartype, several.ok = TRUE)
    vartype <- unique(vartype)
    w.quantiles <- function(d, w, y, probs) {
    #########################################################
    #  Calcolo della stima campionaria dei quantili di una  #
    #  variabile numerica.                                  #
    #  NOTA: La funzione interpola linearmente fra i valori #
    #        osservati piu' prossimi a FCUM^-1(probs) per   #
    #        ottenere i quantili corrispondenti.            #
    #########################################################
        FCUM <- function(d, w, y) {
        #############################################
        #  Distribuzione di Frequenza Cumulativa.   #
        #  Stima, per ogni valore osservato della   #
        #  variabile y nel campione (yobs), la      #
        #  frazione di popolazione che possiede     #
        #  valori di y minori o uguali a yobs.      #
        #  NOTA: Nelle repliche DAGJK possono       #
        #        esistere unita' con peso zero:     #
        #        il loro contributo deve essere     #
        #        eliminato, in analogia a quanto    #
        #        capiterebbe per un totale.         #
        #############################################
            # rimuovo pesi zero e osservazioni corrispondenti
            yobs <- d[d[, w] > 0, y]
            wobs <- d[d[, w] > 0, w]
            # se il vettore risultante ha zero componenti ritorno NaN
            if (length(yobs) == 0) {
                warning("No positive weights in (sub)population: return NaN for SE")
                return(list(yobs = NaN, F = NaN))
            }
            # altrimenti stimo la distribuzione di frequenza cumulativa
            yobs.o <- sort(yobs)
            wobs.o <- wobs[order(yobs)]
            F <- cumsum(wobs.o)/sum(wobs.o)
            list(yobs = yobs.o, F = F)
        }
        FF <- FCUM(d, w, y)
        yobs <- FF$yobs
        F <- FF$F
        w.qu <- function(p) {
        ########################################################
        #  Inversione (con eventuale interpolazione) di FCUM   #
        #  nel singolo punto p: FCUM^-1(p)                     #
        #  ATTENZIONE: Quando chiamero' w.qu passero' i valori #
        #              di yobs e F con le scoping rules di     #
        #              variabili unbound.                      #
        ########################################################
            # Caso di Warning in FCUM
            if ((length(yobs) == 1) && (length(F) == 1)) {
                 if (is.nan(yobs) && is.nan(F)) {
                     return(NaN)
                    }
                }
            diff <- zapsmall(F - p)
            # 1) Se FCUM e' sempre maggiore o uguale (minore o uguale)
            #    a p stimo qu con il minimo (massimo) dei valori di y
            if (all(diff >= 0)) 
                return(min(yobs))
            if (all(diff <= 0)) 
                return(max(yobs))
            # 2) Altrimenti, se esistono valori y*1<=...<=y*k
            #    tali che FCUM(y*i)=p allora qu=(y*1+y*k)/2
            #    (k puo' essere maggiore di 1 solo a causa di zapsmall)
            if (any(eq <- (diff == 0))) 
                return(mean(range(yobs[eq])))
            # 3)   Altrimenti esistono valori di FCUM
            #      sia maggiori che minori di p:
            # 3.1) individuo il massimo valore y*- tale che FCUM(y*-)
            #      approssima meglio p da sotto
            eps.low <- max(diff[diff < 0])
            eps.low.yobs <- max(yobs[diff == eps.low])
            # 3.2) individuo il minimo valore y*+ tale che FCUM(y*+)
            #      approssima meglio p da sopra
            eps.up <- min(diff[diff > 0])
            eps.up.yobs <- min(yobs[diff == eps.up])
            # 3.3) Se (y*-)=(y*+)=y* (possibile solo per approssimazione numerica)
            #      ritorno y*, altrimenti interpolo linearmente
            if (isTRUE(all.equal(eps.low.yobs, eps.up.yobs))) {
                qu <- eps.low.yobs
            }
            else {
                qu <- eps.low.yobs - eps.low * (eps.up.yobs - 
                  eps.low.yobs)/(eps.up - eps.low)
            }
            qu
        }
        out <- sapply(probs, function(p) w.qu(p))
        names(out) <- paste(100 * probs, "%", sep = "")
        out
    }
    out <- kottby.user(deskott = deskott, by = by, user.estimator = w.quantiles, 
        vartype = vartype, conf.int = conf.int, conf.lev = conf.lev, df = attr(deskott, "nrg") - 1, 
        y = y.char, probs = probs)
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