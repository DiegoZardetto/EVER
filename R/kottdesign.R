`kottdesign` <-
function (data, ids, strata = FALSE, weights, nrg, 
    self.rep.str = FALSE, check.data = FALSE, aux = FALSE)
################################################
#  Definisce un oggetto di classe kott.design  #
################################################    
{
    if (!inherits(data, "data.frame")) 
        stop("Survey data must be supplied as a data frame")
    data.expr <- substitute(data)
    if (!inherits(ids, "formula")) 
        stop("Cluster identifiers must be supplied as a formula")
    ids.charvect <- names(model.frame(ids, data[1, ]))
    na.fail(data, ids.charvect)
    if (!inherits(weights, "formula")) 
        stop("Weigths must be supplied as a formula")
    weights.char <- names(model.frame(weights, data[1, ]))
    if (length(weights.char)>1) 
        stop("Weights formula must reference only one variable")
    na.fail(data, weights.char)
    if (!is.numeric(data[, weights.char])) 
        stop("Variable ", weights.char, " is not numeric")
    if (!is.numeric(nrg)) 
        stop("Parameter 'nrg' must be numeric")
    if (!is.logical(check.data))
        stop("Parameter 'check.data' must be logic")
    if (!is.logical(aux)) 
        stop("Parameter 'aux' must be logic")
    check.nest <- function(data, ids, strata) {
    #################################################
    #  Controlla che le unita' di campionamento di  #
    #  stadio k+1 siano correttamente innestate     #
    #  all'interno di quelle di stadio k.           #
    #  Nota: Il controllo si estende agli strati,   #
    #        che sono trattati come unita' di       #
    #        di stadio k=0.                         #
    #################################################
        nested.unit <- ids
        unit <- c(strata, nested.unit)
        for (i in 1:length(nested.unit)) {
            if (any(tapply(as.numeric(data[, unit[i]]), data[, 
                nested.unit[i]], function(x) length(unique(x)) > 
                1))) {
                if (i == 1) {
                  stop("There is at least one PSU belonging to different strata")
                }
                else {
                  stop("There is at least one stage ", 
                    i, " cluster belonging to different stage ", 
                    i - 1," clusters")
                }
            }
        }
    }
    var.PSU <- function(data, ids.charvect, self.rep.str) {
    ####################################################
    # Crea la variabile var.PSU (identificativo delle  #
    # unita' che forniscono il contributo leading alla #
    # varianza) per disegni di campionamento che       #
    # includono strati auto-rappresentativi.           #
    ####################################################
        var.PSU <- as.character(data[, ids.charvect[1]])
        self <- as.logical(data[, self.rep.str])
        self.ri <- which(self)
        var.PSU[self.ri] <- paste(data[self.ri, ids.charvect[1]],
                                  data[self.ri, ids.charvect[2]], sep=".")
        var.PSU
    }
    if (!identical(strata, FALSE)) {
        if (!inherits(strata, "formula")) 
            stop("Strata must be supplied as a formula")
        strato <- names(model.frame(strata, data[1, ]))
        if (length(strato)>1) 
            stop("Strata formula must reference only one variable")
        na.fail(data, strato)
        if (!is.factor(data[, strato])) 
            stop("Strata variable ", strato, " is not a factor")
        if (check.data) 
            check.nest(data, ids.charvect, strato)
        if (!identical(self.rep.str, FALSE)) {
            if (!inherits(self.rep.str, "formula")) 
                stop("Variable self.rep.str must be supplied as a formula")
            srs.char <- names(model.frame(self.rep.str, data[1, ]))
            if (length(srs.char)>1) 
                stop("self.rep.str formula must reference only one variable")
            na.fail(data, srs.char)
            if (!is.logical(data[, srs.char]) && !is.numeric(data[, srs.char]))
                stop("Variable self.rep.str must be logic or numeric")
            if (length(ids.charvect)==1) {
                self.rep.str <- FALSE
            }
            else {
                data[["var.PSU"]] <- var.PSU(data, ids.charvect, srs.char)
                nvar.PSU <- length(unique(data[, "var.PSU"]))
            }
        }
    }
    else {
        data[["strata.default"]] <- as.factor(1)
        if (check.data && (length(ids.charvect) > 1)) 
            check.nest(data, ids.charvect, "strata.default")
        self.rep.str <- FALSE
    }
    #############################################################
    #  Il parametro need.gc determina se la garbage collection  #
    #  debba, o non debba, essere gestita dal programma.        #
    #  Il valore di soglia per la dimensione di data e' fissata #
    #  ad 1/10 della memoria massima allocabile.                #
    #############################################################
    need.gc <- FALSE
    if (Sys.info()["sysname"] == "Windows"){
        need.gc <- ((8 * prod(dim(data))/(1024^2)) > memory.limit()/10)
    }
    if (need.gc) 
        warning("Survey data frame takes up more than 1/10 of maximum allocable memory", 
                immediate. = TRUE)
    gc.here <- function(doit) {
    #################################################
    #  Se doit=TRUE effettua la garbage collection  #
    #  quando viene invocata.                       #
    #################################################
        if (doit) 
            gc()
    }
    DAG <- function(data, id.PSU, strato, nrg) {
    ###########################################################################
    #  Genera il dataframe dei random groups di PSU in modo conforme alle     #
    #  specifiche del metodo DAGJK stratificato (Kott 1999 e 2001):           #
    #  1. Costruisce il dataframe delle PSU mantenendo per ciascuna solo i    #
    #     valori di strato (rimuove le SSU ed eventuali unita' successive).   #
    #  2. Raggruppa le PSU nei rispettivi strati, ordina i gruppi per strato  #
    #     e le PSU casualmente al loro interno.                               #
    #  3. Crea (per allocazione sistematica) la nuova colonna 'rgi' che       #
    #     indica il random group di appartenenza della singola PSU.           #
    #  4. Aggiunge al dataframe data la colonna 'nh' che, data una generica   #
    #     PSU, indica il numero di PSU presenti nel suo strato.               #
    #  Poi, per r che va da 1 ad nrg:                                         #
    #  5. Crea il vettore colonna colonna 'nhr' che, data una generica PSU,   #
    #     indica il numero di PSU simultaneamente presenti nello strato       #
    #     della PSU e nel random group r.                                     #
    #  6. Crea il vettore colonna colonna 'cwr' che, per una generica PSU,    #
    #     indica il fattore correttivo del peso iniziale della PSU nel        #
    #     random sample r.                                                    #
    ###########################################################################
        nPSU <- length(unique(data[, id.PSU]))
        if (nrg <= 1 || nrg > nPSU) 
            stop("The number of random groups must be greater than 1 and less than or equal to the number of PSU (", 
                nPSU, ")")
        # 1.  costruisco il dataframe PSU-strato
        PSUframe <- unique(data[, c(id.PSU, strato)])
        # 2.1 elimino gli eventuali empty levels di strato in PSUframe
        stratofact <- PSUframe[, strato] <- factor(PSUframe[, strato])
        # 2.2 ordino PSUframe per strato
        PSUframe <- PSUframe[order(stratofact), ]
        stratodist <- table(stratofact)
        lpsustrata <- names(stratodist)[stratodist == 1]
        if (length(lpsustrata) > 0) 
            warning("Lonely PSUs in strata: ", paste(lpsustrata, 
                collapse = ", "), "\n")
        # 2.3 permuto casualmente le PSU di PSUframe nei rispettivi strati (sfrutto PSUframe ordinato per strato)
        stratocum <- c(0, cumsum(stratodist))
        permuta <- rep(NA, nrow(PSUframe))
        sapply(1:length(stratodist), function(i) {
            permuta[(stratocum[i] + 1):stratocum[i + 1]] <<- (stratocum[i] + 
                sample(stratodist[i]))
        })
        PSUframe <- PSUframe[permuta, ]
        # 3.  costruisco la colonna 'rgi' che alloca le PSU nei random groups
        rg <- PSUframe[["rgi"]] <- rep(1:nrg, length.out = nPSU)
        # 4. calcolo il numero 'nh' di PSU per strato (sfrutto PSUframe ordinato per strato)
        nh <- PSUframe[["nh"]] <- rep(stratodist, stratodist)
        for (r in 1:nrg) {
            # 5.1 calcolo il numero 'nhr' di PSU per strato e random group (sfrutto PSUframe ordinato per strato)
            stratofact.r <- PSUframe[rg == r, strato]
            # NOTA: sfrutto il fatto che table su un factor ritorna 0 per i livelli vuoti
            stratodist.r <- table(stratofact.r)
            nhr <- PSUframe[[paste("nh", r, sep = "")]] <- rep(stratodist.r, stratodist)
            # 6.  calcolo il fattore correttivo 'cwr' del peso iniziale delle singole PSU
            Z <- sqrt(nrg/((nrg - 1) * nh * (nh - 1)))
            cwr <- rep(NA, nrow(PSUframe))
            f1 <- nh/(nh - nhr)
            f2 <- 1 + Z
            f3 <- ifelse(nh > 1, 1 - (nh - 1) * Z, 1)
            # NOTA: Se lo strato h contiene una lonely psu (nh=1) allora f3=lim(nh->1)[1-(nh-1)*Z]=1
            # - A) Strati con nh >= nrg: 
            # ---- A1) PSU assenti dal random group r (cioe' presenti nel random sample r):
            cwr[(nh >= nrg) & (rg != r)] <- f1[(nh >= nrg) & (rg != 
                r)]
            # ---- A2) PSU presenti dal random group r (cioe' assenti nel random sample r):
                cwr[(nh >= nrg) & (rg == r)] <- 0
            # - B) Strati con nh < nrg (Modifica di Kott):
            # ---- B1) PSU assenti dal random group r perche' esso non contiene nessuna PSU dello strato:
            cwr[(nh < nrg) & (nhr == 0)] <- 1
            # ---- B2) PSU assenti dal random group r (che pero' contiene alcune PSU dello strato: nhr>0):
            cwr[(nh < nrg) & (rg != r) & (nhr > 0)] <- f2[(nh < 
                nrg) & (rg != r) & (nhr > 0)]
            # ---- B3) PSU presenti nel random group r (che quindi contiene alcune PSU dello strato: nhr>0):
            cwr[(nh < nrg) & (rg == r)] <- f3[(nh < nrg) & (rg == 
                r)]
            PSUframe[[paste("cw", r, sep = "")]] <- cwr
        }
        PSUframe
    }
    mergeDAG <- function(data, data.DAG, id.PSU) {
    ############################################################################
    #  Fonde il dataframe delle PSU prodotto dalla funzione DAG con il         #
    #  dataframe di origine (eliminando i duplicati delle colonne comuni).     #
    #  NOTA: Il risultato e' identico (a meno di un riordinamento delle righe  #
    #        e delle colonne) a quello che si otterrebbe invocando:            #
    #        merge(data,data.DAG) MA IL TEMPO DI ELABORAZIONE E' PIU' PICCOLO. #
    ############################################################################
        espandi <- function(data, quanto) {
            if (length(quanto) != nrow(data)) 
                stop("Expansion vector does not conform to data frame")
            data[rep(1:nrow(data), quanto), ]
        }
        data.DAG <- data.DAG[order(data.DAG[, id.PSU]), ]
        data.DAG <- data.DAG[, -which(names(data.DAG) %in% names(data))]
        data <- data[order(data[, id.PSU]), ]
        gc.here(need.gc)
        fsu.psu.dist <- as.integer(table(data[, id.PSU]))
        # NOTA: Se id.PSU e' un factor con empty levels allora devo rimuovere gli zeri da fsu.psu.dist
        if (is.factor(data[, id.PSU]))
            fsu.psu.dist <- fsu.psu.dist[fsu.psu.dist != 0] 
        cbind(data, espandi(data.DAG, fsu.psu.dist))
    }
    if (identical(self.rep.str, FALSE)) {
        id.PSU <- ids.charvect[1]
    }
    else {
        id.PSU <- "var.PSU"
    }
    nrg <- as.integer(nrg[1])
    if (!identical(strata, FALSE)) {
        PSU.DAG <- DAG(data, id.PSU, strato, nrg)
    }
    else {
        PSU.DAG <- DAG(data, id.PSU, "strata.default", nrg)
    }
    gc.here(need.gc)
    deskott <- mergeDAG(data, PSU.DAG, id.PSU)
    gc.here(need.gc)
    for (r in 1:nrg) {
        deskott[[paste(weights.char, r, sep = "")]] <- deskott[, 
            paste("cw", r, sep = "")] * deskott[, weights.char]
        gc.here(need.gc)
    }
    if (!aux) {
        nhnames <- c("nh", sapply(1:nrg, function(r) paste("nh", 
            r, sep = "")))
        cwnames <- sapply(1:nrg, function(r) paste("cw", r, sep = ""))
        auxnames <- c(nhnames, cwnames)
        deskott <- deskott[, !(names(deskott) %in% auxnames)]
    }
    if (identical(strata, FALSE)) 
        deskott <- deskott[, (names(deskott) != "strata.default")]
    if (!identical(self.rep.str, FALSE)) 
        deskott <- deskott[, (names(deskott) != "var.PSU")]
    gc.here(need.gc)
    attr(deskott, "data") <- data.expr
    gc.here(need.gc)
    attr(deskott, "ids") <- ids
    gc.here(need.gc)
    attr(deskott, "strata") <- strata
    gc.here(need.gc)
    attr(deskott, "weights") <- weights
    gc.here(need.gc)
    attr(deskott, "self.rep.str") <- self.rep.str
    gc.here(need.gc)
    if (!identical(self.rep.str, FALSE)) {
        attr(deskott, "nvar.PSU") <- nvar.PSU
        gc.here(need.gc)
    }
    attr(deskott, "nrg") <- nrg
    gc.here(need.gc)
    attr(deskott, "call") <- sys.call()
    gc.here(need.gc)
    if (!inherits(data, "kott.design")) {
        class(deskott) <- c("kott.design", class(data))
    }
    else class(deskott) <- class(data)
    gc.here(need.gc)
    deskott
}