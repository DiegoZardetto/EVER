`bounds.hint` <-
function (deskott, df.population, 
    calmodel = if (inherits(df.population, "pop.totals"))
                   attr(df.population, "calmodel"), 
    partition = if (inherits(df.population, "pop.totals")) 
                    attr(df.population, "partition") else FALSE)
############################################################################
# Dati (i) un oggetto da calibrare, (ii) i totali noti ed (iii) i metadati #
# che descrivono il modello di calibrazione, la funzione calcola un        #
# intervallo di valori [smallest, greatest] che deve essere                #
# NECESSARIAMENTE interno al range specificato dal parametro 'bounds' di   #
# kottcalibrate (in caso contrario certamente l'algoritmo di calibrazione  #
# non converge).                                                           #
# Sulla base di tale intervallo la funzione costruisce un nuovo intervallo #
# con il medesimo punto medio ed ampiezza doppia [L.sugg, U.sugg]. Tale    #
# intervallo viene suggerito come possibile valore di 'bounds'.            #
# Nota: L'intervallo suggerito NON garantisce la convergenza di            #
#       kottcalibrate ma puo' rivelarsi un buon guess iniziale.            #
# NOTA: La funzione stampa entrambi gli intervalli ma ritorna              #
#       [L.sugg, U.sugg] con attributi che forniscono informazioni piu'    #
#       ricche.                                                            #
# NOTA: In caso di mancata convergenza puo' aiutare: 1) l'ispezione della  #
#       lista di output, 2) l'ispezione di 'kottcal.status'.               #
#       In alternativa e' possibile 1) specificare inizialmente dei bounds #
#       molto ampi (tali da assicurare la convergenza) e 2) restringere    #
#       progressivamente tali bounds con l'aiuto della funzione g.range.   #
############################################################################
{
    if (!inherits(deskott, "kott.design")) 
        stop("Object ", substitute(deskott), " must be of class kott.design")
    if (!inherits(df.population, "data.frame")) 
        stop("Object ", substitute(df.population), " must be of class data frame")
    if (!inherits(calmodel, "formula")) 
        stop("Parameter 'calmodel' must be supplied as a formula")
    if (!identical(partition, FALSE)) {
        if (!inherits(partition, "formula")) 
            stop("Parameter 'partition' must be supplied as a formula")
    }
    if (inherits(df.population, "pop.totals")) {
        if (!all.equal(calmodel, attr(df.population, "calmodel"))) 
            warning("'calmodel' formula (1) does not agree with the 'calmodel' attribute of ", 
                substitute(df.population), " (2). Value (2) will be used",
				immediate. = TRUE)
        if (!all.equal(partition, attr(df.population, "partition"))) 
            warning("'partition' formula (1) does not agree with the 'partition' attribute of ", 
                substitute(df.population), " (2). Value (2) will be used",
				immediate. = TRUE)
        if (!identical(attr(deskott, "data"), attr(df.population, 
            "data"))) 
            warning("Data frames used to build objects ", 
                substitute(deskott),
				" and ", substitute(df.population), " have different names",
				immediate. = TRUE)
    }
    else {
        df.population <- population.check(df.population, deskott, 
            calmodel, partition)
    }
    calmodel <- attr(df.population, "calmodel")
    calmodel.vars <- names(model.frame(attr(df.population, "calmodel"), 
        deskott[1, ]))
    na.fail(deskott, calmodel.vars)
    partition <- attr(df.population, "partition")
    ids <- attr(deskott, "ids")
    ids.char <- names(model.frame(ids, deskott[1, ]))
    weights <- attr(deskott, "weights")
    weights.char <- names(model.frame(weights, deskott[1, ]))
    nrg <- attr(deskott, "nrg")
    mk.bounds <- function(nrg, df.population,partition,partition.names=NULL){
    ############################################################
    # Crea la lista 'bounds' con tre componenti:               #
    #  - 'call':                                               #
    #    identifica la chiamata di bounds.hint che ha generato #
    #    lista 'bounds';                                       #
    #  - 'lower' e 'upper':                                    #
    #    matrici che contengono, per ogni singoli sub-task del #
    #    processo di calibrazione complessivo, il valore       #
    #    minimo e massimo dei rapporti fra totali noti e       #
    #    corrispondenti stime dirette.                         #
    ############################################################
        tasks <- nrg+1
        parts <- nrow(df.population)
        rownam <- c("original", paste("replicate.",1:nrg, sep = ""))
        colnam <- if (identical(partition,FALSE)) "global" else partition.names
        lower <- matrix(-Inf, nrow = tasks, ncol = parts, dimnames = list(rownam, colnam))
        upper <- matrix( Inf, nrow = tasks, ncol = parts, dimnames = list(rownam, colnam))
        bounds <- list(call = sys.call(-1), lower = lower, upper = upper)
        bounds
    }
    upd.bounds <- function(n.sub.task,interval){
    ###############################################
    #  Aggiorna la lista 'bounds' nel             #
    #  .GlobalEnv con i valori low e up ritornati #
    #  dal sub-task di ordine 'n.subtask'.        #
    ###############################################
        last.task <- col(t(bounds[["lower"]]))[n.sub.task]
        last.part <- row(t(bounds[["lower"]]))[n.sub.task]
        bounds[["lower"]][last.task,last.part] <<- interval[1]
        bounds[["upper"]][last.task,last.part] <<- interval[2]
    }
    #############################################################
    #  Il parametro need.gc determina se la garbage collection  #
    #  debba, o non debba, essere gestita dal programma.        #
    #  Il valore di soglia per la dimensione di deskott e'      #
    #  fissata ad 1/10 della memoria massima allocabile.        #
    #############################################################
    need.gc <- FALSE
    if (Sys.info()["sysname"] == "Windows"){
        need.gc <- ((8 * prod(dim(deskott))/(1024^2)) > memory.limit()/10)
    }
    if (need.gc) 
        warning("Input kott.design object takes up more than 1/10 of maximum allocable memory", 
                immediate. = TRUE)
    gc.here <- function(doit) {
    #################################################
    #  Se doit=TRUE effettua la garbage collection  #
    #  quando viene invocata.                       #
    #################################################
        if (doit) 
            gc()
    }
    sub.task <- 0
    if (identical(partition,FALSE)) {
    #####################################################
    #  Intervallo minimo globale (senza domini)         #
    #####################################################
        bounds <- mk.bounds(nrg,df.population,partition)
        gc.here(need.gc)
        # intervallo minimo sul campione originale
        des <- z.design(data = deskott, weights = weights, ids = ids, 
            variables = c(ids.char, weights.char, calmodel.vars))
        gc.here(need.gc)
        interval <- z.ranger(design = des, population = as.numeric(df.population), 
            formula = calmodel)
        gc.here(need.gc)
        sub.task <- (sub.task+1)
        upd.bounds(sub.task,interval)
        gc.here(need.gc)
        for (r in 1:nrg) {
        # intervallo minimo sui random samples r=1,...,nrg
            weights.char.r <- paste(weights.char, r, sep = "")
            weights.r <- as.formula(paste("~", weights.char.r, sep = "")) 
            des <- z.design(data = deskott, weights = weights.r, ids = ids, 
                variables = c(ids.char, weights.char.r, calmodel.vars))
            gc.here(need.gc)
            interval.r <- z.ranger(design = des, population = as.numeric(df.population), 
                formula = calmodel)
            gc.here(need.gc)
            sub.task <- (sub.task+1)
            upd.bounds(sub.task,interval.r)
            gc.here(need.gc)
        }
    }
    else {
    ##########################################################
    #  Intervalli minimi sui singoli domini di calibrazione  #
    ##########################################################
        partition.vars <- names(model.frame(partition,deskott[1, ]))
        na.fail(deskott, partition.vars)
        partition.names <- apply(df.population[,partition.vars,drop=FALSE],1,paste,collapse=".")
        bounds <- mk.bounds(nrg,df.population,partition,partition.names)
        #  'interact': factor i cui livelli identificano le partizioni
        interact <- interaction(deskott[, rev(partition.vars), drop = FALSE], drop=TRUE)
        #  'groups': lista che contiene gli indici di riga delle osservazioni nelle diverse partizioni
        groups <- split(1:nrow(deskott), interact)
        range.iter <- function(wname) {
        #############################################
        # Funzione da ripetere su tutte le repliche #
        #############################################
            i.g <- 0
            for (g in groups) {
                i.g <- i.g+1
                weights <- as.formula(paste("~", wname, sep = ""))
                des.g <- z.design(data = deskott[g,], weights = weights, ids = ids, 
                  variables = c(ids.char, wname, calmodel.vars, partition.vars))
                pop.g <- df.population[i.g, which(!(names(df.population) %in% 
                  partition.vars))]
                interval.g <- z.ranger(design = des.g, population = as.numeric(pop.g), 
                  formula = calmodel)
                gc.here(need.gc)
                sub.task <<- (sub.task+1)
                upd.bounds(sub.task,interval.g)
            }
        }
        # intervallo minimo sul campione originale
        gc.here(need.gc)
        range.iter(weights.char)
        gc.here(need.gc)
        # intervallo minimo sui random samples r=1,...,nrg
        for (r in 1:nrg) {
            range.iter(paste(weights.char, r, sep = ""))
            gc.here(need.gc)
        }
    }

    all <- apply(bounds[["lower"]], 2, min)
    smallest <- min(1, all)
    bounds[["lower"]] <- rbind(bounds[["lower"]], all)
    all <- apply(bounds[["upper"]], 2, max)
    greatest <- max(1, all)
    bounds[["upper"]] <- rbind(bounds[["upper"]], all)
    if (is.finite(smallest) && is.finite(greatest)) {
        mid <- (smallest + greatest)/2
        L <- mid - 2 * (mid - smallest)
        U <- mid + 2 * (greatest - mid)
        L.sugg <- round(L, 3)
        if (isTRUE(all.equal(L.sugg, 1))) 
            L.sugg <- (L.sugg - 0.001)
        U.sugg <- round(U, 3)
        if (isTRUE(all.equal(U.sugg, 1))) 
            U.sugg <- (U.sugg + 0.001)
        suggestion <- c(L.sugg, U.sugg)
    }
    else {
        if (is.infinite(smallest) && is.infinite(greatest)) {
            suggestion <- c(smallest, greatest)
        }
        else {
            if (is.finite(smallest)) {
                mid <- 1
                L <- mid - 2 * (mid - smallest)
                L.sugg <- round(L, 3)
                if (isTRUE(all.equal(L.sugg, 1))) 
                  L.sugg <- (L.sugg - 0.001)
                U.sugg <- greatest
            }
            else {
                mid <- 1
                U <- mid + 2 * (greatest - mid)
                U.sugg <- round(U, 3)
                if (isTRUE(all.equal(U.sugg, 1))) 
                  U.sugg <- (U.sugg + 0.001)
                L.sugg <- smallest
            }
            suggestion <- c(L.sugg, U.sugg)
        }
    }
    attr(suggestion, "star.interval") <- c(smallest, greatest)
    attr(suggestion, "bounds") <- bounds
    class(suggestion) <- c("bounds.hint", class(suggestion))
    cat("\n")
    cat(paste("A starting suggestion: try to calibrate with bounds=c(", 
        L.sugg, ", ", U.sugg, ")\n", sep = ""))
    cat("\n")
    cat("Remark: this is just a hint, not an exact result\n")
    cat(paste("Feasible bounds for calibration problem must cover the interval [", 
        round(smallest, 3), ", ", round(greatest, 3), "]\n", 
        sep = ""))
    cat("\n")
    invisible(suggestion)
}