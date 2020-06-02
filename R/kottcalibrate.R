`kottcalibrate` <-
function (deskott, df.population, 
    calmodel = if (inherits(df.population, "pop.totals"))
                   attr(df.population, "calmodel"), 
    partition = if (inherits(df.population, "pop.totals")) 
                    attr(df.population, "partition") else FALSE, 
    calfun = c("linear", "raking", "logit"), bounds = c(-Inf, Inf), 
    aggregate.stage = NULL, maxit = 50, epsilon = 1e-07, force.rep = FALSE)
#######################################################################
#  Definisce un oggetto di classe kott.cal.design                     #
#  NOTA: L'algoritmo e' iterativo (calibrazioni indipendenti sui      #
#        singoli domini) o non iterativo (calibrazione globale) a     #
#        seconda della struttura del dataframe dei totali noti        #
#        'df.population'.                                             #
#  NOTA: Il programma verifica se 'df.population' sia conforme allo   #
#        standard richiesto da kottcalibrate (la struttura standard   #
#        e' quella generata dalla funzione pop.template); se cosi'    #
#        e', allora la struttura di 'df.population' identifica in     #
#        modo univoco gli eventuali domini di calibrazione e la       #
#        formula che definisce il modello di calibrazione.            #
#        Il programma e', quindi, in grado di generare                #
#        automaticamente, a partire da 'df.population', i parametri   #
#        'formula e 'population' della z.calibrate.                   #
#  NOTA: L'argomento 'force.rep' consente di richiedere, SOLO PER I   #
#        PESI REPLICATI, di accettare l'approssimazione finale dei    #
#        pesi "calibrati" in caso di mancata convergenza dello        #
#        algoritmo di calibrazione.                                   #
#######################################################################
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
                substitute(df.population), " (2). Value (2) will be used", immediate. = TRUE)
        if (!all.equal(partition, attr(df.population, "partition"))) 
            warning("'partition' formula (1) does not agree with the 'partition' attribute of ", 
                substitute(df.population), " (2). Value (2) will be used", immediate. = TRUE)
        if (!identical(attr(deskott, "data"), attr(df.population, 
            "data"))) 
            warning("Data frames used to build objects ", 
                substitute(deskott), " and ", substitute(df.population), 
                " have different names", immediate. = TRUE)
    }
    else {
        df.population <- population.check(df.population, deskott, 
            calmodel, partition)
    }
    if ( !is.numeric(bounds) || (bounds[1] > 1) || (1 > bounds[2]) )
        stop("Bounds must be numeric and must satisfy bounds[1] <= 1 <= bounds[2]")
    if (!is.logical(force.rep)) 
        stop("Parameter 'force.rep' must be logical")
    if (is.character(calfun)) 
        calfun <- match.arg(calfun)
    if (calfun=="logit") {
        if (any(is.infinite(bounds)))
            stop("Bounds must be finite for logit calibration")
        if (any(zapsmall(bounds)==1))
            stop("Logit distance is not defined for bounds[1]=1 or bounds[2]=1")
    }
    calmodel <- attr(df.population, "calmodel")
    calmodel.vars <- names(model.frame(attr(df.population, "calmodel"), 
        deskott[1, ]))
    na.fail(deskott, calmodel.vars)
    partition <- attr(df.population, "partition")
    ids <- attr(deskott, "ids")
    ids.char <- names(model.frame(ids, deskott[1, ]))
    stages <- length(ids.char)
    weights <- attr(deskott, "weights")
    weights.char <- names(model.frame(weights, deskott[1, ]))
    nrg <- attr(deskott, "nrg")
    mk.kottcal.status <- function(nrg,df.population,partition,partition.names=NULL){
    ############################################################
    #  Diagnostica sul processo di calibrazione.               #
    #  Crea nel .GlobalEnv la lista a due componenti           #
    #  'kottcal.status':                                       #
    #  - 'call':                                               #
    #    identifica la chiamata di kottcalibrate che ha        #
    #    generato la lista 'kottcal.status';                   #
    #  - 'return.code':                                        #
    #    matrice che contiene i codici di ritorno dei singoli  #
    #    sub-task del processo di calibrazione complessivo:    #
    #     -1 -> task non ancora affrontato;                    #
    #     0  -> convergenza ottenuta;                          #
    #     1  -> convergenza NON ottenuta ma force=TRUE.        #
    ############################################################
        tasks <- nrg+1
        parts <- nrow(df.population)
        rownam <- c("original", paste("replicate.",1:nrg, sep = ""))
        colnam <- if (identical(partition,FALSE)) "global" else partition.names
        ret.cod <- matrix(-1, nrow = tasks, ncol = parts, dimnames = list(rownam, colnam))
        assign2GE("kottcal.status", list(call = sys.call(-1), return.code = ret.cod))
    }
    upd.kottcal.status <- function(n.sub.task,code){
    ###############################################
    #  Aggiorna la lista 'kottcal.status' nel     #
    #  .GlobalEnv con il codice 'code' ritornato  #
    #  dal sub-task di ordine 'n.subtask'.        #
    ###############################################
        last.task <- col(t(kottcal.status[["return.code"]]))[n.sub.task]
        last.part <- row(t(kottcal.status[["return.code"]]))[n.sub.task]
        kottcal.status[["return.code"]][last.task, last.part] <<- code
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
    cal.sub.task <- 0
    data.out <- deskott
    if (!is.null(aggregate.stage)) {
        if (aggregate.stage < 1 || aggregate.stage > stages) 
            stop("aggregate.stage must be between 1 and the number of sampling stages (", 
            stages, ")")
    }
    if (identical(partition,FALSE)) {
    #####################################################
    #  Calibrazione globale (senza domini)              #
    #####################################################
        mk.kottcal.status(nrg,df.population,partition)
        # calibrazione standard: genera i pesi finali w.cal
        gc.here(need.gc)
        des <- z.design(data = deskott, weights = weights, ids = ids, 
            variables = c(ids.char, weights.char, calmodel.vars))
        gc.here(need.gc)
        w.cal <- z.calib(design = des, population = as.numeric(df.population), 
            formula = calmodel, calfun = calfun, bounds = bounds, 
            aggregate.stage = aggregate.stage, maxit = maxit, epsilon = epsilon, 
            force = FALSE)
        gc.here(need.gc)
        cal.sub.task <- (cal.sub.task+1)
        upd.kottcal.status(cal.sub.task,attr(w.cal,"ret.code"))
        data.out[[paste(weights.char, ".cal",sep = "")]] <- as.numeric(w.cal)
        gc.here(need.gc)
        for (r in 1:nrg) {
        # calibrazione sui random samples r=1,...,nrg: genera i pesi finali w.calr
            weights.char.r <- paste(weights.char, r, sep = "")
            weights.r <- as.formula(paste("~", weights.char.r, sep = "")) 
            des <- z.design(data = deskott, weights = weights.r, ids = ids, 
                variables = c(ids.char, weights.char.r, calmodel.vars))
            gc.here(need.gc)
            w.calr <- z.calib(design = des, population = as.numeric(df.population), 
                formula = calmodel, calfun = calfun, bounds = bounds, 
                aggregate.stage = aggregate.stage, maxit = maxit, epsilon = epsilon, 
                force = force.rep)
            gc.here(need.gc)
            cal.sub.task <- (cal.sub.task+1)
            upd.kottcal.status(cal.sub.task,attr(w.calr,"ret.code"))
            data.out[[paste(weights.char, ".cal", r, sep = "")]] <- as.numeric(w.calr)
            gc.here(need.gc)
        }
    }
    else {
    #####################################################
    #  Calibrazione sui singoli domini di calibrazione  #
    #####################################################
        partition.vars <- names(model.frame(partition,deskott[1, ]))
        na.fail(deskott, partition.vars)
        partition.names <- apply(df.population[,partition.vars,drop=FALSE],1,paste,collapse=".")
        mk.kottcal.status(nrg,df.population,partition,partition.names)
        #  'interact': factor i cui livelli identificano le partizioni
        interact <- interaction(deskott[, rev(partition.vars), drop = FALSE], drop=TRUE)
        #  'groups': lista che contiene gli indici di riga delle osservazioni nelle diverse partizioni
        groups <- split(1:nrow(deskott), interact)
        kottcaliter <- function(wname, force = FALSE) {
        ######################################################################
        #  ATTENZIONE: Quando chiamero' kottcaliter passero' tutto il resto  #
        #              con le scoping rules di variabili unbound.            #
        ######################################################################
            w.cal <- rep(NA,nrow(deskott))
            i.g <- 0
            for (g in groups) {
                i.g <- i.g+1
                weights <- as.formula(paste("~", wname, sep = ""))
                des.g <- z.design(data = deskott[g,], weights = weights, ids = ids, 
                  variables = c(ids.char, wname, calmodel.vars, partition.vars))
                pop.g <- df.population[i.g, which(!(names(df.population) %in% 
                  partition.vars))]
                w.cal.g <- z.calib(design = des.g, population = as.numeric(pop.g), 
                  formula = calmodel, calfun = calfun, bounds = bounds, 
                  aggregate.stage = aggregate.stage, maxit = maxit, epsilon = epsilon,
                  force = force)
                gc.here(need.gc)
                cal.sub.task <<- (cal.sub.task+1)
                upd.kottcal.status(cal.sub.task,attr(w.cal.g,"ret.code"))
                w.cal[g] <- as.numeric(w.cal.g)
            }
            w.cal
        }
        #  calibrazione standard: genera i pesi finali w.cal
        gc.here(need.gc)
        data.out[[paste(weights.char, ".cal", sep = "")]] <- kottcaliter(weights.char, force = FALSE)
        gc.here(need.gc)
        #  calibrazione sui random samples r=1,...,nrg: genera i pesi finali w.calr
        for (r in 1:nrg) {
            data.out[[paste(weights.char, ".cal", r, sep = "")]] <- kottcaliter(paste(weights.char, 
              r, sep = ""), force = force.rep)
            gc.here(need.gc)
        }
    }
    attr(data.out, "data") <- attr(deskott, "data")
    gc.here(need.gc)
    attr(data.out, "ids") <- ids
    gc.here(need.gc)
    attr(data.out, "strata") <- attr(deskott, "strata")
    gc.here(need.gc)
    attr(data.out, "weights") <- as.formula(paste("~", paste(weights.char, 
        ".cal", sep = ""), sep = ""))
    gc.here(need.gc)
    attr(data.out, "self.rep.str") <- attr(deskott, "self.rep.str") 
    gc.here(need.gc)
    attr(data.out, "nrg") <- nrg
    gc.here(need.gc)
    if (inherits(deskott, "kott.cal.design")) {
        attr(data.out, "call") <- c(sys.call(), attr(deskott, 
            "call"))
        class(data.out) <- class(deskott)
        gc.here(need.gc)
    }
    else {
        attr(data.out, "call") <- sys.call()
        class(data.out) <- c("kott.cal.design", class(deskott))
        gc.here(need.gc)
    }
    data.out
}