`z.calib` <-
function (design, formula, population, aggregate.stage = NULL, 
          bounds = c(-Inf, Inf), calfun = c("linear", "raking", "logit"),
          maxit = 50, epsilon = 1e-07, force = FALSE)
########################################################################
#  Versione modificata della funzione calibrate del package survey.    #
#  NOTA: La funzione ritorna DIRETTAMENTE i pesi calibrati e non un    #
#        oggetto design (questo fa risparmiare memoria).               #
#  NOTA: Se aggregate.stage non e' NULL, per diminuire l'uso di        #
#        memoria e la complessita' computazionale, la funzione:        #
#        1) media la model matrix (mm) e somma i pesi diretti (ww)     #
#           nei cluster di stadio aggregate.stage, RIDUCENDO LA        #
#           DIMENSIONE DELL'INPUT di z.grake;                          #
#        2) calibra i DATI RIDOTTI;                                    #
#        3) espande gli g-weights calcolati in 2) sulle unita' finali. #
#  NOTA: Per motivi di efficienza sono state eliminate alcune          #
#        funzionalita' originali NON NECESSARIE per il package EVER    #
#        (ad esempio la QR decomposition).                             #
#  NOTA: Per consentire la diagnosi dei processi di calibrazione       #
#        complessi e strutturati in molti sub-tasks (tipici di         #
#        kottcalibrate) e' stato introdotto il NUOVO attributo         #
#        'ret.code' per l'oggetto (i pesi calibrati) ritornato dalla   #
#        funzione:                                                     #
#        - ret.code=0 se la convergenza E' ottenuta                    #
#        - ret.code=1 se la convergenza NON E' ottenuta MA force=TRUE  #
#  NOTA: z.calibrate chiama z.grake.                                   #
########################################################################
{
    if (is.character(calfun)) 
        calfun <- match.arg(calfun)
    if (is.character(calfun) && calfun == "linear" && (bounds == 
        c(-Inf, Inf))) {
        cal.weights <- z.regcalib(design, formula, population, aggregate.stage = aggregate.stage)
        attr(cal.weights,"ret.code") <- 0
        return(cal.weights)
    }
    if (is.character(calfun)) 
        calfun <- switch(calfun, linear = cal.linear, raking = cal.raking, 
            logit = cal.logit)
    else if (!inherits(calfun, "calfun")) 
        stop("'calfun' must be a string or of class 'calfun'.")
    mm <- model.matrix(formula, model.frame(formula, data = design$variables))
    ww <- design$dir.weights
    sample.total <- colSums(mm * ww)
    if (length(sample.total) != length(population)) 
        stop("Population and sample totals are not the same length.")
    if (any(sample.total == 0)) {
        zz <- (population == 0) & (apply(mm, 2, function(x) all(x == 
            0)))
        mm <- mm[, !zz, drop = FALSE]
        population <- population[!zz]
        sample.total <- sample.total[!zz]
    }
    if (!is.null(aggregate.stage)) {
        aggindex <- design$cluster[[aggregate.stage]]
        mm <- apply(mm, 2, function(col) tapply(col, aggindex, mean))
        ww <- as.numeric(tapply(ww, aggindex, sum))
    }
    g <- z.grake(sample.total, mm, ww, calfun, bounds = bounds, population = population, 
        epsilon = epsilon, maxit = maxit)
    if (!force && !is.null(attr(g, "failed"))) 
        stop("Calibration failed")
    ret.code <- if (is.null(attr(g, "failed"))) 0 else 1
    if (!is.null(aggregate.stage)) {
        g <- g[aggindex]
    }
    cal.weights <- g*design$dir.weights
    attr(cal.weights,"ret.code") <- ret.code  
    cal.weights
}