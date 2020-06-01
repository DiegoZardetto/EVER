`z.regcalib` <-
function (design, formula, population, aggregate.stage = NULL)
##########################################################################
#  Versione modificata della funzione regcalibrate del package survey.   #
#  NOTA: La funzione ritorna DIRETTAMENTE i pesi calibrati e non un      #
#        oggetto design (questo fa risparmiare memoria).                 #
#  NOTA: Per motivi di efficienza sono state eliminate alcune            #
#        funzionalita' originali NON NECESSARIE per il package EVER      #
#        (ad esempio la QR decomposition).                               #
##########################################################################
{
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
    whalf <- sqrt(ww)
    g <- rep(1, NROW(mm))
    Tmat <- crossprod(mm * whalf)
    tT <- solve(Tmat, population - sample.total)
    g <- drop(1 + mm %*% tT)
    if (!is.null(aggregate.stage)) {
        g <- g[aggindex]
    }
    cal.weights <- g*design$dir.weights
    cal.weights
}