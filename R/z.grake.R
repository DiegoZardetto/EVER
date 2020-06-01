`z.grake` <-
function (sample.total, mm, ww, calfun, eta = rep(0, NCOL(mm)), bounds, population, 
    epsilon, maxit)
###################################################################
#  Versione modificata della funzione grake del package survey.   #
#  NOTA: Per motivi di efficienza sono state eliminate alcune     #
#        funzionalita' originali NON NECESSARIE per il package    #
#        EVER (require di MASS e attr(g,"eta")).                  #
###################################################################
{
    if (!inherits(calfun, "calfun")) 
        stop("'calfun' must be of class 'calfun'")
    Fm1 <- calfun$Fm1
    dF <- calfun$dF
    xeta <- drop(mm %*% eta)
    g <- 1 + Fm1(xeta, bounds)
    iter <- 1
    repeat ({
        Tmat <- crossprod(mm * ww * dF(xeta, bounds), mm)
        misfit <- (population - sample.total - colSums(mm * ww * 
            Fm1(xeta, bounds)))
        deta <- MASS::ginv(Tmat, tol = 256 * .Machine$double.eps) %*% 
            misfit
        eta <- eta + deta
        xeta <- drop(mm %*% eta)
        g <- 1 + Fm1(xeta, bounds)
        misfit <- (population - sample.total - colSums(mm * ww * 
            Fm1(xeta, bounds)))
        if (all(abs(misfit)/(1 + abs(population)) < epsilon)) 
            break
        iter <- iter + 1
        if (iter > maxit) {
            achieved <- abs(misfit)/(1 + abs(population))
            worst.achieved <- max(achieved)
            warning("Failed to converge: eps= ", worst.achieved, " in ", 
                iter, " iterations (variable ",names(which.max(achieved)),")")
            attr(g, "failed") <- achieved
            break
        }
    })
    g
}