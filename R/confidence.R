`confidence` <-
function (estim, se, df, alpha)
################################################
#  Calcolo dell'intervallo di confidenza       #
#  derivato dalla distribuzione t di Student   #
#  (vettorizzato su 'estim' e 'se').           #
################################################
{
    conf.int <- function(estim, se, df, alpha) {
    ################################################
    #  Calcolo dell'intervallo di confidenza       #
    #  derivato dalla distribuzione t di Student.  #
    ################################################
        upper.prob <- (1 - 0.5 * (1 - alpha))
        q <- qt(upper.prob, df)
        l <- (estim - q * se)
        u <- (estim + q * se)
        out <- c(l, u)
        names(out) <- c("l.conf", "u.conf")
        out
    }
    if (length(estim) == 1 && length(se) == 1) 
        return(conf.int(estim, se, df, alpha))
    out.mat <- t(mapply(FUN = conf.int, estim = estim, se = se, 
        MoreArgs = list(df = df, alpha = alpha)))
    dimnames(out.mat) <- list(NULL, c("l.conf", "u.conf"))
    out.mat
}