`cal.raking` <-
list(
        Fm1=function (u, bounds){
            pmin(pmax(exp(u), bounds[1]), bounds[2]) - 1
        },
        dF= function (u, bounds){
            ifelse(u < bounds[2] - 1 & u > bounds[1] - 1, exp(u), 0)
        },
        name="raking"
    )
class(cal.raking) <- "calfun"