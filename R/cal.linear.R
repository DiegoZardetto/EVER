`cal.linear` <-
list(
        Fm1=function (u, bounds){
            pmin(pmax(u + 1, bounds[1]), bounds[2]) - 1
        },
        dF=function (u, bounds){
            as.numeric(u < bounds[2] - 1 & u > bounds[1] - 1)
        },
        name="linear calibration"
    )
class(cal.linear) <- "calfun"