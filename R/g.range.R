`g.range` <-
function (cal.deskott)
#################################################################################
#  Dato un oggetto di classe kott.cal.design, calcola il range degli g-weights  #
#  sul dataframe originario e su tutti i random samples.                        #
#################################################################################
{
    if (!inherits(cal.deskott, "kott.cal.design")) 
        stop("Object ", substitute(cal.deskott), " must be of class kott.cal.design")
    nrg <- attr(cal.deskott, "nrg")
    w.cal.char <- as.character(attr(cal.deskott, "weights"))[2]
    w.char <- substr(w.cal.char, 0, nchar(w.cal.char) - 4)
    range.r <- function(cal.deskott, w, w.cal, r) {
    ################################################
    #  Dato un oggetto di classe kott.cal.design,  #
    #  calcola il range degli g-weights sul suo    #
    #  r-esimo random sample.                      #
    ################################################
        wr <- cal.deskott[, paste(w, r, sep = "")]
        w.calr <- cal.deskott[, paste(w.cal, r, sep = "")]
        g <- w.calr/wr
        g.nan <- g[is.nan(g)]
        if (length(g.nan) == 0) {
            ranger <- range(g)
        }
        else {
            ranger <- range(g[!is.nan(g)])
            if (!all(wr[is.nan(g)] == w.calr[is.nan(g)] & wr[is.nan(g)] == 
                0)) {
                warning("NaN g-weights not arising from 0/0 in random group ", 
                  r)
            }
            else {
                ranger <- range(ranger, 1)
            }
        }
        ranger
    }
    # original data
    out.matrix <- rbind(range.r(cal.deskott, w.char, w.cal.char,
        NULL))
    # replicated data
    lapply(1:nrg, function(r) out.matrix <<- rbind(out.matrix,
        range.r(cal.deskott, w.char, w.cal.char, r)))
    dimnames(out.matrix) <- list(NULL, c("g.min", "g.max"))
    names.col <- data.frame(sample = c("original", paste("replicate.",
        1:nrg, sep = "")))
    cbind(names.col, out.matrix)
}