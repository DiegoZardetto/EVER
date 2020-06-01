`population.check` <-
function (df.population, data, calmodel, partition = FALSE)
####################################################################
#  Verifica se un dataframe di totali noti di popolazione sia,     #
#  o non sia, conforme allo standard richiesto da kottcalibrate    #
#  (la struttura standard e' quella generata dalla funzione        #
#  pop.template).                                                  #
#  Nel primo caso ritorna (senza stamparlo a video) il dataframe   #
#  dei totali noti dopo averlo convertito in un oggetto di classe  #
#  pop.totals, nel secondo caso interrompe l'elaborazione e stampa #
#  un messaggio di errore.                                         #
####################################################################
{
    if (!inherits(df.population, "data.frame")) 
        stop("Known population totals must be supplied as a data frame")
    if (!inherits(data, "data.frame")) 
        stop("Survey data must be supplied as a data frame")
    if (!inherits(calmodel, "formula")) 
        stop("Calibration model must be supplied as a formula")
    if (!identical(partition, FALSE)) {
        if (!inherits(partition, "formula")) 
            stop("Partition variables must be supplied as a formula")
    }
    template <- pop.template(data, calmodel, partition)
    if (!identical(dim(df.population), dim(template))){
        stop.dim <- paste("Dimension of dataframe ", substitute(df.population)," does not agree with 'calmodel' and 'partition' formulas\n(to solve the problem use pop.template)")
        stop(stop.dim)
    }
    if (!identical(names(df.population),names(template))){
        stop.names <- paste("Columns names of data frame ", substitute(df.population)," does not agree with 'calmodel' and 'partition' formulas\n(to solve the problem use pop.template)")
        stop(stop.names)
    }
    if (!identical(partition, FALSE)) {
        test.var.class <- function(df, class) sapply(names(df), 
            function(v) inherits(df[, v], class))
        template.factor <- data.frame(template[, test.var.class(template, 
            "factor"), drop = FALSE])
        df.population.factor <- data.frame(df.population[, test.var.class(df.population, 
            "factor"), drop = FALSE])
        if (!identical(as.matrix(df.population.factor), as.matrix(template.factor))){
            stop.fact <- paste("Columns of data frame ", substitute(df.population)," defining calibration domains\ndoes not agree with 'calmodel' and 'partition' formulas\n(to solve the problem use pop.template)")
            stop(stop.fact)
        }
    }
    template[, ] <- df.population
    print("OK")
    invisible(template)
}