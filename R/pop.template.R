`pop.template` <-
function (data, calmodel, partition = FALSE)
###################################################################
#  Costruisce il template del dataframe dei totali noti per un    #
#  determinato modello di calibrazione.                           #
#  La funzione kottcalibrate() esige dataframe dei totali noti    #
#  conformi a tale template.                                      #
#  NOTA: 'calmodel' specifica la formula che definisce il         #
#         modello di calibrazione (in particolare le variabili    #
#         ausiliarie).                                            #
#  NOTA: 'partition' specifica la formula che definisce i factor  #
#        la cui interazione individua le sottopopolazioni per le  #
#        quali sono noti i totali di TUTTE le variabili           #
#        ausiliarie (i "domini di calibrazione").                 #
#        La formula deve essere del tipo:                         #
#          partition=~var1:...:varn                               #
#        (ogni operatore nella formula verra' comunque            #
#        interpretato come ":")                                   #
#        Il valore FALSE per il parametro 'partition' (opzione di #
#        default) indica un modello di calibrazione globale.      #
#  NOTA: Se l'argomento 'partition' viene fornito in modo         #
#        esplicito e' ERRATO riportarne le variabili nella        #
#        formula 'calmodel'. Il parametro 'calmodel' deve         #
#        indicare la formula che definisce il modello di          #
#        calibrazione comune a tutti i domini di calibrazione.    #
###################################################################
{
    if (!inherits(data, "data.frame")) 
        stop("Survey data must be supplied as a data frame")
    if (!inherits(calmodel, "formula")) 
        stop("Parameter 'calmodel' must be supplied as a formula")
    calmodel.vars <- names(model.frame(calmodel, data[1, ]))
    na.fail(data, calmodel.vars)
    calmodel.names <- colnames(model.matrix(calmodel, model.frame(calmodel, 
        data[1, ])))
    if (!identical(partition, FALSE)) {
        if (!inherits(partition, "formula")) 
            stop("Parameter 'partition' must be supplied as a formula")
        partition.vars <- names(model.frame(partition, data[1, 
            ]))
        na.fail(data, partition.vars)
        typetest <- sapply(partition.vars, function(v) is.factor(data[, 
            v]))
        if (!all(typetest)) 
            stop("Partition variables must be factors")
        if (any(partition.vars %in% calmodel.vars)) 
            stop("Calibration model formula cannot reference partition variables")
        u <- unique(data[, partition.vars, drop = FALSE])
        us <- u[do.call(order, u), , drop = FALSE]
        calmodel.tot <- as.data.frame(matrix(data = as.numeric(NA), 
            nrow = nrow(us), ncol = length(calmodel.names), dimnames = list(NULL, 
                calmodel.names)))
        template <- cbind(us, calmodel.tot)
    }
    else {
        template <- as.data.frame(matrix(data = as.numeric(NA), 
            nrow = 1, ncol = length(calmodel.names), dimnames = list(NULL, 
                calmodel.names)))
    }
    row.names(template) <- NULL
    attr(template, "calmodel") <- calmodel
    attr(template, "partition") <- partition
    attr(template, "data") <- substitute(data)
    class(template) <- c("pop.totals", class(template))
    template
}