`global` <-
function(deskott)
###########################################################
#  Consente di referenziare in modo GLOBALE il dataframe  #    
#  associato a 'deskott' nel body di uno stimatore        #
#  'user.estimator' definito dall'utente.                 #
#  Qualora venga richiesta la stima di 'user.estimator'   #
#  in sottopopolazioni, il dataframe ritornato da global  #
#  NON VERRA' SEZIONATO secondo i livelli delle variabili #
#  di 'by'.                                               #
###########################################################
{
if (!inherits(deskott, "kott.design"))
    stop("Object ", substitute(deskott), " must be of class kott.design")
if (!any( k <- sapply(sys.calls(),function(call) call[[1]]=="kottby.user") ))
   {
    stop("There is not any kottby.user call in the call stack")
   }
pos<-match(TRUE,k)   
last.env<-sys.frames()[[pos]]
deskott.name<-sys.calls()[[pos]][[2]]   
eval(deskott.name,envir=last.env)
}