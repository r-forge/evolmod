"db.default" <-
function(x, ...)
{
tbr<-dualbrothers(123,x,format="interleaved",window.size=400,par_lambda=5,top_lambda=0.693)
tbr
}

