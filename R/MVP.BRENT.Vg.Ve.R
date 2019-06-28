MVP.BRENT.Vg.Ve <-
function(y, X, eigenK)
{
    reml <- brent(Y=y, X=X, eigenK=eigenK, verbose=FALSE)
    vg <- reml[[2]]
    ve <- reml[[1]]
    delta <- ve / vg
    return(list(vg=vg, ve=ve, delta=delta))
}
