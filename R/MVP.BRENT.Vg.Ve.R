MVP.BRENT.Vg.Ve <-
function(y, X, eigenK)
{
    p=0
    Sigma <- eigenK$values
    w <- which(Sigma < 1e-6)
    Sigma[w] <- 1e-6
    U <- eigenK$vectors
    min_h2 = 0
    max_h2 = 1
    tol = .Machine$double.eps^0.25
    reml <- brent(y, X, p, Sigma, U, min_h2, max_h2, tol, FALSE, TRUE)
    vg <- reml[[2]]
    ve <- reml[[1]]
    delta <- ve / vg
    return(list(vg=vg, ve=ve, delta=delta))
}
