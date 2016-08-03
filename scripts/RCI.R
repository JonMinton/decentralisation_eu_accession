RCI <- function(povvec, popvec, order)
{
    
    #Storage vectors for X = ordered vector of people "in poverty"; Y = "not in pov"
    X <- 1:length(povvec)
    Y <- X
    
    #povvec[i] is the number of people in poverty (on JSA, Income Support or whatever) in area i
    #X[i] is total number of people in poverty in the i units closest to city centre,
    #divided by overall number of people in poverty.
    for(i in 1:length(X))
    {
        X[i] <- sum(povvec[order[1:i]])/sum(povvec)
        Y[i] <- sum(popvec[order[1:i]] - povvec[order[1:i]])/sum(popvec - povvec)
    }
    
    #Variable to hold result
    res.rce = 0
    for(i in 2:length(X))
    {
        res.rce <- res.rce + X[i-1] * Y[i] - X[i] * Y[i-1]
    }
    
    return(res.rce)
}





RCI.path <- function(povvec, popvec, order, K)
{
    
    #Storage vectors for X = ordered vector of people "in poverty"; Y = "not in pov"
    X <- 1:length(povvec)
    Y <- X
    
    #povvec[i] is the number of people in poverty (on JSA, Income Support or whatever) in area i
    #X[i] is total number of people in poverty in the i units closest to city centre,
    #divided by overall number of people in poverty.
    for(i in 1:K)
    {
        X[i] <- sum(povvec[order[1:i]])/sum(povvec[order[1:K]])
        Y[i] <- sum(popvec[order[1:i]] - povvec[order[1:i]])/sum(popvec[order[1:K]] - povvec[order[1:K]])
    }
    
    #Variable to hold result
    res.rce = 0
    for(i in 2:K)
    {
        res.rce <- res.rce + X[i-1] * Y[i] - X[i] * Y[i-1]
    }
    
    return(res.rce)
}
