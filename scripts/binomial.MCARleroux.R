binomial.MCARleroux <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, verbose=TRUE)
{
    #### Check on the verbose option
    if(is.null(verbose)) verbose=TRUE     
    if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)
    
    if(verbose)
    {
        cat("Setting up the model\n")
        a<-proc.time()
    }else{}
    
    
    
    ##############################################
    #### Format the arguments and check for errors
    ##############################################
    #### Overall formula object
    frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
    if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths or the data object has not been specified.", call.=FALSE)
    
    
    #### Design matrix
    ## Create the matrix
    X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
    N.all <- nrow(X)
    p <- ncol(X)
    
    ## Check for linearly related columns
    cor.X <- suppressWarnings(cor(X))
    diag(cor.X) <- 0
    
    if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    
    if(p>1)
    {
        if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
    }else
    {
    }
    
    ## Standardise the matrix
    X.standardised <- X
    X.sd <- apply(X, 2, sd)
    X.mean <- apply(X, 2, mean)
    X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
    
    for(j in 1:p)
    {
        if(length(table(X[ ,j]))>2)
        {
            X.indicator[j] <- 1
            X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
        }else if(length(table(X[ ,j]))==1)
        {
            X.indicator[j] <- 2
        }else
        {
            X.indicator[j] <- 0
        }
    }
    
    
    #### Response variable
    Y <- model.response(frame)
    
    ## Check for errors
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
    int.check <- N.all-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
    
    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
    int.check <- N.all-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)
    if(sum(Y>trials)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)
    Y <- as.numeric(Y)
    failures <- trials - Y
    
    #### Offset variable
    offset <- try(model.offset(frame), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    
    #### Format and check the neighbourhood matrix W
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    K <- nrow(W)
    if(ncol(W)!= K) stop("W has the wrong number of columns.", call.=FALSE)
    if(nrow(W)!= K) stop("W has the wrong number of rows.", call.=FALSE)
    if(floor(N.all/K)!=ceiling(N.all/K)) stop("The number of spatial areas is not a multiple of the number of data points.", call.=FALSE)
    N <- N.all / K
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
    
    #### Specify the initial parameter values
    dat <- cbind(Y, failures)
    beta <- glm(dat~X.standardised-1, offset=offset, family=binomial)$coefficients
    theta.hat <- Y / trials
    theta.hat[theta.hat==0] <- 0.01
    theta.hat[theta.hat==1] <- 0.99
    regression.vec <- X.standardised %*% beta
    res.temp <- log(theta.hat / (1 - theta.hat)) - regression.vec - offset
    phi <- rnorm(n=N.all, mean=0, sd = 2*sd(res.temp))
    phi.mat <- matrix(phi, nrow=K, byrow=TRUE)
    Sigma <- cov(phi.mat)
    Sigma.inv <- solve(Sigma)
    rho <- runif(1)
    
    
    #### Check and specify the priors
    ## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- N+1
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- diag(rep(1,N))
    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)       
    
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)    
    
    #### Format and check the MCMC quantities
    if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
    if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(n.sample <= thin)  stop("thin is greater than n.sample.", call.=FALSE)
    if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
    if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
    if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
    
    
    ## Compute the blocking structure for beta     
    blocksize.beta <- 5 
    if(blocksize.beta >= p)
    {
        n.beta.block <- 1
        beta.beg <- 1
        beta.fin <- p
    }else
    {
        n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
        remainder <- p - n.standard * blocksize.beta
        
        if(remainder==0)
        {
            beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
            beta.fin <- seq(blocksize.beta, p, blocksize.beta)
            n.beta.block <- length(beta.beg)
        }else
        {
            beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
            beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
            n.beta.block <- length(beta.beg)
        }
    }         
    
    
    #### Set up matrices to store samples
    n.keep <- floor((n.sample - burnin)/thin)
    samples.beta <- array(NA, c(n.keep, p))
    samples.phi <- array(NA, c(n.keep, N.all))
    samples.Sigma <- array(NA, c(n.keep, N^2))
    samples.rho <- array(NA, c(n.keep, 1))  
    samples.fitted <- array(NA, c(n.keep, N.all))
    samples.deviance <- array(NA, c(n.keep, 1))
    
    
    #### Specify the Metropolis quantities
    accept.all <- rep(0,6)
    accept <- accept.all
    proposal.sd.phi <- 0.1
    proposal.sd.rho <- 0.05
    proposal.sd.beta <- 0.01
    proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
    chol.proposal.corr.beta <- chol(proposal.corr.beta)     
    Sigma.post.df <- prior.Sigma.df + K
    
    
    
    #### Spatial quantities
    ## Create the triplet object
    W.triplet <- c(NA, NA, NA)
    for(i in 1:K)
    {
        for(j in 1:K)
        {
            if(W[i,j]>0)
            {
                W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
            }else{}
        }
    }
    W.triplet <- W.triplet[-1, ]     
    W.n.triplet <- nrow(W.triplet) 
    W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
    W.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
    
    
    ## Create the start and finish points for W updating
    W.begfin <- array(NA, c(K, 2))     
    temp <- 1
    for(i in 1:K)
    {
        W.begfin[i, ] <- c(temp, (temp + W.neighbours[i]-1))
        temp <- temp + W.neighbours[i]
    }
    
    
    ## Create the determinant     
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    Q.spat <- rho * (diag(apply(W, 1, sum)) - W) + diag(rep(1-rho), K)
    det.Q.spat <-  sum(log((rho * Wstar.val + (1-rho))))    
    
    
    #### Specify quantities that do not change
    offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=TRUE) 
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=TRUE)
    Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=TRUE)
    failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=TRUE)
    
    
    
    ###########################
    #### Run the Bayesian model
    ###########################
    ## Start timer
    if(verbose)
    {
        cat("Collecting", n.sample, "samples\n", sep = " ")
        progressBar <- txtProgressBar(style = 3)
        percentage.points<-round((1:100/100)*n.sample)
    }else
    {
        percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
    
    for(j in 1:n.sample)
    {
        ####################
        ## Sample from beta
        ####################
        proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
        proposal.beta <- beta
        offset.temp <- offset + phi     
        
        for(r in 1:n.beta.block)
        {
            proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
            prob <- binomialbetaupdate(X.standardised, N.all, p, beta, proposal.beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta)
            if(prob > runif(1))
            {
                beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
                accept[1] <- accept[1] + 1  
            }else
            {
                proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
            }
        }        
        accept[2] <- accept[2] + n.beta.block
        regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=TRUE)   
        
        
        
        
        ####################
        ## Sample from phi
        ####################
        chol.Sigma <- t(chol(proposal.sd.phi * Sigma))
        rand <- matrix(rnorm(n=N.all), nrow=K)
        den.offset <- rho * W.triplet.sum + 1 - rho
        phi.offset <- regression.mat + offset.mat
        temp1 <- binomialmcarupdate(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, Y.mat, failures.mat, phi.offset, den.offset, Sigma, Sigma.inv, rho,  chol.Sigma, rand)      
        phi.mat <- temp1[[1]] - mean(temp1[[1]])
        phi <- as.numeric(t(phi.mat))
        accept[3] <- accept[3] + temp1[[2]]
        accept[4] <- accept[4] + K    
        
         
        #####################
        ## Samples from Sigma
        #####################
        Sigma.post.scale <- t(phi.mat) %*% Q.spat %*% phi.mat + prior.Sigma.scale
        Sigma <- riwish(Sigma.post.df, Sigma.post.scale)
        Sigma.inv <- solve(Sigma)

        
        
        ##################
        ## Sample from rho
        ##################
        ## Propose a new value
        proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)
        Q.spat.prop <- proposal.rho * (diag(apply(W, 1, sum)) - W) + diag(rep(1-proposal.rho), K)
        det.Q.spat.prop <-  sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))    
        
        ## Compute the acceptance rate
        logprob.current <- det.Q.spat - 0.5 * sum(diag(t(phi.mat) %*% Q.spat %*% phi.mat %*% Sigma.inv))
        logprob.proposal <- det.Q.spat.prop - 0.5 * sum(diag(t(phi.mat) %*% Q.spat.prop %*% phi.mat %*% Sigma.inv))
        prob <- exp(logprob.proposal - logprob.current)
        if(prob > runif(1))
        {
            rho <- proposal.rho
            det.Q.spat <- det.Q.spat.prop
            Q.spat <- Q.spat.prop
            accept[5] <- accept[5] + 1           
        }else
        {
        }              
        accept[6] <- accept[6] + 1       
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        lp <- offset +  X.standardised %*% beta + phi
        prob <- exp(lp) / (1+exp(lp))
        fitted <- trials * prob
        deviance.all <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)
        deviance <- -2 * sum(deviance.all) 
        
        
        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.beta[ele, ] <- beta
            samples.phi[ele, ] <- phi
            samples.rho[ele, ] <- rho
            samples.Sigma[ele, ] <- as.numeric(Sigma)
            samples.deviance[ele, ] <- deviance
            samples.fitted[ele, ] <- fitted
        }else
        {
        }
        
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
        k <- j/100
        if(ceiling(k)==floor(k))
        {
            #### Determine the acceptance probabilities
            accept.beta <- 100 * accept[1] / accept[2]
            accept.phi <- 100 * accept[3] / accept[4]
            accept.rho <- 100 * accept[5] / accept[6]
            accept.all <- accept.all + accept
            accept <- rep(0,6)
            
            #### phi tuning parameter
            if(accept.phi > 50)
            {
                proposal.sd.phi <- 2 * proposal.sd.phi
            }else if(accept.phi < 40)              
            {
                proposal.sd.phi <- 0.5 * proposal.sd.phi
            }else
            {
            }
            
            #### beta tuning parameter
            if(accept.beta > 40)
            {
                proposal.sd.beta <- 2 * proposal.sd.beta
            }else if(accept.beta < 20)              
            {
                proposal.sd.beta <- 0.5 * proposal.sd.beta
            }else
            {
            }
            
            #### rho tuning parameter
            if(accept.rho > 50)
            {
                proposal.sd.rho <- 2 * min(proposal.sd.rho, 0.5)
            }else if(accept.rho < 40)              
            {
                proposal.sd.rho <- 0.5 * proposal.sd.rho
            }else
            {
            }      
        }else
        {   
        }
        
        
        
        ################################       
        ## print progress to the console
        ################################
        if(j %in% percentage.points & verbose)
        {
            setTxtProgressBar(progressBar, j/n.sample)
        }
    }
    
    # end timer
    if(verbose)
    {
        cat("\nSummarising results")
        close(progressBar)
    }else
    {}
    
    
    ###################################
    #### Summarise and save the results 
    ###################################
    ## Compute the acceptance rates
    accept.beta <- 100 * accept.all[1] / accept.all[2]
    accept.phi <- 100 * accept.all[3] / accept.all[4]
    accept.rho <- 100 * accept.all[5] / accept.all[6]
    accept.final <- c(accept.beta, accept.phi, accept.rho)
    names(accept.final) <- c("beta", "phi", "rho")
    
    
    ## Compute information criterion (DIC, DIC3, WAIC)
    median.beta <- apply(samples.beta,2,median)
    median.phi <- apply(samples.phi, 2, median)
    lp.median <- offset + median.phi + X.standardised %*% median.beta   
    median.prob <- exp(lp.median)  / (1 + exp(lp.median))
    fitted.median <- trials * median.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=median.prob, log=TRUE))
    p.d <- median(samples.deviance) - deviance.fitted
    DIC <- 2 * median(samples.deviance) - deviance.fitted     
    
    ## Compute the LMPL
    CPO <- rep(NA, N.all)
    for(j in 1:N.all)
    {
        CPO[j] <- 1/median((1 / dbinom(x=Y[j], size=trials[j], prob=(samples.fitted[ ,j] / trials[j]))))    
    }
    LMPL <- sum(log(CPO))  
    
    
    ## Create the Fitted values
    fitted.values <- apply(samples.fitted, 2, median)
    residuals <- as.numeric(Y) - fitted.values
    
    
    #### transform the parameters back to the origianl covariate scale.
    samples.beta.orig <- samples.beta
    number.cts <- sum(X.indicator==1)     
    if(number.cts>0)
    {
        for(r in 1:p)
        {
            if(X.indicator[r]==1)
            {
                samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
            }else if(X.indicator[r]==2 & p>1)
            {
                X.transformed <- which(X.indicator==1)
                samples.temp <- as.matrix(samples.beta[ ,X.transformed])
                for(s in 1:length(X.transformed))
                {
                    samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
                }
                intercept.adjustment <- apply(samples.temp, 1,sum) 
                samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
            }else
            {
            }
        }
    }else
    {
    }
    
    
    #### Create a summary object
    samples.beta.orig <- mcmc(samples.beta.orig)
    summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
    summary.hyper <- array(NA, c(4, 7))     
    summary.hyper[1,1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[2,1:3] <- quantile(samples.Sigma[ ,1], c(0.5, 0.025, 0.975))
    summary.hyper[3,1:3] <- quantile(samples.Sigma[ ,4], c(0.5, 0.025, 0.975))
    summary.hyper[4,1:3] <- quantile(samples.Sigma[ ,2], c(0.5, 0.025, 0.975))
    rownames(summary.hyper) <- c("rho", "Sigma11", "Sigma22", "Sigma12")     
    summary.hyper[1, 4:7] <- c(n.keep, accept.rho, effectiveSize(mcmc(samples.rho)), geweke.diag(mcmc(samples.rho))$z)   
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.Sigma[ ,1])), geweke.diag(mcmc(samples.Sigma[ ,1]))$z)   
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.Sigma[ ,4])), geweke.diag(mcmc(samples.Sigma[ ,4]))$z)   
    summary.hyper[4, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.Sigma[ ,2])), geweke.diag(mcmc(samples.Sigma[ ,2]))$z)   
    
    
    summary.results <- rbind(summary.beta, summary.hyper)
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
    ## Compile and return the results
    modelfit <- c(DIC, p.d, LMPL)
    names(modelfit) <- c("DIC", "p.d", "LMPL")
    samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  rho=mcmc(samples.rho), Sigma=mcmc(samples.Sigma), fitted=mcmc(samples.fitted))
    model.string <- c("Likelihood model - binomial (logit link function)", "\nLatent structure model - Autoregressive MCAR model\n")
    results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  X=X)
    class(results) <- "carbayesMV"
    if(verbose)
    {
        b<-proc.time()
        cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
    return(results)
}
