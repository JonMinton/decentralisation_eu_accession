#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, 
                          NumericVector beta, NumericVector offset)
{
//Create new objects
// Compute the linear predictor
NumericVector linpred(nsites);
double temp; 


//  Compute the linear predictor via a double for loop
     for(int j = 0; j < nsites; j++)
     {
     temp = 0;
      
          for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];     
          
     linpred[j] = temp + offset[j];  
     }


// Return the result
return linpred;
}




// [[Rcpp::export]]
double binomialbetaupdate(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                          NumericVector proposal, NumericVector offset, NumericVector y, 
                          NumericVector failures, NumericVector prior_meanbeta,
                          NumericVector prior_varbeta)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
    
    
    // Create the log likelihood acceptance probability component
    lp_current = linpredcompute(X, nsites, p, beta, offset);
    lp_proposal = linpredcompute(X, nsites, p, proposal, offset);     
    for(int j = 0; j < nsites; j++)     
    {
        p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
        p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
        oldlikebit = oldlikebit + y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
        newlikebit = newlikebit + y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    
    // Create the prior acceptance component
    for(int j = 0; j < p; j++)     
    {
        priorbit = priorbit + 0.5 * pow((beta[j]-prior_meanbeta[j]),2) / prior_varbeta[j] - 0.5 * pow((proposal[j]-prior_meanbeta[j]),2) / prior_varbeta[j];
    }
    
    
    // Compute the acceptance probability and return the value
    acceptance = exp(likebit + priorbit);
    return acceptance;
}


// [[Rcpp::export]]
List binomialmcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum,const int nsites,  const int nvar, NumericMatrix phi, 
     NumericMatrix Y, NumericMatrix failures, NumericMatrix phioffset, 
     NumericVector denoffset, NumericMatrix Sigma, NumericMatrix Sigmainv, double rho, 
     NumericMatrix cholsigma, NumericMatrix rand)
{
// Update the spatially correlated random effects 
//Create new objects
NumericMatrix phinew(nsites, nvar);
NumericMatrix fcvar(nvar, nvar), fcprec(nvar, nvar);
NumericVector sumphi(nvar), priormean(nvar), propphi(nvar);
NumericVector diff(nvar), diffprop(nvar);
NumericVector pold(nvar), pnew(nvar);
double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
double priorvardenom;
int rowstart=0, rowend=0, accept=0;

//  Update each random effect in turn
     for(int j = 0; j < nsites; j++)
     {
     // Calculate prior variance and precision
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     fcvar(_,0) = Sigma(_,0) / priorvardenom;
     fcvar(_,1) = Sigma(_,1) / priorvardenom;
     fcprec(_,0) = Sigmainv(_,0) * priorvardenom;
     fcprec(_,1) = Sigmainv(_,1) * priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = rep(0,2);
     for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
     priormean = rho * sumphi / priorvardenom; 
      
      // propose a value  
      propphi[0] = phi(j, 0) + sum(cholsigma(0, _) * rand(j,_));
      propphi[1] = phi(j, 1) + sum(cholsigma(1, _) * rand(j,_));

      
      // Accept or reject it
      diff = priormean - phi(j,_);
      diffprop = priormean - propphi;
      oldpriorbit = 0.5 * (diff[0] * diff[0] * fcprec(0,0) + 2 * diff[0] * diff[1] * fcprec(0,1) + diff[1] * diff[1] * fcprec(1,1));
      newpriorbit = 0.5 * (diffprop[0] * diffprop[0] * fcprec(0,0) + 2 * diffprop[0] * diffprop[1] * fcprec(0,1) + diffprop[1] * diffprop[1] * fcprec(1,1));
      pold = exp(phioffset(j,_) + phi(j,_)) / (1 + exp(phioffset(j,_) + phi(j,_)));
      pnew = exp(phioffset(j,_) + propphi) / (1 + exp(phioffset(j,_) + propphi));
      oldlikebit = sum(Y(j,_) * log(pold) + failures(j,_) * log(1 - pold));
      newlikebit = sum(Y(j,_) * log(pnew) + failures(j,_) * log(1 - pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            if(runif(1)[0] <= acceptance) 
            {
            phi(j,_) = propphi;
            accept = accept + 1;
            }
            else
            { 
            }
    }


List out(2);
out[0] = phi;
out[1] = accept;
return out;
}
