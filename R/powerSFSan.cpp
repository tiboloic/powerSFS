// LT 19/08/2020

// version for variable allele numbers
// use asymptotic approximation 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // DATA
  
  // observations: counts per category per gene
  DATA_MATRIX(obs);
  
  // Average allele number for each gene
  DATA_IVECTOR(an);
  
  // average slope beta
  PARAMETER(b);
  
  //log-stdev of random effect on slope
  PARAMETER(lsigb);
  
  // random effects on slope
  PARAMETER_VECTOR(bdevs);
  
  // Bernouilli polynomials at 0
  //const Type Bs[8] = {Type(1/6), Type(-1/30), Type(1/42), Type(-1/30), Type(5/66), Type(-691/2730), Type(7/6), Type(-3617/510)};
  //const Type Bs[8] = {1.0/6, -1.0/30, 1.0/42, -1.0/30, 5.0/66, -691.0/2730, 7.0/6, -3617.0/510};
  //Type Bs[8];
  
  vector<Type> Bs(8);Bs << 1.0/6, -1.0/30, 1.0/42, -1.0/30, 5.0/66, -691.0/2730, 7.0/6, -3617.0/510;
  //Bs(0) = 1/6; Bs(1) = -1/30; Bs(2) = 1/42; Bs(3) = -1/30; Bs(4) = 5/66; Bs(5) = -691/2730; Bs(6) = 7/6;Bs(7) = -3617/510;
  //REPORT(Bs[0]);REPORT(Bs[1]);REPORT(Bs[2]);REPORT(Bs[3]);REPORT(Bs[4]);REPORT(Bs[5]);
  //REPORT(Bs);
  
  
  Type nll=0.0;
  
  Type beta;
  
 /* Type zeta(Type n, Type beta) {
    // sum_{i=n}^\infty 1/i^beta
    // using asymptotic approximation
    Type ans = pow(n,1-beta)/(beta  - 1) + 0.5 / pow(n, beta);
    Type term = beta/2 / pow(n,beta-1);
    ans = ans + term * Bs[1];
    
    for (int i=2;i<9;i++) {
      term = term * (beta + 2*i-3) * (beta + 2*i -2) / ((2*i -1) * 2*i * pow(n, 2));
      ans = ans + term * Bs[i];
    }
    return(ans);
  }
*/
  //parallel_accumulator<Type> nll(this);
  
  // add likelihood contribution for gene level random effects
  nll-= dnorm(bdevs, Type(0), Type(1), true).sum();
  
  // likelihood contribution per gene
  for(int i=0;i<obs.rows();i++){
    

    // build beta
    beta = b + exp(lsigb) * bdevs(i);
    
    // vector of probabilities
    vector<Type> ps(obs.cols());
    //fill(ps.begin(), ps.end(), Type(0.0));
    // initialize to 0
    for (int j=0; j<obs.cols(); j++) {
      //ps(j) = Type(0.0);
      ps(j) = Type(1.0e-8);
    }

    // calculate the first terms up to 15
    // we assume safely that AN is bigger than 15
    for (int j=0; j<15; j++) {
      ps(floor(log2(j+1))) +=  Type(pow(j+1, -beta));
    }
    
    // use asympotic approximation from 16 to AN. 16 = 2^4
    
    // build vector of bin limits on which zeta will be applied
    int neval = floor(log2(an(i))) - 2;
    vector<int> lims(neval);
    for (int j=0; j<neval-1; j++) {
      lims(j) = pow(2,j+4);
    }
    lims(neval-1) = floor(an(i)) + 1;
    REPORT(lims);
        
    // now eval zeta on each element of the vector
    vector<Type> zetas(neval);
    for (int j=0; j<neval; j++) {
      
      //vector<Type> zs(9);
      
      Type n = Type(lims(j));
      Type zeta = pow(n,1-beta)/(beta  - 1) + 0.5 / pow(n, beta);
      
      //zs(0) = zeta;
      
      Type term = beta/2 / pow(n, beta+1);
      zeta = zeta + term * Bs[0];
      
      //zs(1) = zeta;
      
      for (int k=2;k<9;k++) {
        term = term * (beta + 2*k-3) * (beta + 2*k -2) / ((2*k -1) * 2*k * n * n);
        zeta = zeta + term * Bs[k-1];
        //zs(k) = zeta;
      }
      zetas(j) = zeta;
      //REPORT(zs);
    }
    REPORT(zetas);
    REPORT(beta);REPORT(bdevs);REPORT(lsigb);
    
    // get adjacent differences and update probabilities
    for (int j=0; j<neval-1; j++) {
      ps(j+4) = zetas(j) - zetas(j+1);
    }
    
    Type psum = ps.sum();
    ps = ps / psum;
    
    REPORT(ps);
    
    vector<Type> counts = obs.row(i);
    
    // multinomial sampling
    nll -= dmultinom(counts, ps, true);
    
    REPORT(dmultinom(counts, ps, false));
    REPORT(counts);
  }
  return nll;
}
