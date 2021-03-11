// C++ input/output includes
#include <iostream>
#include <fstream>

// ROOT includes (for special functions)
#include "TMath.h"

// our local includes
#include "ProbDist.hh"

// generic constructor for ProbDist class
ProbDist::ProbDist(){

}

// destructor for ProbDist
ProbDist::~ProbDist(){

}

// constructor for Gaussian class
// sigma_match: # sigma to match slope of
// exponential (for generating random
// samples with rejection sampling)
Gaussian::Gaussian(double sigma_match)
  : ProbDist() {
  if(sigma_match > 0)
    m_match = sigma_match;
  else
    m_match = 1;

  m_tailprob = 1./m_match*exp(-m_match*m_match/2.);
  m_tailprob = m_tailprob / (2.*m_tailprob + 2*m_match);
}

// destructor for Gaussian class
Gaussian::~Gaussian(){

}

// probability distribution for point x
// x: pointer to 1D array which will be filled with
//    generated value
// param: pinter to 2D array with mean and sigma
double Gaussian::Prob(const double* x, const double* param){
  double mean  = param[0];
  double sigma = param[1];
  if(sigma <= 0)
    sigma = 1.;
      
  return 1./sqrt(2.*acos(-1))/sigma*exp(-(x[0]-mean)*(x[0]-mean)/2./sigma/sigma);
}

// sample random value x
// x: pointer to 1D array (or value) to be filled with random
// param: pinter to 2D array with mean and sigma
void Gaussian::Rand(double* x, const double* param){
  double mean  = param[0];
  double sigma = param[1];
  if(sigma <= 0)
    sigma = 1.;

  bool accept = false;

  // put return value here
  double X;
  while(!accept){
    // Get a random number between [0,1]
    // to decide which piece (expo or flat)
    // to sample from
    double R1 = m_random.rand();
    
    // Get a 2nd random number between [0,1]
    // to decide value from segment
    double R2 = m_random.rand();
    
    double ratio;
    
    if(R1 < m_tailprob){ // get from lower exponential
      X = -m_match + log(R2)/m_match;
      ratio = exp(-X*X/2.) / (exp(m_match*m_match/2.)*exp(X*m_match));
    } else if(R1 > 1-m_tailprob) { // get from upper exponential
      X = m_match - log(R2)/m_match;
      ratio = exp(-X*X/2.) / (exp(m_match*m_match/2.)*exp(-X*m_match));
    } else {
      X = -m_match + R2*2.*m_match;
      ratio = exp(-X*X/2.);
    }

    // Get a 3nd random number between [0,1]
    // to decide whether to accept X
    double R3 = m_random.rand();
    if(R3 <= ratio)
      accept = true;
  }

  x[0] = mean + sigma*X;
}

// constructor for Binomial class
// Nmax: maximum number of steps in either direction
Binomial::Binomial(int Nmax)
  : ProbDist() {
  if(Nmax > 0)
    m_Nmax = Nmax;
  else
    m_Nmax = 1;

  double integral = 0;
  for(int n = 0; n < 2*m_Nmax+1; n++){
    m_prob.push_back(TMath::Binomial(2*m_Nmax, n));
    integral += m_prob[n];
    m_cdf.push_back(integral);
  }
  for(int n = 0; n < 2*m_Nmax+1; n++){
    m_prob[n] /= integral;
    m_cdf[n] /= integral;
  }
}

// destructor for Binomial
Binomial::~Binomial(){

}

// probability distribution for point x
// x: pointer to 1D array which will be filled with
//    generated value
// param: ignored
double Binomial::Prob(const double* x, const double* param){
  if(abs(int(x[0])) > m_Nmax)
    return 0;

  return m_prob[int(x[0])+m_Nmax];
}

// sample random value x
// x: pointer to 1D array (or value) to be filled with random
// param: ignored
void Binomial::Rand(double* x, const double* param){
  // get a random number [0,1] and see which cdf bin the value falls in
  double R = m_random.rand();

  for(int i = 0; i < 2*m_Nmax+1; i++){
    if(R <= m_cdf[i] || i == 2*m_Nmax){
      x[0] = i - m_Nmax;
      break;
    }
  }
}

// constructor for Poisson class
// Nmax: maximum number of steps in either direction
//       in transition step
Poisson::Poisson(int Nmax)
  : ProbDist() {
  m_Nprev = 0; // start at zero for previous sample

  m_binomial = new Binomial(Nmax);
}

// destructor for Poisson
Poisson::~Poisson(){
  if(m_binomial)
    delete m_binomial;
}

// probability distribution for point x
// x: pointer to 1D array which will be filled with
//    generated value
// param: pointer to 1D array with rate parameter
double Poisson::Prob(const double* x, const double* param){
  double rate = param[0];
  if(rate <= 0)
    rate = 1.;

  return TMath::Poisson(x[0], rate);
}

// sample random value x
// x: pointer to 1D array (or value) to be filled with random
// param: pointer to 1D array with rate parameter
void Poisson::Rand(double* x, const double* param){

  int Nnew = -1;

  // get a new trial value of N using binomial transition
  while(Nnew < 0){
    m_binomial->Rand(x, param);
    Nnew = m_Nprev + x[0];
  }

  // calculate Poisson prob of new value
  x[0] = Nnew;
  double Ppois_new = Prob(x, param);
  // calculate prob of old -> new transition
  x[0] = Nnew - m_Nprev;
  double Psample_new = m_binomial->Prob(x, param);
  // calculate Poisson prob of old value
  x[0] = m_Nprev;
  double Ppois_old = Prob(x, param);
  // calculate prob of new -> old transition
  x[0] = m_Nprev - Nnew;
  double Psample_old = m_binomial->Prob(x, param);
  
  // calculate acceptance probability
  double A = std::min(1., Ppois_new*Psample_old / (Ppois_old*Psample_new));

  // get random number [0,1]
  double R = m_random.rand();
  if(R <= A){
    x[0] = Nnew;
    m_Nprev = Nnew;
  } else {
    x[0] = m_Nprev;
  }
}

// constructor for Gamma class
// sigma: width of Gaussian used for
//        transition distribution
Gamma::Gamma(int sigma)
  : ProbDist() {
  m_Xprev = 1; // start at 1 for previous sample

  m_sigma = sigma;
  m_gaussian = new Gaussian();
}

// destructor for Gamma
Gamma::~Gamma(){
  if(m_gaussian)
    delete m_gaussian;
}

// probability distribution for point x
// x: pointer to 1D array which will be filled with
//    generated value
// param: pointer to 2D array with alpha (shape) and beta (rate)
double Gamma::Prob(const double* x, const double* param){
  double alpha = param[0];
  double beta  = param[1];

  if(alpha <= 0)
    alpha = 1.;

  if(beta <= 0)
    beta = 1.;

  return pow(beta, alpha)/TMath::Gamma(alpha)*pow(x[0], alpha-1.)*exp(-beta*x[0]);
}

// sample random value x
// x: pointer to 1D array (or value) to be filled with random
// param: pointer to 2D array with alpha (shape) and beta (rate)
void Gamma::Rand(double* x, const double* param){
  double Xnew = -1;

  double param_gaus[2];
  param_gaus[0] = 0; // Gaussian mean;
  param_gaus[1] = m_sigma; // Gaussian sigma
  
    // get a new trial value of X using gaussian transition
    while(Xnew < 0){
      m_gaussian->Rand(x, param_gaus);
      Xnew = m_Xprev + x[0];
    }

  // calculate Gamma prob of new value
  x[0] = Xnew;
  double Pgam_new = Prob(x, param);
  // calculate prob of old -> new transition
  x[0] = Xnew - m_Xprev;
  double Psample_new = m_gaussian->Prob(x, param_gaus);
  // calculate Gamma prob of old value
  x[0] = m_Xprev;
  double Pgam_old = Prob(x, param);
  // calculate prob of new -> old transition
  x[0] = m_Xprev - Xnew;
  double Psample_old = m_gaussian->Prob(x, param_gaus);
  
  // calculate acceptance probability
  double A = std::min(1., Pgam_new*Psample_old / (Pgam_old*Psample_new));

  // get random number [0,1]
  double R = m_random.rand();
  if(R <= A){
    x[0] = Xnew;
    m_Xprev = Xnew;
  } else {
    x[0] = m_Xprev;
  }
}

// constructor for MyPrior class
// sigma: width of Gaussian used for
//        transition distribution
MyPrior::MyPrior(int sigma)
  : ProbDist() {

  m_sigma = sigma;
  m_gaussian = new Gaussian();
  
  // initial parameter values
  m_alpha_prev = 1.;
  m_beta_prev = 1.;
}
  
// destructor for MyPrior
MyPrior::~MyPrior(){
  if(m_gaussian)
    delete m_gaussian;
}

// probability distribution for point x
// x: pointer to 2D array which will be filled with
//    generated values of alpha and beta
// param: pointer to 2D array with alpha0 (shape) and beta0 (rate)
// NOTE: this probability is not normalized!
double MyPrior::Prob(const double* x, const double* param){
  double Na = 3.;
  double Nb = 1.;
  double a = x[0];
  double b = x[1];

  double alpha0 = param[0];
  double beta0  = param[1];

  if(alpha0 <= 0)
    alpha0 = 1;
  if(beta0 <= 0)
    beta0 = 1;

  return pow(alpha0, Na*(a-1.))*exp(-Nb*b*beta0)/pow(TMath::Gamma(a), Na)/pow(b, -Nb*a);
  
}

// sample random value x
// x: pointer to 2D array which will be filled with
//    generated values of alpha and beta
// param: pointer to 2D array with alpha0 (shape) and beta0 (rate)
void MyPrior::Rand(double* x, const double* param){
  double alpha_new = -1;
  double beta_new  = -1;

  double param_gaus[2];
  param_gaus[0] = 0; // Gaussian mean;
  param_gaus[1] = m_sigma; // Gaussian sigma
  
  // get a new trial value of alpha using gaussian transition
  while(alpha_new < 0){
    m_gaussian->Rand(x, param_gaus);
    alpha_new = m_alpha_prev + x[0];
  }
  // get a new trial value of beta using gaussian transition
  while(beta_new < 0){
    m_gaussian->Rand(x, param_gaus);
    beta_new = m_beta_prev + x[0];
  }

  // calculate prob of old -> new transition
  x[0] = alpha_new - m_alpha_prev;
  double Psample_new = m_gaussian->Prob(x, param_gaus);
  x[0] = beta_new - m_beta_prev;
  Psample_new *= m_gaussian->Prob(x, param_gaus);

  // calculate prob of new -> old transition
  x[0] = m_alpha_prev - alpha_new;
  double Psample_old = m_gaussian->Prob(x, param_gaus);
  x[0] =  m_beta_prev - beta_new;
  Psample_old *= m_gaussian->Prob(x, param_gaus);
  
  // calculate prob of new values for our function
  x[0] = alpha_new;
  x[1] = beta_new;
  double Pfunc_new = Prob(x, param);
  
  // calculate prob of old values for our function
  x[0] = m_alpha_prev;
  x[1] = m_beta_prev;
  double Pfunc_old = Prob(x, param);
  
  // calculate acceptance probability
  double A = std::min(1., Pfunc_new*Psample_old / (Pfunc_old*Psample_new));

  // get random number [0,1]
  double R = m_random.rand();
  if(R <= A){
    x[0] = alpha_new;
    x[1] = beta_new;
    m_alpha_prev = alpha_new;
    m_beta_prev  = beta_new;
  } else {
    x[0] = m_alpha_prev;
    x[1] = m_beta_prev;
  }
}

