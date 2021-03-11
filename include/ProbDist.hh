#ifndef ProbDist_HH
#define ProbDist_HH

// C++ input/output includes
#include <iostream>
#include <fstream>

// ROOT includes (special functions)
#include "TMath.h"

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TRandom.h"

// our local includes
#include "Random.hh"

using namespace std;

/////////////////////////////
// ProbDist class
/////////////////////////////
// abstract base class for
// probability distributions
/////////////////////////////
class ProbDist {

public:
  // generic constructor for ProbDist class
  ProbDist();

  // destructor for ProbDist
  virtual ~ProbDist();

  // probability distribution for point x (pointer to an array)
  // given some parameters;
  virtual double Prob(const double* x, const double* param) = 0;

  // generate random point x (pointer to an array) and fill them
  // into array (better have the right dimension declared...)
  virtual void Rand(double* x, const double* param) = 0;

protected:
  // random number generator for class instance
  Random m_random;
  
};

/////////////////////////////
// Gaussian class
/////////////////////////////
// derived class of ProbDist
// implements a Gaussian function
// (samples using rejection sampling)
/////////////////////////////
class Gaussian : public ProbDist {

public:
  // constructor for Gaussian class
  // sigma_match: # sigma to match slope of
  // exponential (for generating random
  // samples with rejection sampling)
  Gaussian(double sigma_match = 2.);

  // destructor for ProbDist
  virtual ~Gaussian();

  // probability distribution for point x
  // x: pointer to 1D array which will be filled with
  //    generated value
  // param: pinter to 2D array with mean and sigma
  double Prob(const double* x, const double* param);

  // sample random value x
  // x: pointer to 1D array (or value) to be filled with random
  // param: pinter to 2D array with mean and sigma
  void Rand(double* x, const double* param);

private:
  // location for matching exponential with slope of Gaussian
  // (# of sigma) for piecewise sampling function in rejection sampling
  double m_match;
  // fraction of total probability in exponential tail from 
  // sampling function (for rejection sampling)
  double m_tailprob;
  
};

/////////////////////////////
// Binomial class
/////////////////////////////
// derived class of ProbDist
// implements a shifted Binomial distribution
// (samples using pre-calculated categorical distribution)
/////////////////////////////
class Binomial : public ProbDist {

public:
  // constructor for Binomial class
  // Nmax: maximum number of steps in either direction
  Binomial(int Nmax = 5);

  // destructor for Binomial
  virtual ~Binomial();

  // probability distribution for point x
  // x: pointer to 1D array which will be filled with
  //    generated value
  // param: ignored
  double Prob(const double* x, const double* param);

  // sample random value x
  // x: pointer to 1D array (or value) to be filled with random
  // param: ignored
  void Rand(double* x, const double* param);

private:
  // maximum number of steps in either direction
  int m_Nmax;
  // pre-calculated probabilities
  vector<double> m_prob;
  // pre-calculated cumlative distribution function
  vector<double> m_cdf;
  
};

/////////////////////////////
// Poisson class
/////////////////////////////
// derived class of ProbDist
// implements a Poisson distribution
// (samples using Metropolis-Hastings)
/////////////////////////////
class Poisson : public ProbDist {

public:
  // constructor for Poisson class
  // Nmax: maximum number of steps in either direction
  //       in transition step
  Poisson(int Nmax = 5);

  // destructor for Poisson
  virtual ~Poisson();

  // probability distribution for point x
  // x: pointer to 1D array which will be filled with
  //    generated value
  // param: pointer to 1D array with rate parameter
  double Prob(const double* x, const double* param);

  // sample random value x
  // x: pointer to 1D array (or value) to be filled with random
  // param: pointer to 1D array with rate parameter
  void Rand(double* x, const double* param);

private:
  Binomial* m_binomial;

  int m_Nprev;
  
};

/////////////////////////////
// Gamma class
/////////////////////////////
// derived class of ProbDist
// implements a Gamma distribution
// (samples using Metropolis-Hastings)
/////////////////////////////
class Gamma : public ProbDist {

public:
  // constructor for Gamma class
  // sigma: width of Gaussian used for
  //        transition distribution
  Gamma(int sigma = 1);

  // destructor for Gamma
  virtual ~Gamma();

  // probability distribution for point x
  // x: pointer to 1D array which will be filled with
  //    generated value
  // param: pointer to 2D array with alpha (shape) and beta (rate)
  double Prob(const double* x, const double* param);

  // sample random value x
  // x: pointer to 1D array (or value) to be filled with random
  // param: pointer to 2D array with alpha (shape) and beta (rate)
  void Rand(double* x, const double* param);

private:
  double m_sigma;
  Gaussian* m_gaussian;

  double m_Xprev;
  
};

/////////////////////////////
// MyPrior class
/////////////////////////////
// derived class of ProbDist
// implements a custom prior distribution
// for the parameters of a Gamma function
// (samples using Metropolis-Hastings)
/////////////////////////////
class MyPrior : public ProbDist {

public:
  // constructor for MyPrior class
  // sigma: width of Gaussian used for
  //        transition distribution
  MyPrior(int sigma = 1);

  // destructor for MyPrior
  virtual ~MyPrior();

  // probability distribution for point x
  // x: pointer to 2D array which will be filled with
  //    generated values of alpha and beta
  // param: pointer to 2D array with alpha0 (shape) and beta0 (rate)
  double Prob(const double* x, const double* param);

  // sample random value x
  // x: pointer to 2D array which will be filled with
  //    generated values of alpha and beta
  // param: pointer to 2D array with alpha0 (shape) and beta0 (rate)
  void Rand(double* x, const double* param);

private:
  double m_sigma;
  Gaussian* m_gaussian;

  double m_alpha_prev;
  double m_beta_prev;
  
};



#endif
