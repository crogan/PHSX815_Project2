#include <iostream>
#include <fstream>

#include "ProbDist.hh"

bool g_model0 = true;

Poisson g_poisson;
Gamma   g_gamma;
MyPrior g_prior;

int GetModelObservation();

using namespace std;

// main function for compiled executable ModelSim.x
int main(int argc, char* argv[]){
  bool printhelp = false;
  long seed = 5555;

  // Number of time measurements (number of missing cookies) - per experiment
  int Nmeas = 1;

  // Number of experiments
  int Nexp = 1;

  // Number of MCMC events to "burn-in" (generate and discard)
  int Nburn = 0;

  // Number of MCMC events to skip in between keeping events
  int Nskip = 0;

  bool doOutputFile = false;
  string OutputFileName;
  

  // parsing any command line options added by the user
  for(int i = 0; i < argc; i++){
    if(strncmp(argv[i],"--help", 6) == 0){
      printhelp = true;
    }
    if(strncmp(argv[i],"-h", 2) == 0){
      printhelp = true;
    }
    if(strncmp(argv[i],"-seed", 5) == 0){
      i++;
      seed = std::stol(argv[i]);
    }
    if(strncmp(argv[i],"--model0", 8) == 0){
      g_model0 = true;
    }
    if(strncmp(argv[i],"--model1", 8) == 0){
      g_model0 = false;
    }
    if(strncmp(argv[i],"-Nmeas", 6) == 0){
      i++;
      int Nt = std::stoi(argv[i]);
      if(Nt > 0)
	Nmeas = Nt;
    }
    if(strncmp(argv[i],"-Nexp", 4) == 0){
      i++;
      int Ne = std::stoi(argv[i]);
      if(Ne > 0)
	Nexp = Ne;
    }
    if(strncmp(argv[i],"-Nburn", 6) == 0){
      i++;
      int Nb = std::stoi(argv[i]);
      if(Nb > 0)
	Nburn = Nb;
    }
    if(strncmp(argv[i],"-Nskip", 6) == 0){
      i++;
      int Ns = std::stoi(argv[i]);
      if(Ns > 0)
	Nskip = Ns;
    }
    if(strncmp(argv[i],"-output", 7) == 0){
      i++;
      OutputFileName = string(argv[i]);
      doOutputFile = true;
    }
  }

  // print the executable usage options if the user adds -h or --help
  if(printhelp){
    cout << "Usage: " << argv[0] << " [options]" << endl;
    cout << "  options:" << endl;
    cout << "   --help(-h)          print options" << endl;
    cout << "   -seed [number]      random seed to use" << endl;
    cout << "   --model0            simulate model 0 (null)" << endl;
    cout << "   --model1            simulate model 1 (aliens)" << endl;
    cout << "   -Nmeas [number]     number of measurements per experiment" << endl;
    cout << "   -Nexp [number]      number of experiments" << endl;
    cout << "   -Nburn [number]     number of MCMC events to discard initially" << endl;
    cout << "   -Nskip [number]     number of MCMC events to skip in between acceptance" << endl;
    cout << "   -output [filename]  name of ouptut file" << endl;
   
    return 0;
  }

  // MCMC burn-in
  for(int i = 0; i < Nburn; i++)
    GetModelObservation();
  
  if(doOutputFile){ // write experiments to file
    ofstream outfile;
    outfile.open(OutputFileName.c_str());
    for(int e = 0; e < Nexp; e++){
      for(int t = 0; t < Nmeas; t++){
	for(int i = 0; i < Nskip; i++)
	  GetModelObservation();
	outfile << GetModelObservation() << " ";
      }
      outfile << endl;
    }
    outfile.close();
  } else { // write experiments to stdout
    for(int e = 0; e < Nexp; e++){
      for(int t = 0; t < Nmeas; t++){
	for(int i = 0; i < Nskip; i++)
	  GetModelObservation();
	cout << GetModelObservation() << " ";
      }
      cout << endl;
    }
  }
  
}

int GetModelObservation(){
  // Generate an observation from our models;

  double x[2];
  double param[2];
  
  // parameters for model 0
  double alpha = 2.;
  double beta  = 1.;

  // if model 1 get random alpha/beta from prior
  if(!g_model0){
    param[0] = 3.;
    param[1] = 4.;
    g_prior.Rand(x, param);
    alpha = x[0];
    beta  = x[1];
  }

  // sample a random weight parameter from
  // Gamma distribution
  param[0] = alpha;
  param[1] = beta;
  g_gamma.Rand(x, param);

  double rate = x[0];

  // sample a random cookie count from Poisson with rate
  param[0] = rate;
  g_poisson.Rand(x, param);

  return x[0];
}
