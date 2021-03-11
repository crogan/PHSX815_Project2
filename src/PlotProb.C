// C++ input/output includes
#include <iostream>
#include <fstream>
#include <string>

// ROOT includes
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom.h"

// our local includes
#include "ProbDist.hh"

// canvas counter
int g_ican = 0;
// set some default colors
int g_color[4] = {kBlue+2,kGreen+2,kRed+2,kMagenta+2}; 

TCanvas* DrawHists(vector<TH1D*>& hists,
		   vector<string>& labels,
		   string title,
		   string xlabel,
		   string ylabel);

TCanvas* Draw2DHist(TH2D* hist,
		    string title,
		    string xlabel,
		    string ylabel,
		    string zlabel);

// main function for compiled executable PlotProb.x
// takes no command line arguments
int main(int argc, char* argv[]){

  // some formating settings
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  string title;
  vector<string> labels;
  string xlabel;
  string ylabel;
  
  ///////////////////////////////////////
  // make Gaussian transition probability
  // plots for different sigma
  ///////////////////////////////////////
  Gaussian GausDist;
  double param_gaus[2];
  vector<double> sigmas = {1, 2, 3};
  double x;
  param_gaus[0] = 0; // mean
  param_gaus[1] = 1; // sigma

  vector<TH1D*> hists_gaus;
  for(int i = 0; i < int(sigmas.size()); i++){
    hists_gaus.push_back(new TH1D(Form("h_gaus_%d",i),
				  Form("h_gaus_%d",i),
				  100, -5., 5.));

    param_gaus[1] = sigmas[i]; // sigma
    for(int j = 0; j < 1000000; j++){
      GausDist.Rand(&x, param_gaus);
      hists_gaus[i]->Fill(x);
    }
    
    hists_gaus[i]->Scale(1./hists_gaus[i]->Integral());
  }

  title = "Gaussian transition probabilities";
  labels.clear();
  for(int i = 0; i < int(sigmas.size()); i++)
    labels.push_back(string(Form("#sigma = %d",int(sigmas[i]))));
  xlabel = "X";
  ylabel = "#propto P(X | #sigma)";
 
  TCanvas* can0 = DrawHists(hists_gaus, labels,
			    title, xlabel, ylabel);
  can0->SaveAs("Gaussians.pdf");

  ///////////////////////////////////////
  // make Binomial transition probability
  // plots for different max steps
  ///////////////////////////////////////
  vector<Binomial*> BinomialDist;
  vector<double> Nmax = {2, 4, 6, 8};
  
  vector<TH1D*> hists_binom;
  for(int i = 0; i < int(Nmax.size()); i++){
    hists_binom.push_back(new TH1D(Form("h_binom_%d",i),
				  Form("h_binom_%d",i),
				  17, -8.5, 8.5));
    BinomialDist.push_back(new Binomial(Nmax[i]));
    for(int j = 0; j < 1000000; j++){
      BinomialDist[i]->Rand(&x, nullptr);
      hists_binom[i]->Fill(x);
    }
    
    hists_binom[i]->Scale(1./hists_binom[i]->Integral());
  }

  title = "Binomial transition probabilities";
  labels.clear();
  for(int i = 0; i < int(Nmax.size()); i++)
    labels.push_back(string(Form("N_{max} = %d",int(Nmax[i]))));
  xlabel = "X";
  ylabel = "#propto P(X | N_{max})";
 
  TCanvas* can1 = DrawHists(hists_binom, labels,
			    title, xlabel, ylabel);
  can1->SaveAs("Binomials.pdf");

  ///////////////////////////////////////
  // make Poisson distribution
  // plots for different rate parameters
  ///////////////////////////////////////
  Poisson PoissonDist;
  vector<double> rate = {1, 2, 4, 6};
  double param_pois;
  vector<TH1D*> hists_pois;
  for(int i = 0; i < int(rate.size()); i++){
    hists_pois.push_back(new TH1D(Form("h_pois_%d",i),
				  Form("h_pois_%d",i),
				  12, -0.5, 11.5));
    param_pois = rate[i];
    for(int j = 0; j < 1000000; j++){
      PoissonDist.Rand(&x, &param_pois);
      hists_pois[i]->Fill(x);
    }
    
    hists_pois[i]->Scale(1./hists_pois[i]->Integral());
  }

  title = "Poisson Distributions";
  labels.clear();
  for(int i = 0; i < int(Nmax.size()); i++)
    labels.push_back(string(Form("#lambda = %d",int(rate[i]))));
  xlabel = "X";
  ylabel = "#propto P(X | #lambda)";
 
  TCanvas* can2 = DrawHists(hists_pois, labels,
			    title, xlabel, ylabel);
  can2->SaveAs("Poissons.pdf");

  ///////////////////////////////////////
  // make Gamma distribution
  // plots for different parameters
  ///////////////////////////////////////
  Gamma GammaDist;
  vector<double> alpha = {1, 1.5, 3.5, 6};
  vector<double> beta = {1., 0.5, 0.8, 0.9};
  double param_gam[2];
  vector<TH1D*> hists_gam;
  for(int i = 0; i < int(alpha.size()); i++){
    hists_gam.push_back(new TH1D(Form("h_gam_%d",i),
				  Form("h_gam_%d",i),
				  100, 0., 12.));
    param_gam[0] = alpha[i];
    param_gam[1] = beta[i];
    for(int j = 0; j < 1000000; j++){
      GammaDist.Rand(&x, param_gam);
      hists_gam[i]->Fill(x);
    }
    
    hists_gam[i]->Scale(1./hists_gam[i]->Integral());
  }

  title = "Gamma Distributions";
  labels.clear();
  for(int i = 0; i < int(alpha.size()); i++)
    labels.push_back(string(Form("#alpha = %0.1f, #beta = %0.1f", alpha[i], beta[i])));
  xlabel = "X";
  ylabel = "#propto P(X | #alpha, #beta)";
 
  TCanvas* can3 = DrawHists(hists_gam, labels,
			    title, xlabel, ylabel);
  can3->SaveAs("Gammas.pdf");

  ///////////////////////////////////////
  // make distribution of alpha/beta Gamma
  // parameters (prior)
  ///////////////////////////////////////
  MyPrior PriorDist;
  double x_prior[2];
  double param_prior[2];
  param_prior[0] = 3.; // alpha0 parameter for func
  param_prior[1] = 4.; // beta0 parameter for func
  TH2D* hist2D = new TH2D("hist2D", "hist2D",
			  50, 0., 8.,
			  50, 0., 4.);
  
  for(int i = 0; i < 1000000; i++){
    PriorDist.Rand(x_prior, param_prior);
    hist2D->Fill(x_prior[0], x_prior[1]);
  }

  // // uncomment below to evaluate PDF analytically
  // for(int ix = 0; ix < 50; ix++){
  //   x_prior[0] = hist2D->GetXaxis()->GetBinCenter(ix+1);
  //   for(int iy = 0; iy < 50; iy++){
  //     x_prior[1] = hist2D->GetYaxis()->GetBinCenter(iy+1);
  //     hist2D->SetBinContent(ix+1,iy+1,PriorDist.Prob(x_prior, param_prior));
  //   }
  // }
  
  hist2D->Scale(1./hist2D->Integral());


  title = "Prior Distribution for #alpha/#beta";
  xlabel = "#alpha";
  ylabel = "#beta";
  string zlabel = "#propto P(#alpha, #beta | #alpha_{0}, #beta_{0})";
 
  TCanvas* can4 = Draw2DHist(hist2D, title,
			     xlabel, ylabel, zlabel);
  can4->SaveAs("Prior2D.pdf");
  
 
}

TCanvas* DrawHists(vector<TH1D*>& hists,
		   vector<string>& labels,
		   string title,
		   string xlabel,
		   string ylabel){
  g_ican++;
  TCanvas* can = (TCanvas*) new TCanvas(Form("can_%d",g_ican),
					Form("can_%d",g_ican),
					450, 400);
  double hlo = 0.15;
  double hhi = 0.08;
  double hbo = 0.15;
  double hto = 0.07;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridx();
  can->SetGridy();
  can->Draw();
  can->cd();

  for(int i = 0; i < int(hists.size()); i++){
    hists[i]->SetLineColor(g_color[i]);
    hists[i]->SetFillColor(g_color[i]-2);
    hists[i]->SetFillStyle(3004);
  }
    
  hists[0]->Draw("hist");
  hists[0]->GetXaxis()->CenterTitle();
  hists[0]->GetXaxis()->SetTitleFont(42);
  hists[0]->GetXaxis()->SetTitleSize(0.05);
  hists[0]->GetXaxis()->SetTitleOffset(1.1);
  hists[0]->GetXaxis()->SetLabelFont(42);
  hists[0]->GetXaxis()->SetLabelSize(0.04);
  hists[0]->GetXaxis()->SetTitle(xlabel.c_str());
  hists[0]->GetXaxis()->SetTickSize(0.);
  hists[0]->GetYaxis()->CenterTitle();
  hists[0]->GetYaxis()->SetTitleFont(42);
  hists[0]->GetYaxis()->SetTitleSize(0.05);
  hists[0]->GetYaxis()->SetTitleOffset(1.2);
  hists[0]->GetYaxis()->SetLabelFont(42);
  hists[0]->GetYaxis()->SetLabelSize(0.035);
  hists[0]->GetYaxis()->SetTitle(ylabel.c_str());

  for(int i = int(hists.size())-1; i >= 0; i--)
    hists[i]->Draw("same hist");

  TLegend* leg = new TLegend(0.62,0.72,0.82,0.925);
  leg->SetTextFont(132);
  leg->SetTextSize(0.045);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int i = 0; i < int(hists.size()); i++)
    leg->AddEntry(hists[i], labels[i].c_str());
  leg->Draw();

  TLatex l;
  l.SetTextFont(42);
  l.SetTextAlign(21);
  l.SetTextSize(0.04);
  l.SetNDC();

  l.DrawLatex((1.-hhi+hlo)/2., 1.-hto+0.012, title.c_str());

  return can;
}

TCanvas* Draw2DHist(TH2D* hist,
		    string title,
		    string xlabel,
		    string ylabel,
		    string zlabel){
  g_ican++;
  TCanvas* can = (TCanvas*) new TCanvas(Form("can_%d",g_ican),
					Form("can_%d",g_ican),
					500, 400);

  double hlo = 0.15;
  double hhi = 0.23;
  double hbo = 0.15;
  double hto = 0.07;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridx();
  can->SetGridy();
  //can->SetLogz();
  can->Draw();
  can->cd();

  hist->Draw("colz");
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitle(xlabel.c_str());
  hist->GetXaxis()->SetTickSize(0.);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTitle(ylabel.c_str());
  hist->GetZaxis()->CenterTitle();
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelSize(0.035);
  hist->GetZaxis()->SetTitle(zlabel.c_str());
  
  TLatex l;
  l.SetTextFont(42);
  l.SetTextAlign(21);
  l.SetTextSize(0.04);
  l.SetNDC();

  l.DrawLatex((1.-hhi+hlo)/2., 1.-hto+0.012, title.c_str());
  
  return can;
}
