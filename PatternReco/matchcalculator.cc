#include "TMath.h"
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGaxis.h"
#include <fstream>
#include <iostream>
#include "TMath.h"
#include "plotHist.cc"
#include <string>
#include <map>



void matchcalculator(){


  gROOT->Reset();
  
  gROOT->SetStyle("Plain");

  gStyle->SetCanvasColor(kWhite);

  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptTitle(1);
  
  // For publishing:
  gStyle->SetLineWidth(1);
  gStyle->SetTextSize(1.1);
  gStyle->SetLabelSize(0.06,"xy");
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.0,"y");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.12);
  
  
  int ncut1=72;
  int ncut2=108;


  TCanvas* c1 = new TCanvas("c1","Track performance",200,10,700,800);
  c1->Divide(2,2);
  c1->SetFillColor(0);
  c1->SetGrid();
  
  TCanvas* c2 = new TCanvas("c2","Track performance",200,10,700,800);
  c2->Divide(2,3);
  c2->SetFillColor(0);
  c2->SetGrid();

  TCanvas* c3 = new TCanvas("c3","Track performance",200,10,700,800);
  c3->Divide(2,3);
  c3->SetFillColor(0);
  c3->SetGrid();

  double max=128.0;
 
  TH1 *hist1 = new TH1F("h1","MatchesAll Layer",max,0.5,max+0.5);
  TH1 *hist2 = new TH1F("h2","MatchesPass Layer",max,0.5,max+0.5);

  TH1 *hist3 = new TH1F("h3","MatchesAll Disk",max,0.5,max+0.5);
  TH1 *hist4 = new TH1F("h4","MatchesPass Disk",max,0.5,max+0.5);

  TH1 *hist_L1 = new TH1F("L1","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L1pass = new TH1F("L1","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L2 = new TH1F("L2","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L2pass = new TH1F("L2","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L3 = new TH1F("L3","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L3pass = new TH1F("L3","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L4 = new TH1F("L4","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L4pass = new TH1F("L4","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L5 = new TH1F("L5","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L5pass = new TH1F("L5","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L6 = new TH1F("L6","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_L6pass = new TH1F("L6","Matches per VM",max,0.5,max+0.5);

  TH1 *hist_D1 = new TH1F("D1","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D1pass = new TH1F("D1","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D2 = new TH1F("D2","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D2pass = new TH1F("D2","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D3 = new TH1F("D3","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D3pass = new TH1F("D3","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D4 = new TH1F("D4","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D4pass = new TH1F("D4","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D5 = new TH1F("D5","Matches per VM",max,0.5,max+0.5);
  TH1 *hist_D5pass = new TH1F("D5","Matches per VM",max,0.5,max+0.5);
  
  ifstream in("matchcalculator.txt");

  int count=0;

  std::map<TString, std::pair<TH1*,TH1*> > hists;

  while (in.good()){

    TString name;
    int matchesAll,matchesPass;
  
    in >>name>>matchesAll>>matchesPass;

    if (!in.good()) continue;
    //cout << name <<" "<< countall <<" "<< countpass << endl;

    if (matchesAll>max) matchesAll=max;
    if (matchesPass>max) matchesPass=max;

    if (name[3]=='L'){
      hist1->Fill(matchesAll);
      hist2->Fill(matchesPass);
    }
    
    if (name[3]=='D'){
      hist3->Fill(matchesAll);
      hist4->Fill(matchesPass);
    }

    if (name[3]=='L' && name[4]=='1'){
      hist_L1->Fill(matchesAll);
      hist_L1pass->Fill(matchesPass);
    }

    if (name[3]=='L' && name[4]=='2'){
      hist_L2->Fill(matchesAll);
      hist_L2pass->Fill(matchesPass);
    }

    if (name[3]=='L' && name[4]=='3'){
      hist_L3->Fill(matchesAll);
      hist_L3pass->Fill(matchesPass);
    }

    if (name[3]=='L' && name[4]=='4'){
      hist_L4->Fill(matchesAll);
      hist_L4pass->Fill(matchesPass);
    }

    if (name[3]=='L' && name[4]=='5'){
      hist_L5->Fill(matchesAll);
      hist_L5pass->Fill(matchesPass);
    }

    if (name[3]=='L' && name[4]=='6'){
      hist_L6->Fill(matchesAll);
      hist_L6pass->Fill(matchesPass);
    }

    if (name[3]=='D' && name[4]=='1'){
      hist_D1->Fill(matchesAll);
      hist_D1pass->Fill(matchesPass);
    }

    if (name[3]=='D' && name[4]=='2'){
      hist_D2->Fill(matchesAll);
      hist_D2pass->Fill(matchesPass);
    }

    if (name[3]=='D' && name[4]=='3'){
      hist_D3->Fill(matchesAll);
      hist_D3pass->Fill(matchesPass);
    }

    if (name[3]=='D' && name[4]=='4'){
      hist_D4->Fill(matchesAll);
      hist_D4pass->Fill(matchesPass);
    }

    if (name[3]=='D' && name[4]=='5'){
      hist_D5->Fill(matchesAll);
      hist_D5pass->Fill(matchesPass);
    }

    std::map<TString, std::pair<TH1*,TH1*> >::iterator it=hists.find(name);
   
    if (it==hists.end()) {
      TString name1=name+"_All";
      TH1 *histAll = new TH1F(name1,name1,max,+0.5,max+0.5);
      histAll->Fill(matchesAll);
      name1=name+"_Pass";
      TH1 *histPass = new TH1F(name1,name1,max,+0.5,max+0.5);
      histPass->Fill(matchesPass);
      std::pair<TH1*,TH1*> histpair(histAll,histPass);
      hists[name]=histpair;
    } else {
      hists[name].first->Fill(matchesAll);
      hists[name].second->Fill(matchesPass);
   }

    count++;

  }

  cout << "count = "<<count<<endl;

  c1->cd(1);
  plotHist(hist1,0.05,ncut1,ncut2);
  hist2->SetLineColor(kBlue);
  hist2->Draw("same");

  c1->cd(2);
  plotHist(hist3,0.05,ncut1,ncut2);
  hist4->SetLineColor(kBlue);
  hist4->Draw("same");

  c1->Print("matchcalculator.pdf(","pdf");

  c2->cd(1);
  plotHist(hist_L1,0.05,ncut1,ncut2);
  hist_L1pass->SetLineColor(kBlue);
  hist_L1pass->Draw("same");

  c2->cd(2);
  plotHist(hist_L2,0.05,ncut1,ncut2);
  hist_L2pass->SetLineColor(kBlue);
  hist_L2pass->Draw("same");

  c2->cd(3);
  plotHist(hist_L3,0.05,ncut1,ncut2);
  hist_L3pass->SetLineColor(kBlue);
  hist_L3pass->Draw("same");

  c2->cd(4);
  plotHist(hist_L4,0.05,ncut1,ncut2);
  hist_L4pass->SetLineColor(kBlue);
  hist_L4pass->Draw("same");

  c2->cd(5);
  plotHist(hist_L5,0.05,ncut1,ncut2);
  hist_L5pass->SetLineColor(kBlue);
  hist_L5pass->Draw("same");

  c2->cd(6);
  plotHist(hist_L6,0.05,ncut1,ncut2);
  hist_L6pass->SetLineColor(kBlue);
  hist_L6pass->Draw("same");

  c2->Print("matchcalculator.pdf(","pdf");

  c3->cd(1);
  plotHist(hist_D1,0.05,ncut1,ncut2);
  hist_D1pass->SetLineColor(kBlue);
  hist_D1pass->Draw("same");

  c3->cd(2);
  plotHist(hist_D2,0.05,ncut1,ncut2);
  hist_D2pass->SetLineColor(kBlue);
  hist_D2pass->Draw("same");

  c3->cd(3);
  plotHist(hist_D3,0.05,ncut1,ncut2);
  hist_D3pass->SetLineColor(kBlue);
  hist_D3pass->Draw("same");

  c3->cd(4);
  plotHist(hist_D4,0.05,ncut1,ncut2);
  hist_D4pass->SetLineColor(kBlue);
  hist_D4pass->Draw("same");

  c3->cd(5);
  plotHist(hist_D5,0.05,ncut1,ncut2);
  hist_D5pass->SetLineColor(kBlue);
  hist_D5pass->Draw("same");

  c3->Print("matchcalculator.pdf(","pdf");

  int pages=0;

  std::map<TString, std::pair<TH1*,TH1*> >::iterator it=hists.begin();

  TCanvas* c=0;

  while(it!=hists.end()) {

    if (pages%4==0) {
     
      c = new TCanvas(it->first,"Track performance",200,50,600,700);
      c->Divide(2,2);
      c->SetFillColor(0);
      c->SetGrid();

    }

    c->cd(pages%4+1);
    //gPad->SetLogy();
    plotHist(it->second.first,0.05,ncut1,ncut2);
    it->second.second->SetLineColor(kBlue);
    it->second.second->Draw("same");
    
    pages++;

    if (pages%4==0) {
      c->Print("matchcalculator.pdf","pdf");
    }

    ++it;


 }

  c->Print("matchcalculator.pdf)","pdf");

}
