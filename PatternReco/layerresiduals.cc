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

class PlotResiduals {

public:
  PlotResiduals(int layer, int seedindex){
    layer_=layer;
    seedindex_=seedindex;

    double drphimax= 2.0; //0.2;
    double dzmax=15.0; // 10.0 ;
    if ((layer_>4)&&(seedindex_==1)) drphimax=1.0;
    if (seedindex_==1) dzmax=5.0; //10
    //if ((layer_>4)&&(seedindex_==1)) drphimax=0.1;

    string name="r*phi residual in layer "+std::to_string(layer)+", "+seedname(seedindex)
      +" seed, pt<3 GeV";
    hist16l_ = new TH1F(name.c_str(),name.c_str(),25,-drphimax,drphimax);
    name="r*iphi residual in layer "+std::to_string(layer)+", "+seedname(seedindex)
      +" seed, pt<3 GeV";
    hist116l_ = new TH1F(name.c_str(),name.c_str(),25,-drphimax,drphimax);
    name="r*phi residual in layer "+std::to_string(layer)+", "+seedname(seedindex)
      +" seed, 3<pt<8 GeV";
    hist16m_ = new TH1F(name.c_str(),name.c_str(),25,-drphimax,drphimax);
    name="r*iphi residual in layer "+std::to_string(layer)+", "+seedname(seedindex)
      +" seed, 3<pt<8 GeV";
    hist116m_ = new TH1F(name.c_str(),name.c_str(),25,-drphimax,drphimax);
    name="r*phi residual in layer "+std::to_string(layer)+", "+seedname(seedindex)
      +" seed, pt>8 GeV";
    hist16h_ = new TH1F(name.c_str(),name.c_str(),25,-drphimax,drphimax);
    name="r*iphi residual in layer "+std::to_string(layer)+", "+seedname(seedindex)
      +" seed, pt>8 GeV";
    hist116h_ = new TH1F(name.c_str(),name.c_str(),25,-drphimax,drphimax);
    
    name="z residual in layer "+std::to_string(layer)+", "+seedname(seedindex);
    hist16_ = new TH1F(name.c_str(),name.c_str(),10,-dzmax,dzmax);
    name="iz residual in layer "+std::to_string(layer)+", "+seedname(seedindex);
    hist116_ = new TH1F(name.c_str(),name.c_str(),10,-dzmax,dzmax);
}

  string seedname(int seedindex) {

    if (seedindex==0) return "L1L2";
    if (seedindex==1) return "L2L3";
    if (seedindex==2) return "L3L4";
    if (seedindex==3) return "L5L6";
    if (seedindex==4) return "D1D2";
    if (seedindex==5) return "D3D4";
    if (seedindex==6) return "L1D1";
    if (seedindex==7) return "L2D1";
    if (seedindex==8) return "L3L4L2";
    if (seedindex==9) return "L5L6L4";
    if (seedindex==10) return "L2L3D1";
    if (seedindex==11) return "D1D2L2"; 


    return "Unkown seedindex";
    
  }
  
  void addResid(int layer, int seedindex, double pt, double idphi, double dphi, double dphicut, double idz, double dz, double dzcut){

    bool lowpt=(pt<3.0);
    bool highpt=(pt>8.0);
    bool medpt=(!lowpt)&&(!highpt);
    
    if (layer==layer_&&seedindex==seedindex_) {
      dphicut_=dphicut;
      dzcut_=dzcut;
      if (lowpt) {
	hist116l_->Fill(idphi);
	hist16l_->Fill(dphi);
      }
      if (medpt) {
	hist116m_->Fill(idphi);
	hist16m_->Fill(dphi);
      }
      if (highpt) {
	hist116h_->Fill(idphi);
	hist16h_->Fill(dphi);
      }
      hist116_->Fill(idz);
      hist16_->Fill(dz);

    }


  }

  void Draw(TCanvas* c){

    cout << "layer seed phicut : "<<layer_<<" "<<seedindex_<<" "<<dphicut_<<endl;

    TLegend* leg3 = new TLegend(0.37,0.64,0.59,0.85);
    leg3->SetTextSize(0.08);
    leg3->SetFillColor(0);
    leg3->SetBorderSize(0);


    c->cd(1);
    hist16l_->SetLineColor(kBlue);

    hist116l_->GetXaxis()->SetNdivisions(6);
    hist116l_->SetLineColor(kRed);
    hist116l_->Draw();
    hist116l_->SetMaximum(2.61 *hist16l_->GetMaximum() );
    hist16l_->Draw("Same");


    TLine* ll1 = new TLine(-dphicut_,0,-dphicut_,0.5*hist116l_->GetMaximum());
    ll1->Draw();
    TLine* ll2 = new TLine(dphicut_,0,dphicut_,0.5*hist116l_->GetMaximum());
    ll2->Draw();

    leg3->AddEntry( hist16l_, "float", "L" );
    leg3->AddEntry( hist116l_, "int", "L");
    leg3->Draw();
    c->Update();


    c->cd(2);
    hist16m_->SetLineColor(kBlue);

    hist116m_->GetXaxis()->SetNdivisions(6);
    hist116m_->SetLineColor(kRed);
    hist116m_->Draw();
    hist116m_->SetMaximum(1.61 *hist116m_->GetMaximum() );
    hist16m_->Draw("Same");

    TLine* lm1 = new TLine(-dphicut_,0,-dphicut_,0.5*hist116m_->GetMaximum());
    lm1->Draw();
    TLine* lm2 = new TLine(dphicut_,0,dphicut_,0.5*hist116m_->GetMaximum());
    lm2->Draw();

    c->cd(3);
    hist16h_->SetLineColor(kBlue);


    hist116h_->GetXaxis()->SetNdivisions(6);
    hist116h_->SetLineColor(kRed);
    hist116h_->Draw();
    hist116h_->SetMaximum(1.61 *hist116h_->GetMaximum() );
    hist16h_->Draw("Same");

    TLine* lh1 = new TLine(-dphicut_,0,-dphicut_,0.5*hist116h_->GetMaximum());
    lh1->Draw();
    TLine* lh2 = new TLine(dphicut_,0,dphicut_,0.5*hist116h_->GetMaximum());
    lh2->Draw();

    c->cd(4);
    hist16_->SetLineColor(kBlue);
    hist116_->GetXaxis()->SetNdivisions(15);
    hist116_->SetLineColor(kRed);
    hist116_->Draw();
    hist116_->SetMaximum(1.61 *hist116_->GetMaximum() );
    hist16_->Draw("Same");
    //TLine* l1 = new TLine(-drcut_,0,-drcut_,0.5*hist116_->GetMaximum());
    //l1->Draw();
    //TLine* l2 = new TLine(drcut_,0,drcut_,0.5*hist116_->GetMaximum());
    //l2->Draw();

  }
  
private:

  int layer_;
  int seedindex_;
  double dphicut_;
  double dzcut_;
  
  TLegend *leg3;

  TH1 *hist16l_;
  TH1 *hist116l_;
  TH1 *hist16m_;
  TH1 *hist116m_;
  TH1 *hist16h_;
  TH1 *hist116h_;
  TH1 *hist16_;
  TH1 *hist116_;


};


void layerresiduals(){
//
// To see the output of this macro, click here.

//

#include "TMath.h"

gROOT->Reset();

gROOT->SetStyle("Plain");

gStyle->SetCanvasColor(kWhite);

gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
gStyle->SetPadBorderMode(0);
gStyle->SetOptStat(0);
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




 TCanvas* c1 = new TCanvas("c1","Track performance",200,10,700,800);
 c1->Divide(2,2);
 c1->SetFillColor(0);
 c1->SetGrid();

 PlotResiduals Resid_L3_L1L2(3,0);
 PlotResiduals Resid_L4_L1L2(4,0);
 PlotResiduals Resid_L5_L1L2(5,0);
 PlotResiduals Resid_L6_L1L2(6,0);

 PlotResiduals Resid_L4_L2L3(4,1);
 PlotResiduals Resid_L5_L2L3(5,1);
 PlotResiduals Resid_L6_L2L3(6,1);

 PlotResiduals Resid_L1_L3L4(1,2);
 PlotResiduals Resid_L2_L3L4(2,2);
 PlotResiduals Resid_L5_L3L4(5,2);
 PlotResiduals Resid_L6_L3L4(6,2);

 PlotResiduals Resid_L1_L5L6(1,3);
 PlotResiduals Resid_L2_L5L6(2,3);
 PlotResiduals Resid_L3_L5L6(3,3);
 PlotResiduals Resid_L4_L5L6(4,3);

 PlotResiduals Resid_L1_L2L3(1,1);
 PlotResiduals Resid_L1_D1D2(1,4);
 PlotResiduals Resid_L1_D3D4(1,5);
 PlotResiduals Resid_L1_L2D1(1,7);
 
 PlotResiduals Resid_L1_L3L4L2(1,8);
 PlotResiduals Resid_L1_L5L6L4(1,9);
 PlotResiduals Resid_L1_L2L3D1(1,10); 
 PlotResiduals Resid_L1_D1D2L2(1,11);

 PlotResiduals Resid_L2_D1D2(2,4);
 
 PlotResiduals Resid_L2_L5L6L4(2,9);

 PlotResiduals Resid_L3_L5L6L4(3,9);
 //PlotResiduals Resid_L3_D1D2L2(3,11);
 PlotResiduals Resid_L4_L2L3D1(4,10); 
 //PlotResiduals Resid_L4_D1D2L2(4,11);

 PlotResiduals Resid_L6_L3L4L2(6,8);

 ifstream in("layerresiduals.txt");

 int count=0;

 while (in.good()) {

   double layer,seedindex,pt,idphi,dphi,dphicut,idz,dz,dzcut, rinv;
   
   in>>layer>>seedindex>>pt>>idphi>>dphi>>dphicut>>idz>>dz>>dzcut>>rinv;

   if (!in.good()) continue;

   Resid_L3_L1L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L4_L1L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L5_L1L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L6_L1L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);

   Resid_L4_L2L3.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L5_L2L3.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L6_L2L3.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);

   Resid_L1_L3L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L2_L3L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L5_L3L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L6_L3L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);

   Resid_L1_L5L6.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L2_L5L6.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L3_L5L6.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L4_L5L6.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   
   Resid_L1_L2L3.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L1_D1D2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L1_D3D4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L1_L2D1.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
 
   Resid_L1_L3L4L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut); // added = done
   Resid_L1_L5L6L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L1_L2L3D1.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut); 
   Resid_L1_D1D2L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);

   Resid_L2_D1D2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
 
   Resid_L2_L5L6L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
 
   Resid_L3_L5L6L4.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   //Resid_L3_D1D2L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L4_L2L3D1.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);
   Resid_L6_L3L4L2.addResid(layer, seedindex, pt, idphi, dphi, dphicut, idz, dz, dzcut);

   count++;

 }

//cout << "Processed: "<<count<<" events"<<endl;
 Resid_L1_L2L3.Draw(c1);
 c1->Print("layerresiduals.pdf(","pdf");
 Resid_L1_L2L3D1.Draw(c1);
 c1->Print("layerresiduals.pdf(","pdf"); 
 Resid_L1_L3L4.Draw(c1);
 c1->Print("layerresiduals.pdf(","pdf");
 Resid_L1_L3L4L2.Draw(c1);
 c1->Print("layerresiduals.pdf(","pdf");
 Resid_L1_L5L6.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L1_L5L6L4.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L1_D1D2.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L1_D1D2L2.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf"); 
 Resid_L1_D3D4.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L1_L2D1.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L2_L3L4.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L2_L5L6.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L2_L5L6L4.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L2_D1D2.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L3_L1L2.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L3_L5L6.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L3_L5L6L4.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L4_L1L2.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L4_L5L6.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L4_L2L3.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L4_L2L3D1.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L5_L1L2.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L5_L3L4.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L5_L2L3.Draw(c1);
 c1->Print("layerresiduals.pdf","pdf");
 Resid_L6_L1L2.Draw(c1);
 c1->Print("layerresiduals_ALLprompt.pdf","pdf");
 Resid_L6_L3L4.Draw(c1);
 c1->Print("layerresiduals.pdf)","pdf");
 //Resid_L6_L3L4L2.Draw(c1);
 //c1->Print("layerresiduals.pdf)","pdf");
 //Resid_L6_L2L3.Draw(c1);
 //c1->Print("layerresiduals.pdf","pdf");





}

