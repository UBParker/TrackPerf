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
  PlotResiduals(int disk, int isPS,  int seedindex){
    disk_=disk;
    isPS_=isPS;
    seedindex_=seedindex;

    double drphimax=0.6;
    double drmax=5.0;
    int phibins = 15 ; //50;
    int rbins = 20 ; //20;
    if (seedindex>7) { //displaced seeds
      drmax=10.0;
      drphimax=1.5;
      phibins = 15;
      rbins = 20;
    };
    if (isPS==1) drmax=3.0;
    if (disk ==1 && isPS ==0 && seedindex==2){
    drphimax=1.5; 
    };

    string name="r*phi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt<3 GeV";
    hist16l_ = new TH1F(name.c_str(),name.c_str(),phibins,-drphimax,drphimax);
    name="r*iphi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt<3 GeV";
    hist116l_ = new TH1F(name.c_str(),name.c_str(),phibins,-drphimax,drphimax);
    name="r*phi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, 3<pt<8 GeV";
    hist16m_ = new TH1F(name.c_str(),name.c_str(),phibins,-drphimax,drphimax);
    name="r*iphi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, 3<pt<8 GeV";
    hist116m_ = new TH1F(name.c_str(),name.c_str(),phibins,-drphimax,drphimax);
    name="r*phi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt>8 GeV";
    hist16h_ = new TH1F(name.c_str(),name.c_str(),phibins,-drphimax,drphimax);
    name="r*iphi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt>8 GeV";
    hist116h_ = new TH1F(name.c_str(),name.c_str(),phibins,-drphimax,drphimax);
    
    name="r residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex);
    hist16_ = new TH1F(name.c_str(),name.c_str(),rbins,-drmax,drmax);
    name="ir residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex);
    hist116_ = new TH1F(name.c_str(),name.c_str(),rbins,-drmax,drmax);


    // for + charge tracks

    string namep="p r*phi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt<3 GeV";
    hist16lp_ = new TH1F(namep.c_str(),namep.c_str(),phibins,-drphimax,drphimax);
    namep="p r*iphi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt<3 GeV";
    hist116lp_ = new TH1F(namep.c_str(),namep.c_str(),phibins,-drphimax,drphimax);
    namep="p r*phi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, 3<pt<8 GeV";
    hist16mp_ = new TH1F(namep.c_str(),namep.c_str(),phibins,-drphimax,drphimax);
    namep=" p r*iphi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, 3<pt<8 GeV";
    hist116mp_ = new TH1F(namep.c_str(),namep.c_str(),phibins,-drphimax,drphimax);
    namep="p r*phi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt>8 GeV";
    hist16hp_ = new TH1F(namep.c_str(),namep.c_str(),phibins,-drphimax,drphimax);
    namep="p r*iphi residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex)+" seed, pt>8 GeV";
    hist116hp_ = new TH1F(namep.c_str(),namep.c_str(),phibins,-drphimax,drphimax);
    
    namep="p r residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex);
    hist16p_ = new TH1F(namep.c_str(),namep.c_str(),rbins,-drmax,drmax);
    namep="p ir residual in D"+std::to_string(disk)+ps2s(isPS)+", "+seedname(seedindex);
    hist116p_ = new TH1F(namep.c_str(),namep.c_str(),rbins,-drmax,drmax);    



}

  int addResid(int disk, int isPS, int seedindex, double pt, double idphi, double dphi, double dphicut, double idr, double dr, double drcut, double rinv){

    bool lowpt=(pt<3.0);
    bool highpt=(pt>8.0);
    bool medpt=(!lowpt)&&(!highpt);
    
    if (disk==disk_&&seedindex==seedindex_&&isPS==isPS_) {
      dphicut_=dphicut;
      drcut_=drcut;

      if (rinv < 0) {
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
      hist116_->Fill(idr);
      hist16_->Fill(dr);
    }
      if (rinv > 0) {
      if (lowpt) {
  hist116lp_->Fill(idphi);
  hist16lp_->Fill(dphi);
      }
      if (medpt) {
  hist116mp_->Fill(idphi);
  hist16mp_->Fill(dphi);
      }
      if (highpt) {
  hist116hp_->Fill(idphi);
  hist16hp_->Fill(dphi);
      }
      hist116p_->Fill(idr);
      hist16p_->Fill(dr);
    }
      return 1;
    }

    return 0;

  }

  void Draw(TCanvas* c){

    cout << "disk seedindex phicut : "<<disk_<<" "<<seedindex_<<" "<<dphicut_<<endl;


    TLegend* leg3 = new TLegend(0.37,0.64,0.59,0.85);
    leg3->SetTextSize(0.08);
    leg3->SetFillColor(0);
    leg3->SetBorderSize(0);



    hist16m_->SetLineColor(kBlue);
    hist16m_->GetXaxis()->SetNdivisions(6);
    hist16m_->Draw("hist");
    hist16m_->SetMaximum(2.21 *hist16m_->GetMaximum() );


    hist16mp_->SetLineColor(kRed);

    hist16mp_->Draw("hist Same");


    TLine* lm1 = new TLine(-dphicut_,0,-dphicut_,1.7*hist16m_->GetMaximum());
    lm1->Draw();
    TLine* lm2 = new TLine(dphicut_,0,dphicut_,1.7*hist16m_->GetMaximum());
    lm2->Draw();


    leg3->AddEntry( hist16m_, "- charge ", "L" );
    leg3->AddEntry( hist16mp_, "+ charge ", "L" );
    leg3->Draw();
    c->Update();
    leg3->Draw();
    c->Update();

   /*


    c->cd(2);
=======
    c->cd(3);
    hist16h_->SetLineColor(kBlue);
    hist16h_->Draw();
    hist16h_->SetMaximum(1.61 *hist16h_->GetMaximum() );
    hist116h_->SetLineColor(kRed);
    hist116h_->Draw("Same");
    TLine* lh1 = new TLine(-dphicut_,0,-dphicut_,0.5*hist116h_->GetMaximum());
    lh1->Draw();
    TLine* lh2 = new TLine(dphicut_,0,dphicut_,0.5*hist116h_->GetMaximum());
    lh2->Draw();

    c->cd(4);
>>>>>>> 99300ab49fe70d09593989b91f012cc2cd2c0383
    hist16_->SetLineColor(kBlue);

    hist16_->GetXaxis()->SetNdivisions(10);
    hist16_->Draw();
<<<<<<< HEAD
    hist16_->SetMaximum(3.61 *hist16_->GetMaximum() );
    hist116p_->SetMarkerColor(kRed);
    hist116p_->SetMarkerStyle(20);
    hist116p_->SetMarkerSize(0.5);

    hist116p_->Draw("l Same");

    hist16p_->SetMarkerColor(kBlue);
    hist16p_->SetMarkerStyle(20);
    hist16p_->SetMarkerSize(0.5);

    hist16p_->Draw("l Same");
    //hist16l_->SetMaximum(3.61 *hist16l_->GetMaximum() );
=======
    hist16_->SetMaximum(1.61 *hist16_->GetMaximum() );
>>>>>>> 99300ab49fe70d09593989b91f012cc2cd2c0383
    hist116_->SetLineColor(kRed);
    //
    hist116_->Draw("l Same");

    TLine* l1 = new TLine(-drcut_,0,-drcut_,0.5*hist116p_->GetMaximum());
    l1->Draw();
    TLine* l2 = new TLine(drcut_,0,drcut_,0.5*hist116p_->GetMaximum());
    l2->Draw();

    leg3->AddEntry( hist16l_, "- charge float", "L" );
    leg3->AddEntry( hist116lp_, "+ charge int", "p");
    leg3->AddEntry( hist16lp_, "+ charge float", "p" );
    leg3->AddEntry( hist116l_, "- charge int", "L");
    leg3->Draw();
    c->Update();
    leg3->Draw();
    c->Update();

    */


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

  string ps2s(int isPS){

    if (isPS) return "PS";

    return "2S";

  }
  
private:

  int disk_;
  int isPS_;
  int seedindex_;
  double dphicut_;
  double drcut_;
  
  TLegend *leg3;  

  TH1 *hist16l_;
  TH1 *hist116l_;
  TH1 *hist16m_;
  TH1 *hist116m_;
  TH1 *hist16h_;
  TH1 *hist116h_;
  TH1 *hist16_;
  TH1 *hist116_;

  TH1 *hist16lp_;
  TH1 *hist116lp_;
  TH1 *hist16mp_;
  TH1 *hist116mp_;
  TH1 *hist16hp_;
  TH1 *hist116hp_;
  TH1 *hist16p_;
  TH1 *hist116p_;

};


void diskresiduals(){
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
  gStyle->SetLineWidth(2);
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
 //c1->Divide(1);
 c1->SetFillColor(0);
 //c1->SetGrid();

 PlotResiduals Resid_D1PS_L1L2(1,1,0);
 PlotResiduals Resid_D12S_L1L2(1,0,0);
 PlotResiduals Resid_D1PS_L2L3(1,1,1);
 PlotResiduals Resid_D12S_L2L3(1,0,1);

 PlotResiduals Resid_D12S_L3L4(1,0,2);
 PlotResiduals Resid_D1PS_D3D4(1,1,5);
		
 PlotResiduals Resid_D2PS_L1L2(2,1,0);
 PlotResiduals Resid_D22S_L1L2(2,0,0);
 PlotResiduals Resid_D2PS_L2L3(2,1,1);
 PlotResiduals Resid_D22S_L2L3(2,0,1); 
 PlotResiduals Resid_D22S_L3L4(2,0,2);
 PlotResiduals Resid_D2PS_D3D4(2,1,5);
 PlotResiduals Resid_D2PS_L1D1(2,1,6);
 PlotResiduals Resid_D2PS_L2D1(2,1,7);
 PlotResiduals Resid_D22S_L2D1(2,0,7);
		
 PlotResiduals Resid_D3PS_L1L2(3,1,0);
 PlotResiduals Resid_D32S_L1L2(3,0,0);
 PlotResiduals Resid_D3PS_D1D2(3,1,4);
 PlotResiduals Resid_D32S_D1D2(3,0,4);
 PlotResiduals Resid_D3PS_L1D1(3,1,6);
 PlotResiduals Resid_D3PS_L2D1(3,1,7);
 PlotResiduals Resid_D32S_L2D1(3,0,7);
		
 PlotResiduals Resid_D42S_L1L2(4,0,0);
 PlotResiduals Resid_D4PS_D1D2(4,1,4);
 PlotResiduals Resid_D42S_D1D2(4,0,4);
 PlotResiduals Resid_D4PS_L1D1(4,1,6);
 PlotResiduals Resid_D42S_L1D1(4,0,6);
 PlotResiduals Resid_D42S_L2D1(4,0,7);
		
 PlotResiduals Resid_D5PS_D1D2(5,1,4);
 PlotResiduals Resid_D52S_D1D2(5,0,4);
 PlotResiduals Resid_D5PS_D3D4(5,1,5);
 PlotResiduals Resid_D52S_D3D4(5,0,5);
 PlotResiduals Resid_D5PS_L1D1(5,1,6);
 PlotResiduals Resid_D52S_L1D1(5,0,6);
 
 PlotResiduals Resid_D5PS_L2D1(5,1,7);
 PlotResiduals Resid_D52S_L2D1(5,0,7);

 PlotResiduals Resid_D1PS_L3L4L2(1,1,8);

 PlotResiduals Resid_D12S_L3L4L2(1,0,8);

 PlotResiduals Resid_D2PS_L3L4L2(2,1,8);
 PlotResiduals Resid_D22S_L3L4L2(2,0,8);


 PlotResiduals Resid_D3PS_L3L4L2(3,1,8);
 PlotResiduals Resid_D4PS_L3L4L2(4,1,8); 
 PlotResiduals Resid_D5PS_L3L4L2(5,1,8);
 PlotResiduals Resid_D52S_L3L4L2(5,0,8);

 PlotResiduals Resid_D1PS_L5L6L4(1,1,9);
 PlotResiduals Resid_D12S_L5L6L4(1,0,9);

 PlotResiduals Resid_D5PS_L2L3D1(5,1,10);
 PlotResiduals Resid_D52S_L2L3D1(5,0,10);
 PlotResiduals Resid_D4PS_L2L3D1(4,1,10);
 PlotResiduals Resid_D42S_L2L3D1(4,0,10);
 PlotResiduals Resid_D3PS_L2L3D1(3,1,10);
 PlotResiduals Resid_D32S_L2L3D1(3,0,10);
 PlotResiduals Resid_D2PS_L2L3D1(2,1,10);
 PlotResiduals Resid_D22S_L2L3D1(2,0,10);

 PlotResiduals Resid_D5PS_D1D2L2(5,1,11);
 PlotResiduals Resid_D52S_D1D2L2(5,0,11);
   
 PlotResiduals Resid_D4PS_D1D2L2(4,1,11);
 PlotResiduals Resid_D42S_D1D2L2(4,0,11);


 PlotResiduals Resid_D3PS_D1D2L2(3,1,11);
 PlotResiduals Resid_D32S_D1D2L2(3,0,11);

 ifstream in("diskresiduals_false_10000events.txt");
 //"diskresiduals_disp_10000evptl5.txt");


 int count=0;

 while (in.good()) {

   double disk,isPS,pt,idphi,dphi,dphicut,idr,dr,drcut,rinv;
   int seedindex=-1;

   in>>disk>>isPS>>seedindex>>pt>>idphi>>dphi>>dphicut>>idr>>dr>>drcut>>rinv;

   if (!in.good()) continue;

   int added=0;
   

   added+=Resid_D1PS_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D12S_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   
   added+=Resid_D1PS_L2L3.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D12S_L2L3.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D12S_L3L4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D1PS_D3D4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
		
   added+=Resid_D2PS_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D22S_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D2PS_L2L3.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D22S_L2L3.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
      
   added+=Resid_D22S_L3L4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D2PS_D3D4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D2PS_L1D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D2PS_L2D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D22S_L2D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
		
   added+=Resid_D3PS_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D3PS_D1D2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_D1D2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D3PS_L1D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D3PS_L2D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_L2D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
		
   added+=Resid_D42S_L1L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D4PS_D1D2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_D1D2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D4PS_L1D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_L1D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_L2D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );   
		
   added+=Resid_D5PS_D1D2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_D1D2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D5PS_D3D4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_D3D4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D5PS_L1D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_L1D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );


   added+=Resid_D1PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D12S_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   
   added+=Resid_D2PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D22S_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D3PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D4PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D5PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D52S_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D1PS_L5L6L4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D12S_L5L6L4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D5PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D4PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D3PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D2PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D22S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );



   added+=Resid_D5PS_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   
   added+=Resid_D3PS_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D4PS_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );


   added+=Resid_D1PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D2PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D3PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D4PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D5PS_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D52S_L3L4L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D1PS_L5L6L4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D12S_L5L6L4.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D5PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D4PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D3PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D2PS_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D22S_L2L3D1.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );



   added+=Resid_D5PS_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D52S_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   
   added+=Resid_D3PS_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D32S_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   added+=Resid_D4PS_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );
   added+=Resid_D42S_D1D2L2.addResid(disk, isPS, seedindex, pt, idphi, dphi, dphicut, idr, dr, drcut, rinv );

   if (added!=1) {
     cout << "Added = "<<added<<" : disk isPS seedindex "<<disk<<" "<<isPS<<" "<<seedindex<<endl;
   }
   
   count++;

 }

//cout << "Processed: "<<count<<" events"<<endl;

 // this one is empty
 //Resid_D1PS_L3L4L2.Draw(c1);
 //c1->Print("diskresiduals.pdf(","pdf");



 Resid_D12S_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf"); 
 

 Resid_D12S_L3L4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D22S_L3L4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");




 Resid_D2PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D22S_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D2PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf"); 
 Resid_D2PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D22S_L2L3.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D22S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D3PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");


 Resid_D4PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D3PS_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D3PS_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D4PS_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D4PS_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");








 Resid_D52S_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf)","pdf");



}

/*










 Resid_D1PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D12S_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf"); 
 Resid_D2PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D22S_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");



=======
 Resid_D1PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
>>>>>>> 99300ab49fe70d09593989b91f012cc2cd2c0383
 Resid_D12S_L5L6L4.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D1PS_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D12S_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D1PS_L2L3.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D12S_L2L3.Draw(c1);
 c1->Print("diskresiduals.pdf(","pdf");
 Resid_D12S_L3L4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D1PS_D3D4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 
 Resid_D2PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf"); 
 Resid_D2PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
<<<<<<< HEAD
 Resid_D22S_L2L3.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
=======
>>>>>>> 99300ab49fe70d09593989b91f012cc2cd2c0383
 Resid_D22S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D2PS_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D22S_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D22S_L3L4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D2PS_D3D4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D2PS_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D2PS_L2D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D22S_L2D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 
 Resid_D3PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D3PS_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D3PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D3PS_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D3PS_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D3PS_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D3PS_L2D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D32S_L2D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 
 Resid_D4PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D4PS_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_D1D2L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D4PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");

 Resid_D42S_L1L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D4PS_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D4PS_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D42S_L2D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf"); 
 
 Resid_D5PS_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D52S_L3L4L2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D5PS_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D52S_L2L3D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");


 Resid_D5PS_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D52S_D1D2.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D5PS_D3D4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D52S_D3D4.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D5PS_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf","pdf");
 Resid_D52S_L1D1.Draw(c1);
 c1->Print("diskresiduals.pdf)","pdf");
*/




