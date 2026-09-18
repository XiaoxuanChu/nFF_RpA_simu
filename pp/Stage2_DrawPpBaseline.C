#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

// Draw the merged pp baseline. Uses the histogram names written by
// Stage2_pp_Baseline.C / Stage2_MergePpBaseline.C.
void Stage2_DrawPpBaseline(const char* inputFile="Stage2_pp_baseline_full.root",
                           const char* outputDir="Stage2_pp_plots") {
  gSystem->mkdir(outputDir, kTRUE);
  gStyle->SetOptStat(0);
  TFile* input=TFile::Open(inputFile,"READ");
  if(!input || input->IsZombie()) { std::cout << "ERROR: cannot open " << inputFile << std::endl; return; }
  const char* species[3]={"pi","ka","pr"};
  const char* labels[3]={"#pi^{+}+#pi^{-}","K^{+}+K^{-}","p+#bar{p}"};
  const char* jetBins[2]={"jetpt4_6","jetpt6_9"};
  const char* jetLabels[2]={"4 #leq p_{T}^{jet} < 6 GeV/c","6 #leq p_{T}^{jet} < 9 GeV/c"};
  const char* variables[2]={"Z","JT"};
  const char* xTitles[2]={"z = |p_{h}|/|p_{jet}|","j_{T} (GeV/c)"};
  const char* yTitles[2]={"D_{h}(z)","D_{h}(j_{T}) ((GeV/c)^{-1})"};
  int a,b;
  for(a=0;a<2;++a) for(b=0;b<2;++b) {
    TCanvas* canvas=new TCanvas(Form("c_pp_%s_%s",variables[b],jetBins[a]),"",1500,500);
    canvas->Divide(3,1);
    int s;
    for(s=0;s<3;++s) {
      TH1D* raw=(TH1D*)input->Get(Form("hRaw%s_%s_%s",variables[b],species[s],jetBins[a]));
      TH1D* ue=(TH1D*)input->Get(Form("hUe%s_%s_%s",variables[b],species[s],jetBins[a]));
      TH1D* corr=(TH1D*)input->Get(Form("hCorr%s_%s_%s",variables[b],species[s],jetBins[a]));
      canvas->cd(s+1);
      if(!raw || !ue || !corr) { TLatex note; note.SetNDC(); note.DrawLatex(0.15,0.5,"Histogram missing"); continue; }
      raw->SetTitle(Form("%s; %s; %s",labels[s],xTitles[b],yTitles[b]));
      raw->SetLineColor(kBlack); raw->SetMarkerColor(kBlack); raw->SetMarkerStyle(20);
      ue->SetLineColor(kBlue+1); ue->SetMarkerColor(kBlue+1); ue->SetMarkerStyle(24);
      corr->SetLineColor(kRed+1); corr->SetMarkerColor(kRed+1); corr->SetMarkerStyle(21);
      raw->Draw("E1"); ue->Draw("E1 SAME"); corr->Draw("E1 SAME");
      TLatex title; title.SetNDC(); title.SetTextSize(0.040); title.DrawLatex(0.15,0.93,jetLabels[a]);
      TLegend* legend=new TLegend(0.48,0.67,0.88,0.87); legend->SetBorderSize(0); legend->SetFillStyle(0);
      legend->AddEntry(raw,"in jet","lep"); legend->AddEntry(ue,"off-axis UE","lep"); legend->AddEntry(corr,"UE corrected","lep"); legend->Draw();
    }
    canvas->SaveAs(Form("%s/Stage2_pp_%s_%s.png",outputDir,variables[b],jetBins[a]));
    canvas->SaveAs(Form("%s/Stage2_pp_%s_%s.pdf",outputDir,variables[b],jetBins[a]));
    delete canvas;
  }
  input->Close(); delete input;
  std::cout << "Wrote pp baseline plots to " << outputDir << std::endl;
}

