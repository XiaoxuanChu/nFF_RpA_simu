#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include <fstream>
#include <string>
#include <iostream>

// ROOT 5-safe merger. A partial histogram is a per-weighted-jet density;
// it is multiplied by the partial weighted-jet count before summation.
void Stage2_MergePpBaseline(const char* listFile="partial_roots.list",
                            const char* outName="Stage2_pp_baseline_full.root") {
  const char* species[3] = {"pi", "ka", "pr"};
  const char* jetBins[2] = {"jetpt4_6", "jetpt6_9"};
  const char* kinds[6] = {"hRawZ_", "hUeZ_", "hCorrZ_", "hRawJT_", "hUeJT_", "hCorrJT_"};
  TH1D* sumWeight = 0; TH1D* sumCount = 0; TH1D* sumHist[2][3][6];
  int a, b, c;
  for(a=0;a<2;++a) for(b=0;b<3;++b) for(c=0;c<6;++c) sumHist[a][b][c]=0;

  std::ifstream input(listFile);
  if(!input.good()) { std::cout << "ERROR: cannot open " << listFile << std::endl; return; }
  std::string fileName; int fileNumber=0, mergedFiles=0;
  while(input >> fileName) {
    std::cout << "Opening partial " << fileNumber << ": " << fileName << std::endl;
    ++fileNumber;
    TFile* inFile=TFile::Open(fileName.c_str(),"READ");
    if(!inFile || inFile->IsZombie()) { std::cout << "WARNING: cannot open; skipping" << std::endl; if(inFile){inFile->Close();delete inFile;} continue; }
    TH1D* inWeight=(TH1D*)inFile->Get("hJetWeight");
    TH1D* inCount=(TH1D*)inFile->Get("hJetCount");
    if(!inWeight || !inCount) { std::cout << "WARNING: missing jet denominator; skipping" << std::endl; inFile->Close();delete inFile;continue; }

    // Detach every source object while its source TFile is still open.
    TH1D* localWeight=(TH1D*)inWeight->Clone("Stage2_localWeight"); localWeight->SetDirectory(0);
    TH1D* localCount=(TH1D*)inCount->Clone("Stage2_localCount"); localCount->SetDirectory(0);
    TH1D* localHist[2][3][6]; bool complete=true;
    for(a=0;a<2;++a) for(b=0;b<3;++b) for(c=0;c<6;++c) {
      TString histName=TString(kinds[c])+species[b]+"_"+jetBins[a];
      TH1D* sourceHist=(TH1D*)inFile->Get(histName.Data());
      localHist[a][b][c]=0;
      if(!sourceHist) { std::cout << "WARNING: missing " << histName << "; skipping file" << std::endl; complete=false; }
      else { localHist[a][b][c]=(TH1D*)sourceHist->Clone("Stage2_localHist"); localHist[a][b][c]->SetDirectory(0); }
    }
    inFile->Close(); delete inFile;
    if(!complete) {
      delete localWeight; delete localCount;
      for(a=0;a<2;++a) for(b=0;b<3;++b) for(c=0;c<6;++c) if(localHist[a][b][c])delete localHist[a][b][c];
      continue;
    }
    if(!sumWeight) {
      sumWeight=(TH1D*)localWeight->Clone("hJetWeight"); sumWeight->SetDirectory(0); sumWeight->Reset();
      sumCount=(TH1D*)localCount->Clone("hJetCount"); sumCount->SetDirectory(0); sumCount->Reset();
    }
    sumWeight->Add(localWeight); sumCount->Add(localCount);
    for(a=0;a<2;++a) {
      Double_t jetWeight=localWeight->GetBinContent(a+1);
      for(b=0;b<3;++b) for(c=0;c<6;++c) {
        TString histName=TString(kinds[c])+species[b]+"_"+jetBins[a];
        if(!sumHist[a][b][c]) { sumHist[a][b][c]=(TH1D*)localHist[a][b][c]->Clone(histName.Data()); sumHist[a][b][c]->SetDirectory(0); sumHist[a][b][c]->Reset(); }
        sumHist[a][b][c]->Add(localHist[a][b][c],jetWeight);
      }
    }
    delete localWeight; delete localCount;
    for(a=0;a<2;++a) for(b=0;b<3;++b) for(c=0;c<6;++c) delete localHist[a][b][c];
    ++mergedFiles;
  }
  input.close();
  if(!sumWeight) { std::cout << "ERROR: no valid partial ROOT files" << std::endl; return; }
  for(a=0;a<2;++a) {
    Double_t jetWeight=sumWeight->GetBinContent(a+1);
    if(jetWeight<=0.0) { std::cout << "WARNING: non-positive weight for " << jetBins[a] << std::endl; continue; }
    for(b=0;b<3;++b) for(c=0;c<6;++c) sumHist[a][b][c]->Scale(1.0/jetWeight);
  }
  TFile* output=TFile::Open(outName,"RECREATE");
  if(!output || output->IsZombie()) { std::cout << "ERROR: cannot create " << outName << std::endl; return; }
  sumWeight->Write(); sumCount->Write();
  for(a=0;a<2;++a) for(b=0;b<3;++b) for(c=0;c<6;++c) sumHist[a][b][c]->Write();
  output->Close(); delete output;
  std::cout << "Merged " << mergedFiles << " partial ROOT files into " << outName << std::endl;
}

