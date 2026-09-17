// Stage 2: pp identified-hadron-in-jet baseline from existing ting pp trees.
//
// Input manifest columns:
//   pthat_min pthat_max jet_file ue_file
//
// The macro is read-only with respect to all input ROOT files.  It writes only
// outputRoot and outputText, which must point to the caller's own directory.
// Run with ROOT ACLiC in the legacy STAR environment, e.g.
// root4star -b -q 'Stage2_pp_Baseline.C+("Stage2_pp_inputs_absolute.manifest",
//     "Stage2_pp_baseline.root","Stage2_pp_baseline_summary.txt")'

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TMethodCall.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>

namespace Stage2PpBaseline {

const Int_t kNJetBins = 2;
const Int_t kNSpecies = 3;
const char* kSpecies[kNSpecies] = {"pi", "ka", "pr"};
const char* kJetBinName[kNJetBins] = {"jetpt4_6", "jetpt6_9"};
const Double_t kJetPtLow[kNJetBins] = {4.0, 6.0};
const Double_t kJetPtHigh[kNJetBins] = {6.0, 9.0};

const Int_t kNZ[kNJetBins] = {7, 6};
const Double_t kZEdges4_6[8] = {0.0, 0.05, 0.15, 0.30,
                                0.50, 0.70, 0.90, 1.00};
const Double_t kZEdges6_9[7] = {0.0, 0.05, 0.15, 0.30,
                                0.50, 0.70, 1.00};
const Int_t kNJT = 6;
const Double_t kJTEdges[7] = {0.0, 0.10, 0.20, 0.40,
                              0.60, 0.90, 1.00};

// Values transcribed from the PI-supplied pp_jet_lumi record, audited in Stage 1.
const Int_t kHardLow[11] = {2,3,4,5,7,9,11,15,20,25,35};
const Int_t kHardHigh[11] = {3,4,5,7,9,11,15,20,25,35,-1};
const Double_t kHardSigmaMb[11] = {9.0060,1.4620,0.3544,0.1514,
  0.02489,0.005846,0.002305,0.000343,0.00004560,0.000009740,0.000000502};
const Long64_t kHardNGenerated[11] = {3266550,928437,902278,602974,
  371968,348405,571960,339730,188186,151720,30048};

Bool_t Init(TMethodCall& call, const char* className,
            const char* methodName, const char* prototype)
{
  TClass* cl = TClass::GetClass(className);
  if (!cl) {
    std::cout << "ERROR: class unavailable: " << className << std::endl;
    return kFALSE;
  }
  call.InitWithPrototype(cl, methodName, prototype);
  if (!call.IsValid()) {
    std::cout << "ERROR: method unavailable: " << className << "::"
              << methodName << "(" << prototype << ")" << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

Long_t CallLong0(TMethodCall& call, void* object)
{
  Long_t result = 0;
  call.Execute(object, result);
  return result;
}

void* CallObject1(TMethodCall& call, void* object, Long_t index)
{
  Long_t result = 0;
  call.ResetParam();
  call.SetParam((Long_t)index);
  call.Execute(object, result);
  return (void*)result;
}

Double_t CallDouble0(TMethodCall& call, void* object)
{
  Double_t result = 0.0;
  call.Execute(object, result);
  return result;
}

Int_t SpeciesIndex(Long_t pdg)
{
  if (pdg < 0) pdg = -pdg;
  if (pdg == 211)  return 0; // pi+ or pi-
  if (pdg == 321)  return 1; // K+ or K-
  if (pdg == 2212) return 2; // p or pbar
  return -1;
}

Int_t GetHardBinIndex(Int_t low, Int_t high)
{
  for (Int_t i = 0; i < 11; ++i) {
    if (kHardLow[i] == low && kHardHigh[i] == high) return i;
  }
  return -1;
}

void FillParticle(TH1D* hZ, TH1D* hJT, Double_t weight,
                  Double_t particlePt, Double_t particleEta,
                  Double_t particlePhi, Double_t jetPt,
                  Double_t jetEta, Double_t jetPhi)
{
  const Double_t ppx = particlePt * cos(particlePhi);
  const Double_t ppy = particlePt * sin(particlePhi);
  const Double_t ppz = particlePt * sinh(particleEta);
  const Double_t pjx = jetPt * cos(jetPhi);
  const Double_t pjy = jetPt * sin(jetPhi);
  const Double_t pjz = jetPt * sinh(jetEta);
  const Double_t pMag2 = ppx*ppx + ppy*ppy + ppz*ppz;
  const Double_t jetMag2 = pjx*pjx + pjy*pjy + pjz*pjz;

  if (pMag2 <= 0.0 || jetMag2 <= 0.0) return;

  const Double_t pMag = sqrt(pMag2);
  const Double_t jetMag = sqrt(jetMag2);
  const Double_t z = pMag / jetMag; // STAR nFF analysis-note definition
  const Double_t dot = ppx*pjx + ppy*pjy + ppz*pjz;
  Double_t cross2 = pMag2*jetMag2 - dot*dot;
  if (cross2 < 0.0) cross2 = 0.0; // floating-point protection
  const Double_t jT = sqrt(cross2) / jetMag;

  if (z < 0.05) return;
  if (jT >= (0.025 + 0.3295*z)*jetPt) return;

  hZ->Fill(z, weight);
  hJT->Fill(jT, weight);
}

} // namespace Stage2PpBaseline

void Stage2_pp_Baseline(
    const char* manifestFile = "Stage2_pp_inputs_absolute.manifest",
    const char* outputRoot = "Stage2_pp_baseline.root",
    const char* outputText = "Stage2_pp_baseline_summary.txt",
    Int_t maxFilePairs = 0)
{
  using namespace Stage2PpBaseline;

  gSystem->Load("libStJetSkimEvent");
  gSystem->Load("libStJets");
  gSystem->Load("libStJetEvent");
  gSystem->Load("libStUeEvent");

  TMethodCall nJetsCall, jetAtCall, jetPtCall, jetEtaCall, jetPhiCall;
  TMethodCall nJetParticlesCall, jetParticleAtCall;
  TMethodCall nVerticesCall, vertexAtCall, nUeJetsCall, ueJetAtCall;
  TMethodCall nConesCall, coneAtCall, nConeParticlesCall, coneParticleAtCall;
  TMethodCall particlePdgCall, particleStatusCall;
  TMethodCall particlePtCall, particleEtaCall, particlePhiCall;

  Bool_t ok = kTRUE;
  ok = ok && Init(nJetsCall, "StJetEvent", "numberOfJets", "");
  ok = ok && Init(jetAtCall, "StJetEvent", "jet", "Int_t");
  ok = ok && Init(jetPtCall, "StJetCandidate", "pt", "");
  ok = ok && Init(jetEtaCall, "StJetCandidate", "eta", "");
  ok = ok && Init(jetPhiCall, "StJetCandidate", "phi", "");
  ok = ok && Init(nJetParticlesCall, "StJetCandidate", "numberOfParticles", "");
  ok = ok && Init(jetParticleAtCall, "StJetCandidate", "particle", "Int_t");
  ok = ok && Init(nVerticesCall, "StUeOffAxisConesEvent", "numberOfVertices", "");
  ok = ok && Init(vertexAtCall, "StUeOffAxisConesEvent", "vertex", "Int_t");
  ok = ok && Init(nUeJetsCall, "StUeVertex", "numberOfUeJets", "");
  ok = ok && Init(ueJetAtCall, "StUeVertex", "ueJet", "Int_t");
  ok = ok && Init(nConesCall, "StUeOffAxisConesJet", "numberOfCones", "");
  ok = ok && Init(coneAtCall, "StUeOffAxisConesJet", "cone", "Int_t");
  ok = ok && Init(nConeParticlesCall, "StUeOffAxisCones", "numberOfParticles", "");
  ok = ok && Init(coneParticleAtCall, "StUeOffAxisCones", "particle", "Int_t");
  ok = ok && Init(particlePdgCall, "StJetParticle", "pdg", "");
  ok = ok && Init(particleStatusCall, "StJetParticle", "status", "");
  ok = ok && Init(particlePtCall, "StJetParticle", "pt", "");
  ok = ok && Init(particleEtaCall, "StJetParticle", "eta", "");
  ok = ok && Init(particlePhiCall, "StJetParticle", "phi", "");
  if (!ok) return;

  TH1D* hRawZ[kNJetBins][kNSpecies];
  TH1D* hUeZ[kNJetBins][kNSpecies];
  TH1D* hCorrZ[kNJetBins][kNSpecies];
  TH1D* hRawJT[kNJetBins][kNSpecies];
  TH1D* hUeJT[kNJetBins][kNSpecies];
  TH1D* hCorrJT[kNJetBins][kNSpecies];
  TH1D* hJetWeight = new TH1D("hJetWeight", "weighted accepted jets;jet p_{T} bin;sum w", 2, 0.5, 2.5);
  TH1D* hJetCount = new TH1D("hJetCount", "accepted jets;jet p_{T} bin;number of jets", 2, 0.5, 2.5);
  hJetWeight->GetXaxis()->SetBinLabel(1, "4-6");
  hJetWeight->GetXaxis()->SetBinLabel(2, "6-9");
  hJetCount->GetXaxis()->SetBinLabel(1, "4-6");
  hJetCount->GetXaxis()->SetBinLabel(2, "6-9");
  hJetWeight->Sumw2();
  hJetCount->Sumw2();

  Int_t iJetBin = 0;
  Int_t iSpecies = 0;
  for (iJetBin = 0; iJetBin < kNJetBins; ++iJetBin) {
    const Double_t* zEdges = (iJetBin == 0) ? kZEdges4_6 : kZEdges6_9;
    for (iSpecies = 0; iSpecies < kNSpecies; ++iSpecies) {
      TString base = TString(kSpecies[iSpecies]) + "_" + kJetBinName[iJetBin];
      hRawZ[iJetBin][iSpecies] = new TH1D("hRawZ_" + base, "raw z;z;weighted hadron count", kNZ[iJetBin], zEdges);
      hUeZ[iJetBin][iSpecies] = new TH1D("hUeZ_" + base, "UE z;z;weighted hadron count", kNZ[iJetBin], zEdges);
      hRawJT[iJetBin][iSpecies] = new TH1D("hRawJT_" + base, "raw j_{T};j_{T} (GeV/c);weighted hadron count", kNJT, kJTEdges);
      hUeJT[iJetBin][iSpecies] = new TH1D("hUeJT_" + base, "UE j_{T};j_{T} (GeV/c);weighted hadron count", kNJT, kJTEdges);
      hRawZ[iJetBin][iSpecies]->Sumw2();
      hUeZ[iJetBin][iSpecies]->Sumw2();
      hRawJT[iJetBin][iSpecies]->Sumw2();
      hUeJT[iJetBin][iSpecies]->Sumw2();
    }
  }

  Long64_t entriesByHardBin[11];
  Long64_t filesByHardBin[11];
  for (Int_t ih = 0; ih < 11; ++ih) {
    entriesByHardBin[ih] = 0;
    filesByHardBin[ih] = 0;
  }

  Long64_t filesRead = 0;
  Long64_t filesFailed = 0;
  Long64_t eventsRead = 0;
  Long64_t eventsNoUeVertex = 0;
  Long64_t eventsJetUeCountMismatch = 0;
  Long64_t selectedJetsNoCone = 0;
  Long64_t selectedJets = 0;
  Long64_t manifestRowsUsed = 0;

  std::ifstream manifest(manifestFile);
  if (!manifest) {
    std::cout << "ERROR: cannot open manifest " << manifestFile << std::endl;
    return;
  }

  std::string line;
  while (std::getline(manifest, line)) {
    if (line.size() == 0 || line[0] == '#') continue;

    std::istringstream row(line);
    Int_t pthatLow = 0;
    Int_t pthatHigh = 0;
    std::string jetFileName;
    std::string ueFileName;
    if (!(row >> pthatLow >> pthatHigh >> jetFileName >> ueFileName)) {
      std::cout << "WARNING: malformed manifest row skipped: " << line << std::endl;
      continue;
    }
    if (maxFilePairs > 0 && manifestRowsUsed >= maxFilePairs) break;
    ++manifestRowsUsed;

    const Int_t hardIndex = GetHardBinIndex(pthatLow, pthatHigh);
    if (hardIndex < 0) {
      std::cout << "WARNING: unknown pThat bin skipped: " << line << std::endl;
      continue;
    }
    const Double_t weight = kHardSigmaMb[hardIndex] / (Double_t)kHardNGenerated[hardIndex];

    TFile* jetFile = TFile::Open(jetFileName.c_str(), "READ");
    TFile* ueFile = TFile::Open(ueFileName.c_str(), "READ");
    if (!jetFile || jetFile->IsZombie() || !ueFile || ueFile->IsZombie()) {
      std::cout << "WARNING: cannot open pair; skipped:\n  "
                << jetFileName << "\n  " << ueFileName << std::endl;
      ++filesFailed;
      if (jetFile) jetFile->Close();
      if (ueFile) ueFile->Close();
      continue;
    }

    TTree* jetTree = (TTree*)jetFile->Get("jet");
    TTree* ueTree = (TTree*)ueFile->Get("ue");
    TBranch* jetBranch = jetTree ? jetTree->GetBranch("AntiKtR060Particle") : 0;
    TBranch* ueBranch = ueTree ? ueTree->GetBranch("AntiKtR060ParticleOffAxisConesR060") : 0;
    if (!jetBranch || !ueBranch || jetTree->GetEntries() != ueTree->GetEntries()) {
      std::cout << "WARNING: invalid pair/tree/entry mismatch; skipped:\n  "
                << jetFileName << std::endl;
      ++filesFailed;
      jetFile->Close();
      ueFile->Close();
      continue;
    }

    void* jetEvent = 0;
    void* ueEvent = 0;
    jetBranch->SetAddress(&jetEvent);
    ueBranch->SetAddress(&ueEvent);
    const Long64_t nEntries = jetTree->GetEntries();
    entriesByHardBin[hardIndex] += nEntries;
    filesByHardBin[hardIndex] += 1;
    ++filesRead;

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
      jetBranch->GetEntry(iEvent);
      ueBranch->GetEntry(iEvent);
      ++eventsRead;

      const Long_t nJets = CallLong0(nJetsCall, jetEvent);
      const Long_t nVertices = CallLong0(nVerticesCall, ueEvent);
      if (nVertices < 1) {
        ++eventsNoUeVertex;
        continue;
      }
      void* vertex = CallObject1(vertexAtCall, ueEvent, 0);
      const Long_t nUeJets = CallLong0(nUeJetsCall, vertex);
      if (nJets != nUeJets) ++eventsJetUeCountMismatch;

      Long_t nCompare = nJets;
      if (nUeJets < nCompare) nCompare = nUeJets;
      for (Long_t iJet = 0; iJet < nCompare; ++iJet) {
        void* jet = CallObject1(jetAtCall, jetEvent, iJet);
        void* ueJet = CallObject1(ueJetAtCall, vertex, iJet);
        const Double_t jetPt = CallDouble0(jetPtCall, jet);
        const Double_t jetEta = CallDouble0(jetEtaCall, jet);
        const Double_t jetPhi = CallDouble0(jetPhiCall, jet);

        Int_t selectedJetBin = -1;
        for (iJetBin = 0; iJetBin < kNJetBins; ++iJetBin) {
          if (jetPt >= kJetPtLow[iJetBin] && jetPt < kJetPtHigh[iJetBin] &&
              jetEta >= -0.9 && jetEta <= 0.9) {
            selectedJetBin = iJetBin;
            break;
          }
        }
        if (selectedJetBin < 0) continue;

        const Long_t nCones = CallLong0(nConesCall, ueJet);
        if (nCones <= 0) {
          ++selectedJetsNoCone;
          continue;
        }

        ++selectedJets;
        hJetWeight->Fill(selectedJetBin + 1, weight);
        hJetCount->Fill(selectedJetBin + 1, 1.0);

        const Long_t nJetParticles = CallLong0(nJetParticlesCall, jet);
        for (Long_t iParticle = 0; iParticle < nJetParticles; ++iParticle) {
          void* particle = CallObject1(jetParticleAtCall, jet, iParticle);
          if (CallLong0(particleStatusCall, particle) != 1) continue;
          const Int_t species = SpeciesIndex(CallLong0(particlePdgCall, particle));
          if (species < 0) continue;
          FillParticle(hRawZ[selectedJetBin][species], hRawJT[selectedJetBin][species],
                       weight, CallDouble0(particlePtCall, particle),
                       CallDouble0(particleEtaCall, particle),
                       CallDouble0(particlePhiCall, particle),
                       jetPt, jetEta, jetPhi);
        }

        const Double_t ueWeight = weight / (Double_t)nCones;
        for (Long_t iCone = 0; iCone < nCones; ++iCone) {
          void* cone = CallObject1(coneAtCall, ueJet, iCone);
          const Long_t nConeParticles = CallLong0(nConeParticlesCall, cone);
          for (Long_t iParticle = 0; iParticle < nConeParticles; ++iParticle) {
            void* particle = CallObject1(coneParticleAtCall, cone, iParticle);
            if (CallLong0(particleStatusCall, particle) != 1) continue;
            const Int_t species = SpeciesIndex(CallLong0(particlePdgCall, particle));
            if (species < 0) continue;
            // STAR nFF note: evaluate UE particle z and jT relative to original jet.
            FillParticle(hUeZ[selectedJetBin][species], hUeJT[selectedJetBin][species],
                         ueWeight, CallDouble0(particlePtCall, particle),
                         CallDouble0(particleEtaCall, particle),
                         CallDouble0(particlePhiCall, particle),
                         jetPt, jetEta, jetPhi);
          }
        }
      }
    }

    jetFile->Close();
    ueFile->Close();
    if (filesRead % 100 == 0)
      std::cout << "Processed " << filesRead << " file pairs" << std::endl;
  }
  manifest.close();

  for (iJetBin = 0; iJetBin < kNJetBins; ++iJetBin) {
    const Double_t weightedJets = hJetWeight->GetBinContent(iJetBin + 1);
    for (iSpecies = 0; iSpecies < kNSpecies; ++iSpecies) {
      TString base = TString(kSpecies[iSpecies]) + "_" + kJetBinName[iJetBin];
      hCorrZ[iJetBin][iSpecies] = (TH1D*)hRawZ[iJetBin][iSpecies]->Clone("hCorrZ_" + base);
      hCorrJT[iJetBin][iSpecies] = (TH1D*)hRawJT[iJetBin][iSpecies]->Clone("hCorrJT_" + base);
      hCorrZ[iJetBin][iSpecies]->Add(hUeZ[iJetBin][iSpecies], -1.0);
      hCorrJT[iJetBin][iSpecies]->Add(hUeJT[iJetBin][iSpecies], -1.0);
      if (weightedJets > 0.0) {
        hRawZ[iJetBin][iSpecies]->Scale(1.0 / weightedJets, "width");
        hUeZ[iJetBin][iSpecies]->Scale(1.0 / weightedJets, "width");
        hCorrZ[iJetBin][iSpecies]->Scale(1.0 / weightedJets, "width");
        hRawJT[iJetBin][iSpecies]->Scale(1.0 / weightedJets, "width");
        hUeJT[iJetBin][iSpecies]->Scale(1.0 / weightedJets, "width");
        hCorrJT[iJetBin][iSpecies]->Scale(1.0 / weightedJets, "width");
      }
      hRawZ[iJetBin][iSpecies]->SetTitle("pp raw D_{h}(z);z;D_{h}(z)");
      hUeZ[iJetBin][iSpecies]->SetTitle("pp UE D_{h}(z);z;D_{h}^{UE}(z)");
      hCorrZ[iJetBin][iSpecies]->SetTitle("pp UE-corrected D_{h}(z);z;D_{h}^{corr}(z)");
      hRawJT[iJetBin][iSpecies]->SetTitle("pp raw D_{h}(j_{T});j_{T} (GeV/c);D_{h}(j_{T})");
      hUeJT[iJetBin][iSpecies]->SetTitle("pp UE D_{h}(j_{T});j_{T} (GeV/c);D_{h}^{UE}(j_{T})");
      hCorrJT[iJetBin][iSpecies]->SetTitle("pp UE-corrected D_{h}(j_{T});j_{T} (GeV/c);D_{h}^{corr}(j_{T})");
    }
  }

  TFile* out = TFile::Open(outputRoot, "RECREATE");
  if (!out || out->IsZombie()) {
    std::cout << "ERROR: cannot create output ROOT file " << outputRoot << std::endl;
    return;
  }
  TNamed convention("Stage2_particle_level_convention",
    "Existing ting AntiKtR060Particle status-1 jets; neutrinos retained in jet clustering by PI decision. "
    "Species: pi+/pi-, K+/K-, p/pbar. UE is mean of off-axis cones; UE z and jT use original jet momentum.");
  convention.Write();
  hJetWeight->Write();
  hJetCount->Write();
  for (iJetBin = 0; iJetBin < kNJetBins; ++iJetBin) {
    for (iSpecies = 0; iSpecies < kNSpecies; ++iSpecies) {
      hRawZ[iJetBin][iSpecies]->Write(); hUeZ[iJetBin][iSpecies]->Write(); hCorrZ[iJetBin][iSpecies]->Write();
      hRawJT[iJetBin][iSpecies]->Write(); hUeJT[iJetBin][iSpecies]->Write(); hCorrJT[iJetBin][iSpecies]->Write();
    }
  }
  out->Close();

  std::ofstream text(outputText);
  text << "# Stage 2 pp baseline summary\n";
  text << "# Existing ting status-1 particle jets; neutrinos retained in clustering by PI decision.\n";
  text << "# Species: pi=pi+ plus pi-, ka=K+ plus K-, pr=p plus pbar.\n";
  text << "# Analysis cuts: 4<=jetPt<6 or 6<=jetPt<9 GeV/c; |jetEta|<=0.9; z>=0.05; jT<(0.025+0.3295*z)*jetPt.\n";
  text << "# UE = arithmetic mean over stored off-axis cones; z,jT evaluated relative to original jet.\n";
  text << "files_read " << filesRead << "\n";
  text << "manifest_rows_used " << manifestRowsUsed << "\n";
  text << "max_file_pairs " << maxFilePairs << "\n";
  text << "files_failed " << filesFailed << "\n";
  text << "events_read " << eventsRead << "\n";
  text << "events_no_ue_vertex " << eventsNoUeVertex << "\n";
  text << "events_jet_ue_count_mismatch " << eventsJetUeCountMismatch << "\n";
  text << "selected_jets_no_cone " << selectedJetsNoCone << "\n";
  text << "selected_jets " << selectedJets << "\n";
  text << "# hardbin_low hardbin_high files entries expected_Ngenerated sigma_mb\n";
  for (Int_t ih = 0; ih < 11; ++ih) {
    text << kHardLow[ih] << " " << kHardHigh[ih] << " "
         << filesByHardBin[ih] << " " << entriesByHardBin[ih] << " "
         << kHardNGenerated[ih] << " " << std::setprecision(10)
         << kHardSigmaMb[ih] << "\n";
  }
  text << "# jetbin species observable bin_low bin_high raw ue corrected raw_err ue_err corrected_err\n";
  for (iJetBin = 0; iJetBin < kNJetBins; ++iJetBin) {
    for (iSpecies = 0; iSpecies < kNSpecies; ++iSpecies) {
      for (Int_t ib = 1; ib <= hCorrZ[iJetBin][iSpecies]->GetNbinsX(); ++ib) {
        text << kJetBinName[iJetBin] << " " << kSpecies[iSpecies] << " z "
             << hCorrZ[iJetBin][iSpecies]->GetXaxis()->GetBinLowEdge(ib) << " "
             << hCorrZ[iJetBin][iSpecies]->GetXaxis()->GetBinUpEdge(ib) << " "
             << hRawZ[iJetBin][iSpecies]->GetBinContent(ib) << " "
             << hUeZ[iJetBin][iSpecies]->GetBinContent(ib) << " "
             << hCorrZ[iJetBin][iSpecies]->GetBinContent(ib) << " "
             << hRawZ[iJetBin][iSpecies]->GetBinError(ib) << " "
             << hUeZ[iJetBin][iSpecies]->GetBinError(ib) << " "
             << hCorrZ[iJetBin][iSpecies]->GetBinError(ib) << "\n";
      }
      for (Int_t ib = 1; ib <= hCorrJT[iJetBin][iSpecies]->GetNbinsX(); ++ib) {
        text << kJetBinName[iJetBin] << " " << kSpecies[iSpecies] << " jT "
             << hCorrJT[iJetBin][iSpecies]->GetXaxis()->GetBinLowEdge(ib) << " "
             << hCorrJT[iJetBin][iSpecies]->GetXaxis()->GetBinUpEdge(ib) << " "
             << hRawJT[iJetBin][iSpecies]->GetBinContent(ib) << " "
             << hUeJT[iJetBin][iSpecies]->GetBinContent(ib) << " "
             << hCorrJT[iJetBin][iSpecies]->GetBinContent(ib) << " "
             << hRawJT[iJetBin][iSpecies]->GetBinError(ib) << " "
             << hUeJT[iJetBin][iSpecies]->GetBinError(ib) << " "
             << hCorrJT[iJetBin][iSpecies]->GetBinError(ib) << "\n";
      }
    }
  }
  text.close();

  std::cout << "Wrote " << outputRoot << " and " << outputText << std::endl;
  std::cout << "Read " << filesRead << " file pairs, selected " << selectedJets << " jets." << std::endl;
}

