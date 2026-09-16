// ROOT5/STAR macro: particle-level pp hadron-in-jet baseline for one matched
// jet/UE file pair.  It reads ONLY AntiKtR060Particle and its matching
// AntiKtR060ParticleOffAxisConesR060 branch.
//
// Run in the legacy STAR environment selected in Stage 2:
// root4star -b -q 'Stage2_pp_Baseline.C("../pp/Input/pt5_7_16092051_1.jets.root","../pp/Input/pt5_7_16092051_1.ueoc.root","Stage2_pp_pt5_7_test.root",1.0)'
//
// The supplied 196-event 5--7 GeV/c sample is structural validation only.
// Its old jet construction contains status-1 neutrinos and stores jets only
// above 6 GeV/c; it cannot be accepted as the final 4--6 GeV/c baseline.

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "StJetEvent.h"
#include "StJetCandidate.h"
#include "StJetParticle.h"
#include "StUeVertex.h"
#include "StUeJet.h"
#include "StUeOffAxisConesEvent.h"
#include "StUeOffAxisConesJet.h"
#include "StUeOffAxisCones.h"

namespace Stage2PpBaseline {
const Int_t kNJetBins = 2;
const Int_t kNSpecies = 3;
const Double_t kJetMin[kNJetBins] = {4.0, 6.0};
const Double_t kJetMax[kNJetBins] = {6.0, 9.0};
const char* kSpecies[kNSpecies] = {"pi", "ka", "pr"};
const char* kSpeciesLabel[kNSpecies] = {"#pi^{+}+#pi^{-}", "K^{+}+K^{-}", "p+#bar{p}"};
const Double_t kZEdges[] = {0.00, 0.05, 0.15, 0.30, 0.50, 0.70, 0.90, 1.00};
const Double_t kJtEdges[] = {0.00, 0.10, 0.20, 0.40, 0.60, 0.90, 1.00};

Bool_t IsNeutrino(Int_t pdg) {
  const Int_t a = TMath::Abs(pdg);
  return a == 12 || a == 14 || a == 16;
}

Int_t SpeciesIndex(Int_t pdg) {
  const Int_t a = TMath::Abs(pdg);
  if (a == 211) return 0;
  if (a == 321) return 1;
  if (a == 2212) return 2;
  return -1;
}

Int_t JetBin(Double_t pt, Double_t eta) {
  if (TMath::Abs(eta) > 0.9) return -1;
  for (Int_t i = 0; i < kNJetBins; ++i)
    if (pt >= kJetMin[i] && pt < kJetMax[i]) return i;
  return -1;
}

Double_t DeltaPhi(Double_t a, Double_t b) {
  Double_t d = a - b;
  while (d > TMath::Pi()) d -= 2.0*TMath::Pi();
  while (d <= -TMath::Pi()) d += 2.0*TMath::Pi();
  return d;
}

TVector3 Momentum(Double_t pt, Double_t eta, Double_t phi) {
  return TVector3(pt*TMath::Cos(phi), pt*TMath::Sin(phi), pt*TMath::SinH(eta));
}

Bool_t PassCommonCuts(StJetParticle* p, const TVector3& jetMomentum,
                      Double_t& z, Double_t& jt)
{
  if (!p || p->status() != 1 || IsNeutrino(p->pdg())) return kFALSE;
  const TVector3 hMomentum = Momentum(p->pt(), p->eta(), p->phi());
  const Double_t jetP = jetMomentum.Mag();
  if (jetP <= 0.0) return kFALSE;
  z = hMomentum.Mag()/jetP;
  jt = hMomentum.Cross(jetMomentum).Mag()/jetP;
  return z >= 0.05 && jt < (0.025 + 0.3295*z)*jetP;
}

StUeOffAxisConesJet* MatchUeJet(StUeOffAxisConesEvent* ueEvent,
                                StJetCandidate* jet, Double_t& matchDr)
{
  StUeOffAxisConesJet* best = 0;
  matchDr = 1.e9;
  if (!ueEvent || !jet) return best;
  for (Int_t iv = 0; iv < ueEvent->numberOfVertices(); ++iv) {
    StUeVertex* vertex = ueEvent->vertex(iv);
    if (!vertex) continue;
    for (Int_t iu = 0; iu < vertex->numberOfUeJets(); ++iu) {
      StUeJet* candidate = vertex->ueJet(iu);
      StUeOffAxisConesJet* ueJet = dynamic_cast<StUeOffAxisConesJet*>(candidate);
      if (!ueJet) continue;
      const Double_t dEta = ueJet->eta() - jet->eta();
      const Double_t dPhi = DeltaPhi(ueJet->phi(), jet->phi());
      const Double_t dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
      if (dR < matchDr) { matchDr = dR; best = ueJet; }
    }
  }
  return best;
}

void FillParticle(StJetParticle* particle, const TVector3& jetMomentum,
                  TH1D* hZ[kNSpecies], TH1D* hJt[kNSpecies], Double_t weight)
{
  Double_t z = 0.0, jt = 0.0;
  if (!PassCommonCuts(particle, jetMomentum, z, jt)) return;
  const Int_t species = SpeciesIndex(particle->pdg());
  if (species < 0) return;
  hZ[species]->Fill(z, weight);
  hJt[species]->Fill(jt, weight);
}
}

void Stage2_pp_Baseline(const char* jetFile,
                        const char* ueFile,
                        const char* outputFile = "Stage2_pp_Baseline.root",
                        Double_t eventWeight = 1.0)
{
  using namespace Stage2PpBaseline;
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");
  if (gSystem->Load("StJetSkimEvent") < 0 || gSystem->Load("StJets") < 0 ||
      gSystem->Load("StJetEvent") < 0 || gSystem->Load("StUeEvent") < 0) {
    Error("Stage2_pp_Baseline", "Cannot load required persisted-object libraries.");
    return;
  }

  TFile jetInput(jetFile, "READ");
  TFile ueInput(ueFile, "READ");
  TTree* jetTree = dynamic_cast<TTree*>(jetInput.Get("jet"));
  TTree* ueTree = dynamic_cast<TTree*>(ueInput.Get("ue"));
  if (!jetTree || !ueTree) {
    Error("Stage2_pp_Baseline", "Required jet or UE tree is missing.");
    return;
  }
  const char* jetBranchName = "AntiKtR060Particle";
  const char* ueBranchName = "AntiKtR060ParticleOffAxisConesR060";
  if (!jetTree->GetBranch(jetBranchName) || !ueTree->GetBranch(ueBranchName)) {
    Error("Stage2_pp_Baseline", "Required R=0.6 particle branch is missing.");
    return;
  }
  if (jetTree->GetEntries() != ueTree->GetEntries()) {
    Error("Stage2_pp_Baseline", "Jet and UE entry counts differ: %lld versus %lld.",
          jetTree->GetEntries(), ueTree->GetEntries());
    return;
  }

  StJetEvent* jetEvent = 0;
  StUeOffAxisConesEvent* ueEvent = 0;
  jetTree->SetBranchAddress(jetBranchName, &jetEvent);
  ueTree->SetBranchAddress(ueBranchName, &ueEvent);

  TH1D* hRawZ[kNJetBins][kNSpecies]; TH1D* hRawJt[kNJetBins][kNSpecies];
  TH1D* hUeZ[kNJetBins][kNSpecies];  TH1D* hUeJt[kNJetBins][kNSpecies];
  for (Int_t ib = 0; ib < kNJetBins; ++ib) for (Int_t is = 0; is < kNSpecies; ++is) {
    hRawZ[ib][is] = new TH1D(Form("hRawZ_%s_jet%d",kSpecies[is],ib),"",7,kZEdges);
    hRawJt[ib][is] = new TH1D(Form("hRawJt_%s_jet%d",kSpecies[is],ib),"",6,kJtEdges);
    hUeZ[ib][is] = new TH1D(Form("hUeZ_%s_jet%d",kSpecies[is],ib),"",7,kZEdges);
    hUeJt[ib][is] = new TH1D(Form("hUeJt_%s_jet%d",kSpecies[is],ib),"",6,kJtEdges);
    hRawZ[ib][is]->Sumw2(); hRawJt[ib][is]->Sumw2();
    hUeZ[ib][is]->Sumw2(); hUeJt[ib][is]->Sumw2();
  }

  Double_t nJetSelected[kNJetBins] = {0.0,0.0};
  Double_t nJetMatched[kNJetBins] = {0.0,0.0};
  Long64_t nNeutrinoConstituents = 0, nUnmatched = 0;
  Double_t maxUeJetMatchDr = 0.0;
  const Long64_t nEntries = jetTree->GetEntries();
  for (Long64_t ie = 0; ie < nEntries; ++ie) {
    jetTree->GetEntry(ie); ueTree->GetEntry(ie);
    if (!jetEvent || !ueEvent) { Error("Stage2_pp_Baseline", "Null event at %lld",ie); return; }
    for (Int_t ij = 0; ij < jetEvent->numberOfJets(); ++ij) {
      StJetCandidate* jet = jetEvent->jet(ij);
      if (!jet) continue;
      const Int_t jetBin = JetBin(jet->pt(), jet->eta());
      if (jetBin < 0) continue;
      nJetSelected[jetBin] += eventWeight;
      Double_t matchDr = 0.0;
      StUeOffAxisConesJet* ueJet = MatchUeJet(ueEvent, jet, matchDr);
      if (!ueJet || matchDr > 1.e-4) { ++nUnmatched; continue; }
      if (matchDr > maxUeJetMatchDr) maxUeJetMatchDr = matchDr;
      const Int_t nCones = ueJet->numberOfCones();
      if (nCones <= 0) { ++nUnmatched; continue; }
      nJetMatched[jetBin] += eventWeight;
      const TVector3 jetMomentum = Momentum(jet->pt(), jet->eta(), jet->phi());
      for (Int_t ip = 0; ip < jet->numberOfParticles(); ++ip) {
        StJetParticle* particle = jet->particle(ip);
        if (particle && IsNeutrino(particle->pdg())) ++nNeutrinoConstituents;
        FillParticle(particle, jetMomentum, hRawZ[jetBin], hRawJt[jetBin], eventWeight);
      }
      const Double_t coneWeight = eventWeight/Double_t(nCones);
      for (Int_t ic = 0; ic < nCones; ++ic) {
        StUeOffAxisCones* cone = ueJet->cone(ic);
        if (!cone) continue;
        for (Int_t ip = 0; ip < cone->numberOfParticles(); ++ip)
          FillParticle(cone->particle(ip), jetMomentum, hUeZ[jetBin], hUeJt[jetBin], coneWeight);
      }
    }
  }

  TFile output(outputFile, "RECREATE");
  TParameter<Long64_t>("input_entries",nEntries).Write();
  TParameter<Double_t>("event_weight",eventWeight).Write();
  TParameter<Long64_t>("unmatched_ue_jets",nUnmatched).Write();
  TParameter<Long64_t>("neutrino_constituents_in_input_jets",nNeutrinoConstituents).Write();
  TParameter<Double_t>("maximum_ue_jet_match_dR",maxUeJetMatchDr).Write();
  for (Int_t ib = 0; ib < kNJetBins; ++ib) {
    TParameter<Double_t>(Form("njet_selected_%d",ib),nJetSelected[ib]).Write();
    TParameter<Double_t>(Form("njet_matched_%d",ib),nJetMatched[ib]).Write();
    for (Int_t is = 0; is < kNSpecies; ++is) {
      TH1D* hCorrZ = static_cast<TH1D*>(hRawZ[ib][is]->Clone(Form("DcorrZ_%s_jet%d",kSpecies[is],ib)));
      TH1D* hCorrJt = static_cast<TH1D*>(hRawJt[ib][is]->Clone(Form("DcorrJt_%s_jet%d",kSpecies[is],ib)));
      hCorrZ->Add(hUeZ[ib][is],-1.0); hCorrJt->Add(hUeJt[ib][is],-1.0);
      if (nJetMatched[ib] > 0.0) {
        hRawZ[ib][is]->Scale(1.0/nJetMatched[ib],"width"); hUeZ[ib][is]->Scale(1.0/nJetMatched[ib],"width"); hCorrZ->Scale(1.0/nJetMatched[ib],"width");
        hRawJt[ib][is]->Scale(1.0/nJetMatched[ib],"width"); hUeJt[ib][is]->Scale(1.0/nJetMatched[ib],"width"); hCorrJt->Scale(1.0/nJetMatched[ib],"width");
      }
      hRawZ[ib][is]->Write(); hUeZ[ib][is]->Write(); hCorrZ->Write();
      hRawJt[ib][is]->Write(); hUeJt[ib][is]->Write(); hCorrJt->Write();
    }
  }
  output.Close();
  printf("Stage2 pp baseline wrote %s\n",outputFile);
  printf("Selected/matched weighted jets: [4,6) %.6g/%.6g, [6,9) %.6g/%.6g\n",nJetSelected[0],nJetMatched[0],nJetSelected[1],nJetMatched[1]);
  printf("Unmatched UE jets=%lld, max match dR=%.3e, neutrino constituents in input jets=%lld\n",nUnmatched,maxUeJetMatchDr,nNeutrinoConstituents);
}

