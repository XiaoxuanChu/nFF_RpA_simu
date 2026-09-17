#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TMethodCall.h"
#include <iostream>

void Stage2_CheckNeutrinoContent(const char *jetFileName)
{
  gSystem->Load("libStJetSkimEvent");
  gSystem->Load("libStJets");
  gSystem->Load("libStJetEvent");
  gSystem->Load("libStUeEvent");

  TMethodCall nJetsCall, jetAtCall;
  TMethodCall jetPtCall, jetEtaCall;
  TMethodCall nParticlesCall, particleAtCall;
  TMethodCall particlePdgCall, particlePtCall;

  nJetsCall.InitWithPrototype(
      TClass::GetClass("StJetEvent"), "numberOfJets", "");
  jetAtCall.InitWithPrototype(
      TClass::GetClass("StJetEvent"), "jet", "Int_t");

  jetPtCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "pt", "");
  jetEtaCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "eta", "");
  nParticlesCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "numberOfParticles", "");
  particleAtCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "particle", "Int_t");

  particlePdgCall.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "pdg", "");
  particlePtCall.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "pt", "");

  if (!nJetsCall.IsValid() || !jetAtCall.IsValid() ||
      !jetPtCall.IsValid() || !jetEtaCall.IsValid() ||
      !nParticlesCall.IsValid() || !particleAtCall.IsValid() ||
      !particlePdgCall.IsValid() || !particlePtCall.IsValid()) {
    std::cout << "ERROR: reflection setup failed." << std::endl;
    return;
  }

  TFile *file = TFile::Open(jetFileName, "READ");

  if (!file || file->IsZombie()) {
    std::cout << "ERROR: cannot open jet file." << std::endl;
    return;
  }

  TTree *tree = (TTree *)file->Get("jet");
  TBranch *branch =
      tree ? tree->GetBranch("AntiKtR060Particle") : 0;

  if (!tree || !branch) {
    std::cout << "ERROR: jet tree or branch is absent." << std::endl;
    return;
  }

  void *jetEvent = 0;
  branch->SetAddress(&jetEvent);

  Long64_t nEntries = tree->GetEntries();
  Long_t acceptedJets = 0;
  Long_t jetsWithNeutrino = 0;
  Long_t neutrinoParticles = 0;
  Double_t neutrinoPtSum = 0.0;
  Double_t largestNeutrinoPt = 0.0;

  for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
    branch->GetEntry(iEvent);

    Long_t result = 0;
    nJetsCall.Execute(jetEvent, result);
    Long_t nJets = result;

    for (Long_t iJet = 0; iJet < nJets; ++iJet) {
      jetAtCall.ResetParam();
      jetAtCall.SetParam((Long_t)iJet);
      jetAtCall.Execute(jetEvent, result);
      void *jet = (void *)result;

      Double_t jetPt = 0.0;
      Double_t jetEta = 0.0;
      jetPtCall.Execute(jet, jetPt);
      jetEtaCall.Execute(jet, jetEta);

      if (jetPt < 4.0 || jetPt >= 9.0) continue;
      if (jetEta < -0.9 || jetEta > 0.9) continue;

      ++acceptedJets;

      nParticlesCall.Execute(jet, result);
      Long_t nParticles = result;
      Bool_t thisJetHasNeutrino = kFALSE;

      for (Long_t iParticle = 0;
           iParticle < nParticles; ++iParticle) {
        particleAtCall.ResetParam();
        particleAtCall.SetParam((Long_t)iParticle);
        particleAtCall.Execute(jet, result);
        void *particle = (void *)result;

        Long_t pdgCode = 0;
        Double_t particlePt = 0.0;
        particlePdgCall.Execute(particle, pdgCode);
        particlePtCall.Execute(particle, particlePt);

        if (pdgCode < 0) pdgCode = -pdgCode;

        if (pdgCode == 12 || pdgCode == 14 || pdgCode == 16) {
          ++neutrinoParticles;
          neutrinoPtSum += particlePt;
          if (particlePt > largestNeutrinoPt)
            largestNeutrinoPt = particlePt;
          thisJetHasNeutrino = kTRUE;
        }
      }

      if (thisJetHasNeutrino) ++jetsWithNeutrino;
    }
  }

  std::cout << "entries scanned = " << nEntries << std::endl;
  std::cout << "accepted jets (4 <= pT < 9, |eta| <= 0.9) = "
            << acceptedJets << std::endl;
  std::cout << "accepted jets containing >=1 neutrino = "
            << jetsWithNeutrino << std::endl;
  std::cout << "neutrino constituents = "
            << neutrinoParticles << std::endl;
  std::cout << "sum neutrino pT = "
            << neutrinoPtSum << " GeV/c" << std::endl;
  std::cout << "largest neutrino pT = "
            << largestNeutrinoPt << " GeV/c" << std::endl;

  file->Close();
}
