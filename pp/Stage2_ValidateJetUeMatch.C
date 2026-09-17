#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TMethodCall.h"
#include <iostream>

void Stage2_ValidateJetUeMatch(const char *jetFileName,
                               const char *ueFileName)
{
  gSystem->Load("libStJetSkimEvent");
  gSystem->Load("libStJets");
  gSystem->Load("libStJetEvent");
  gSystem->Load("libStUeEvent");

  TMethodCall nJetsCall, jetAtCall;
  TMethodCall jetPtCall, jetEtaCall, jetPhiCall;

  TMethodCall nVerticesCall, vertexAtCall;
  TMethodCall nUeJetsCall, ueJetAtCall;
  TMethodCall ueJetPtCall, ueJetEtaCall, ueJetPhiCall;

  nJetsCall.InitWithPrototype(
      TClass::GetClass("StJetEvent"), "numberOfJets", "");
  jetAtCall.InitWithPrototype(
      TClass::GetClass("StJetEvent"), "jet", "Int_t");

  jetPtCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "pt", "");
  jetEtaCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "eta", "");
  jetPhiCall.InitWithPrototype(
      TClass::GetClass("StJetCandidate"), "phi", "");

  nVerticesCall.InitWithPrototype(
      TClass::GetClass("StUeOffAxisConesEvent"),
      "numberOfVertices", "");
  vertexAtCall.InitWithPrototype(
      TClass::GetClass("StUeOffAxisConesEvent"),
      "vertex", "Int_t");

  nUeJetsCall.InitWithPrototype(
      TClass::GetClass("StUeVertex"), "numberOfUeJets", "");
  ueJetAtCall.InitWithPrototype(
      TClass::GetClass("StUeVertex"), "ueJet", "Int_t");

  ueJetPtCall.InitWithPrototype(
      TClass::GetClass("StUeJet"), "pt", "");
  ueJetEtaCall.InitWithPrototype(
      TClass::GetClass("StUeJet"), "eta", "");
  ueJetPhiCall.InitWithPrototype(
      TClass::GetClass("StUeJet"), "phi", "");

  if (!nJetsCall.IsValid() || !jetAtCall.IsValid() ||
      !jetPtCall.IsValid() || !jetEtaCall.IsValid() ||
      !jetPhiCall.IsValid() || !nVerticesCall.IsValid() ||
      !vertexAtCall.IsValid() || !nUeJetsCall.IsValid() ||
      !ueJetAtCall.IsValid() || !ueJetPtCall.IsValid() ||
      !ueJetEtaCall.IsValid() || !ueJetPhiCall.IsValid()) {
    std::cout << "ERROR: reflection setup failed." << std::endl;
    return;
  }

  TFile *jetFile = TFile::Open(jetFileName, "READ");
  TFile *ueFile = TFile::Open(ueFileName, "READ");

  if (!jetFile || jetFile->IsZombie() ||
      !ueFile || ueFile->IsZombie()) {
    std::cout << "ERROR: cannot open one or both files." << std::endl;
    return;
  }

  TTree *jetTree = (TTree *)jetFile->Get("jet");
  TTree *ueTree = (TTree *)ueFile->Get("ue");

  TBranch *jetBranch =
      jetTree ? jetTree->GetBranch("AntiKtR060Particle") : 0;
  TBranch *ueBranch =
      ueTree ? ueTree->GetBranch(
          "AntiKtR060ParticleOffAxisConesR060") : 0;

  if (!jetBranch || !ueBranch) {
    std::cout << "ERROR: required branch is absent." << std::endl;
    return;
  }

  void *jetEvent = 0;
  void *ueEvent = 0;

  jetBranch->SetAddress(&jetEvent);
  ueBranch->SetAddress(&ueEvent);

  Long64_t nEntries = jetTree->GetEntries();
  if (ueTree->GetEntries() < nEntries) nEntries = ueTree->GetEntries();
  if (nEntries > 100) nEntries = 100;

  Long_t totalJets = 0;
  Long_t unequalJetCounts = 0;
  Long_t noVertexEvents = 0;
  Double_t maxDeltaPt = 0.0;
  Double_t maxDeltaEta = 0.0;
  Double_t maxDeltaPhi = 0.0;

  for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
    jetBranch->GetEntry(iEvent);
    ueBranch->GetEntry(iEvent);

    Long_t result = 0;

    nJetsCall.Execute(jetEvent, result);
    Long_t nJets = result;

    nVerticesCall.Execute(ueEvent, result);
    Long_t nVertices = result;

    if (nVertices < 1) {
      ++noVertexEvents;
      continue;
    }

    vertexAtCall.ResetParam();
    vertexAtCall.SetParam((Long_t)0);
    vertexAtCall.Execute(ueEvent, result);
    void *vertex = (void *)result;

    nUeJetsCall.Execute(vertex, result);
    Long_t nUeJets = result;

    if (nJets != nUeJets) ++unequalJetCounts;

    Long_t nCompare = nJets;
    if (nUeJets < nCompare) nCompare = nUeJets;

    for (Long_t iJet = 0; iJet < nCompare; ++iJet) {
      jetAtCall.ResetParam();
      jetAtCall.SetParam((Long_t)iJet);
      jetAtCall.Execute(jetEvent, result);
      void *jet = (void *)result;

      ueJetAtCall.ResetParam();
      ueJetAtCall.SetParam((Long_t)iJet);
      ueJetAtCall.Execute(vertex, result);
      void *ueJet = (void *)result;

      Double_t jetPt = 0.0;
      Double_t jetEta = 0.0;
      Double_t jetPhi = 0.0;
      Double_t uePt = 0.0;
      Double_t ueEta = 0.0;
      Double_t uePhi = 0.0;

      jetPtCall.Execute(jet, jetPt);
      jetEtaCall.Execute(jet, jetEta);
      jetPhiCall.Execute(jet, jetPhi);

      ueJetPtCall.Execute(ueJet, uePt);
      ueJetEtaCall.Execute(ueJet, ueEta);
      ueJetPhiCall.Execute(ueJet, uePhi);

      Double_t deltaPt = jetPt - uePt;
      Double_t deltaEta = jetEta - ueEta;
      Double_t deltaPhi = jetPhi - uePhi;

      if (deltaPt < 0.0) deltaPt = -deltaPt;
      if (deltaEta < 0.0) deltaEta = -deltaEta;
      if (deltaPhi < 0.0) deltaPhi = -deltaPhi;

      if (deltaPt > maxDeltaPt) maxDeltaPt = deltaPt;
      if (deltaEta > maxDeltaEta) maxDeltaEta = deltaEta;
      if (deltaPhi > maxDeltaPhi) maxDeltaPhi = deltaPhi;

      ++totalJets;
    }
  }

  std::cout << "events scanned = " << nEntries << std::endl;
  std::cout << "compared jets = " << totalJets << std::endl;
  std::cout << "events with no UE vertex = "
            << noVertexEvents << std::endl;
  std::cout << "events with nJets != nUeJets = "
            << unequalJetCounts << std::endl;
  std::cout << "maximum |delta pT|  = "
            << maxDeltaPt << std::endl;
  std::cout << "maximum |delta eta| = "
            << maxDeltaEta << std::endl;
  std::cout << "maximum |delta phi| = "
            << maxDeltaPhi << std::endl;

  jetFile->Close();
  ueFile->Close();
}
