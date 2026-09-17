#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TMethodCall.h"
#include <iostream>

void Stage2_FindUeParticle(char *ueFileName)
{
  gSystem->Load("libStJetSkimEvent");
  gSystem->Load("libStJets");
  gSystem->Load("libStJetEvent");
  gSystem->Load("libStUeEvent");

  TMethodCall numberOfVertices, vertexAt;
  TMethodCall numberOfUeJets, ueJetAt;
  TMethodCall numberOfCones, coneAt;
  TMethodCall numberOfParticles, particleAt;
  TMethodCall pdg, pt, eta, phi, status;

  numberOfVertices.InitWithPrototype(
      TClass::GetClass("StUeOffAxisConesEvent"),
      "numberOfVertices", "");

  vertexAt.InitWithPrototype(
      TClass::GetClass("StUeOffAxisConesEvent"),
      "vertex", "Int_t");

  numberOfUeJets.InitWithPrototype(
      TClass::GetClass("StUeVertex"),
      "numberOfUeJets", "");

  ueJetAt.InitWithPrototype(
      TClass::GetClass("StUeVertex"),
      "ueJet", "Int_t");

  numberOfCones.InitWithPrototype(
      TClass::GetClass("StUeOffAxisConesJet"),
      "numberOfCones", "");

  coneAt.InitWithPrototype(
      TClass::GetClass("StUeOffAxisConesJet"),
      "cone", "Int_t");

  numberOfParticles.InitWithPrototype(
      TClass::GetClass("StUeOffAxisCones"),
      "numberOfParticles", "");

  particleAt.InitWithPrototype(
      TClass::GetClass("StUeOffAxisCones"),
      "particle", "Int_t");

  pdg.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "pdg", "");

  pt.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "pt", "");

  eta.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "eta", "");

  phi.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "phi", "");

  status.InitWithPrototype(
      TClass::GetClass("StJetParticle"), "status", "");

  if (!numberOfVertices.IsValid() || !vertexAt.IsValid() ||
      !numberOfUeJets.IsValid() || !ueJetAt.IsValid() ||
      !numberOfCones.IsValid() || !coneAt.IsValid() ||
      !numberOfParticles.IsValid() || !particleAt.IsValid() ||
      !pdg.IsValid() || !pt.IsValid() || !eta.IsValid() ||
      !phi.IsValid() || !status.IsValid()) {
    std::cout << "ERROR: at least one reflection method is invalid."
              << std::endl;
    return;
  }

  TFile *file = TFile::Open(ueFileName, "READ");

  if (!file || file->IsZombie()) {
    std::cout << "ERROR: cannot open UE file." << std::endl;
    return;
  }

  TTree *tree = (TTree *)file->Get("ue");
  TBranch *branch =
      tree ? tree->GetBranch("AntiKtR060ParticleOffAxisConesR060") : 0;

  if (!tree || !branch) {
    std::cout << "ERROR: UE tree or branch is absent." << std::endl;
    return;
  }

  void *ueEvent = 0;
  branch->SetAddress(&ueEvent);

  Long64_t nEntries = tree->GetEntries();
  Long64_t nToScan = nEntries < 100 ? nEntries : 100;

  std::cout << "Scanning " << nToScan
            << " UE events out of " << nEntries << std::endl;

  for (Long64_t iEvent = 0; iEvent < nToScan; ++iEvent) {
    branch->GetEntry(iEvent);

    Long_t result = 0;
    numberOfVertices.Execute(ueEvent, result);
    Long_t nVertices = result;

    for (Long_t iVertex = 0; iVertex < nVertices; ++iVertex) {
      vertexAt.ResetParam();
      vertexAt.SetParam((Long_t)iVertex);
      vertexAt.Execute(ueEvent, result);
      void *vertex = (void *)result;

      numberOfUeJets.Execute(vertex, result);
      Long_t nUeJets = result;

      for (Long_t iUeJet = 0; iUeJet < nUeJets; ++iUeJet) {
        ueJetAt.ResetParam();
        ueJetAt.SetParam((Long_t)iUeJet);
        ueJetAt.Execute(vertex, result);
        void *ueJet = (void *)result;

        numberOfCones.Execute(ueJet, result);
        Long_t nCones = result;

        for (Long_t iCone = 0; iCone < nCones; ++iCone) {
          coneAt.ResetParam();
          coneAt.SetParam((Long_t)iCone);
          coneAt.Execute(ueJet, result);
          void *cone = (void *)result;

          numberOfParticles.Execute(cone, result);
          Long_t nParticles = result;

          if (nParticles <= 0) continue;

          particleAt.ResetParam();
          particleAt.SetParam((Long_t)0);
          particleAt.Execute(cone, result);
          void *particle = (void *)result;

          Long_t pdgCode = 0;
          Long_t statusCode = 0;
          Double_t particlePt = 0.0;
          Double_t particleEta = 0.0;
          Double_t particlePhi = 0.0;

          pdg.Execute(particle, pdgCode);
          status.Execute(particle, statusCode);
          pt.Execute(particle, particlePt);
          eta.Execute(particle, particleEta);
          phi.Execute(particle, particlePhi);

          std::cout << "FOUND non-empty UE cone" << std::endl;
          std::cout << "event=" << iEvent
                    << " vertex=" << iVertex
                    << " ueJet=" << iUeJet
                    << " cone=" << iCone
                    << " nParticles=" << nParticles
                    << std::endl;

          std::cout << "particle 0:"
                    << " PDG=" << pdgCode
                    << " status=" << statusCode
                    << " pT=" << particlePt
                    << " eta=" << particleEta
                    << " phi=" << particlePhi
                    << std::endl;

          file->Close();
          return;
        }
      }
    }
  }

  std::cout << "No non-empty UE cone was found in first "
            << nToScan << " events." << std::endl;

  file->Close();
}
