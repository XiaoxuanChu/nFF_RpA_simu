// Stage2_ProbeJetUeReflection.C
// Header-free probe for legacy STAR jet and UE-off-axis-cone ROOT trees.
// Run from xchu's own pp directory.  This macro never writes to input files.

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TMethodCall.h"
#include <iostream>

Bool_t Stage2ProbeInit(TMethodCall &call,
                       const char *className,
                       const char *methodName,
                       const char *prototype)
{
  TClass *cl = TClass::GetClass(className);
  if (!cl) {
    std::cout << "ERROR: cannot find class " << className << std::endl;
    return kFALSE;
  }

  call.InitWithPrototype(cl, methodName, prototype);
  if (!call.IsValid()) {
    std::cout << "ERROR: cannot resolve "
              << className << "::" << methodName
              << "(" << prototype << ")" << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

Long_t Stage2ProbeCallLong0(TMethodCall &call, void *object)
{
  Long_t result = 0;
  call.Execute(object, result);
  return result;
}

Long_t Stage2ProbeCallLong1(TMethodCall &call, void *object, Long_t index)
{
  Long_t result = 0;
  call.ResetParam();
  call.SetParam((Long_t)index);
  call.Execute(object, result);
  return result;
}

Double_t Stage2ProbeCallDouble0(TMethodCall &call, void *object)
{
  Double_t result = 0.0;
  call.Execute(object, result);
  return result;
}

void Stage2_ProbeJetUeReflection(
    const char *jetFileName,
    const char *ueFileName)
{
  gSystem->Load("libStJetSkimEvent");
  gSystem->Load("libStJets");
  gSystem->Load("libStJetEvent");
  gSystem->Load("libStUeEvent");

  TMethodCall numberOfJets, jetAt;
  TMethodCall jetPt, jetEta, numberOfJetParticles, jetParticleAt;
  TMethodCall particlePdg, particlePt, particleEta;
  TMethodCall numberOfVertices, vertexAt;
  TMethodCall numberOfUeJets, ueJetAt;
  TMethodCall ueJetPt, ueJetEta, numberOfCones, coneAt;
  TMethodCall numberOfConeParticles, coneParticleAt;

  Bool_t ok = kTRUE;

  ok = ok && Stage2ProbeInit(numberOfJets,
                             "StJetEvent", "numberOfJets", "");
  ok = ok && Stage2ProbeInit(jetAt,
                             "StJetEvent", "jet", "Int_t");

  ok = ok && Stage2ProbeInit(jetPt,
                             "StJetCandidate", "pt", "");
  ok = ok && Stage2ProbeInit(jetEta,
                             "StJetCandidate", "eta", "");
  ok = ok && Stage2ProbeInit(numberOfJetParticles,
                             "StJetCandidate", "numberOfParticles", "");
  ok = ok && Stage2ProbeInit(jetParticleAt,
                             "StJetCandidate", "particle", "Int_t");

  ok = ok && Stage2ProbeInit(particlePdg,
                             "StJetParticle", "pdg", "");
  ok = ok && Stage2ProbeInit(particlePt,
                             "StJetParticle", "pt", "");
  ok = ok && Stage2ProbeInit(particleEta,
                             "StJetParticle", "eta", "");

  ok = ok && Stage2ProbeInit(numberOfVertices,
                             "StUeOffAxisConesEvent",
                             "numberOfVertices", "");
  ok = ok && Stage2ProbeInit(vertexAt,
                             "StUeOffAxisConesEvent",
                             "vertex", "Int_t");

  ok = ok && Stage2ProbeInit(numberOfUeJets,
                             "StUeVertex", "numberOfUeJets", "");
  ok = ok && Stage2ProbeInit(ueJetAt,
                             "StUeVertex", "ueJet", "Int_t");

  ok = ok && Stage2ProbeInit(ueJetPt,
                             "StUeJet", "pt", "");
  ok = ok && Stage2ProbeInit(ueJetEta,
                             "StUeJet", "eta", "");

  ok = ok && Stage2ProbeInit(numberOfCones,
                             "StUeOffAxisConesJet",
                             "numberOfCones", "");
  ok = ok && Stage2ProbeInit(coneAt,
                             "StUeOffAxisConesJet", "cone", "Int_t");

  ok = ok && Stage2ProbeInit(numberOfConeParticles,
                             "StUeOffAxisCones",
                             "numberOfParticles", "");
  ok = ok && Stage2ProbeInit(coneParticleAt,
                             "StUeOffAxisCones", "particle", "Int_t");

  if (!ok) {
    std::cout << "ERROR: reflection API setup failed." << std::endl;
    return;
  }

  TFile *jetFile = TFile::Open(jetFileName, "READ");
  TFile *ueFile  = TFile::Open(ueFileName, "READ");

  if (!jetFile || jetFile->IsZombie() ||
      !ueFile  || ueFile->IsZombie()) {
    std::cout << "ERROR: cannot open one or both input files."
              << std::endl;
    return;
  }

  TTree *jetTree = (TTree *)jetFile->Get("JetTree");
  TTree *ueTree  = (TTree *)ueFile->Get("JetTree");

  if (!jetTree || !ueTree) {
    std::cout << "ERROR: JetTree is absent in one or both files."
              << std::endl;
    return;
  }

  TBranch *jetBranch =
      jetTree->GetBranch("AntiKtR060Particle");
  TBranch *ueBranch =
      ueTree->GetBranch("AntiKtR060ParticleOffAxisConesR060");

  if (!jetBranch || !ueBranch) {
    std::cout << "ERROR: required branch is absent." << std::endl;
    std::cout << "  jet branch found: " << (jetBranch ? "yes" : "no")
              << std::endl;
    std::cout << "  UE branch found:  " << (ueBranch ? "yes" : "no")
              << std::endl;
    return;
  }

  void *jetEvent = 0;
  void *ueEvent = 0;

  jetBranch->SetAddress(&jetEvent);
  ueBranch->SetAddress(&ueEvent);

  jetBranch->GetEntry(0);
  ueBranch->GetEntry(0);

  if (!jetEvent || !ueEvent) {
    std::cout << "ERROR: event 0 could not be read." << std::endl;
    return;
  }

  Long_t nJets = Stage2ProbeCallLong0(numberOfJets, jetEvent);
  std::cout << "event 0: nJets = " << nJets << std::endl;

  if (nJets > 0) {
    void *jet0 = (void *)Stage2ProbeCallLong1(jetAt, jetEvent, (Long_t)0);

    if (jet0) {
      Long_t nParticles =
          Stage2ProbeCallLong0(numberOfJetParticles, jet0);

      std::cout << "jet 0:"
                << " pT = " << Stage2ProbeCallDouble0(jetPt, jet0)
                << " eta = " << Stage2ProbeCallDouble0(jetEta, jet0)
                << " nParticles = " << nParticles
                << std::endl;

      if (nParticles > 0) {
        void *particle0 = (void *)Stage2ProbeCallLong1(
            jetParticleAt, jet0, (Long_t)0);

        if (particle0) {
          std::cout << "jet 0, particle 0:"
                    << " PDG = "
                    << Stage2ProbeCallLong0(particlePdg, particle0)
                    << " pT = "
                    << Stage2ProbeCallDouble0(particlePt, particle0)
                    << " eta = "
                    << Stage2ProbeCallDouble0(particleEta, particle0)
                    << std::endl;
        }
      }
    }
  }

  Long_t nVertices =
      Stage2ProbeCallLong0(numberOfVertices, ueEvent);

  std::cout << "UE event 0: nVertices = "
            << nVertices << std::endl;

  if (nVertices > 0) {
    void *vertex0 = (void *)Stage2ProbeCallLong1(
        vertexAt, ueEvent, (Long_t)0);

    if (vertex0) {
      Long_t nUeJets =
          Stage2ProbeCallLong0(numberOfUeJets, vertex0);

      std::cout << "UE vertex 0: nUeJets = "
                << nUeJets << std::endl;

      if (nUeJets > 0) {
        void *ueJet0 = (void *)Stage2ProbeCallLong1(
            ueJetAt, vertex0, (Long_t)0);

        if (ueJet0) {
          Long_t nCone =
              Stage2ProbeCallLong0(numberOfCones, ueJet0);

          std::cout << "UE jet 0:"
                    << " pT = "
                    << Stage2ProbeCallDouble0(ueJetPt, ueJet0)
                    << " eta = "
                    << Stage2ProbeCallDouble0(ueJetEta, ueJet0)
                    << " nCones = " << nCone
                    << std::endl;

          if (nCone > 0) {
            void *cone0 = (void *)Stage2ProbeCallLong1(
                coneAt, ueJet0, (Long_t)0);

            if (cone0) {
              Long_t nConeParticles =
                  Stage2ProbeCallLong0(numberOfConeParticles, cone0);

              std::cout << "UE jet 0, cone 0: nParticles = "
                        << nConeParticles << std::endl;

              if (nConeParticles > 0) {
                void *coneParticle0 =
                    (void *)Stage2ProbeCallLong1(
                        coneParticleAt, cone0, (Long_t)0);

                if (coneParticle0) {
                  std::cout << "UE jet 0, cone 0, particle 0:"
                            << " PDG = "
                            << Stage2ProbeCallLong0(
                                   particlePdg, coneParticle0)
                            << " pT = "
                            << Stage2ProbeCallDouble0(
                                   particlePt, coneParticle0)
                            << " eta = "
                            << Stage2ProbeCallDouble0(
                                   particleEta, coneParticle0)
                            << std::endl;
                }
              }
            }
          }
        }
      }
    }
  }

  jetFile->Close();
  ueFile->Close();

  std::cout << "Stage2 probe completed." << std::endl;
}
