// Counts UE cone occupancy for STAR-selected jets without interpreting UE particles.
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClass.h"
#include "TMethodCall.h"
#include <iostream>

Bool_t S2CUInit(TMethodCall& c,const char* cl,const char* m,const char* p) {
  TClass* x=TClass::GetClass(cl); if(!x) return kFALSE;
  c.InitWithPrototype(x,m,p); return c.IsValid();
}
Long_t S2CULong(TMethodCall& c,void* o) { Long_t r=0;c.Execute(o,r);return r; }
void* S2CUObject(TMethodCall& c,void* o,Long_t i) {
  Long_t r=0;c.ResetParam();c.SetParam((Long_t)i);c.Execute(o,r);return (void*)r;
}
Double_t S2CUDouble(TMethodCall& c,void* o) { Double_t r=0;c.Execute(o,r);return r; }

void Stage2_CountSelectedUe(const char* jetName,const char* ueName) {
  gSystem->Load("libStJetSkimEvent"); gSystem->Load("libStJets");
  gSystem->Load("libStJetEvent"); gSystem->Load("libStUeEvent");
  TMethodCall njet,jetat,jpt,jeta,nv,vat,nue,ueat,ncone,coneat,npart;
  Bool_t ok=kTRUE;
  ok=ok&&S2CUInit(njet,"StJetEvent","numberOfJets","");
  ok=ok&&S2CUInit(jetat,"StJetEvent","jet","Int_t");
  ok=ok&&S2CUInit(jpt,"StJetCandidate","pt","");
  ok=ok&&S2CUInit(jeta,"StJetCandidate","eta","");
  ok=ok&&S2CUInit(nv,"StUeOffAxisConesEvent","numberOfVertices","");
  ok=ok&&S2CUInit(vat,"StUeOffAxisConesEvent","vertex","Int_t");
  ok=ok&&S2CUInit(nue,"StUeVertex","numberOfUeJets","");
  ok=ok&&S2CUInit(ueat,"StUeVertex","ueJet","Int_t");
  ok=ok&&S2CUInit(ncone,"StUeOffAxisConesJet","numberOfCones","");
  ok=ok&&S2CUInit(coneat,"StUeOffAxisConesJet","cone","Int_t");
  ok=ok&&S2CUInit(npart,"StUeOffAxisCones","numberOfParticles","");
  if(!ok) { std::cout<<"ERROR: reflection setup failed"<<std::endl; return; }
  TFile* fj=TFile::Open(jetName,"READ"); TFile* fu=TFile::Open(ueName,"READ");
  TTree* tj=fj?(TTree*)fj->Get("jet"):0; TTree* tu=fu?(TTree*)fu->Get("ue"):0;
  TBranch* bj=tj?tj->GetBranch("AntiKtR060Particle"):0;
  TBranch* bu=tu?tu->GetBranch("AntiKtR060ParticleOffAxisConesR060"):0;
  if(!bj||!bu) {std::cout<<"ERROR: missing input branch"<<std::endl;return;}
  void* je=0;void* ue=0;bj->SetAddress(&je);bu->SetAddress(&ue);
  Long64_t n=tj->GetEntries(); if(tu->GetEntries()<n)n=tu->GetEntries();
  Long64_t jets=0,cones=0,nonempty=0,particles=0;
  for(Long64_t e=0;e<n;++e) {
    bj->GetEntry(e);bu->GetEntry(e); if(S2CULong(nv,ue)<1)continue;
    void* v=S2CUObject(vat,ue,0);Long_t nj=S2CULong(njet,je),nu=S2CULong(nue,v);
    if(nu<nj)nj=nu;
    for(Long_t j=0;j<nj;++j) {
      void* x=S2CUObject(jetat,je,j);Double_t pt=S2CUDouble(jpt,x),eta=S2CUDouble(jeta,x);
      if(pt<4||pt>=9||eta<-0.9||eta>0.9)continue; ++jets;
      void* u=S2CUObject(ueat,v,j);Long_t nc=S2CULong(ncone,u);cones+=nc;
      for(Long_t c=0;c<nc;++c){void* q=S2CUObject(coneat,u,c);Long_t np=S2CULong(npart,q);particles+=np;if(np>0)++nonempty;}
    }
  }
  std::cout<<"selected_jets "<<jets<<"\ncones "<<cones<<"\nnonempty_cones "<<nonempty<<"\ncone_particles "<<particles<<std::endl;
  fj->Close();fu->Close();
}

