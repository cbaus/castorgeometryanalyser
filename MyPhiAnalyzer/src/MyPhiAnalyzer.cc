// -*- C++ -*-
//
// Package:    MyPhiAnalyzer
// Class:      MyPhiAnalyzer
//
/**\class MyPhiAnalyzer MyPhiAnalyzer.cc Geometry/MyPhiAnalyzer/src/MyPhiAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Colin Baus, 133, 26229
//         Created:  Fri Sep  7 18:48:50 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Castor

#include "CondFormats/CastorObjects/interface/AllObjects.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "CalibFormats/CastorObjects/interface/CastorDbRecord.h"
#include "CalibFormats/CastorObjects/interface/CastorDbService.h"
#include "CalibFormats/CastorObjects/interface/CastorCalibrations.h"
#include "CalibFormats/CastorObjects/interface/CastorCalibrationWidths.h"
#include "CalibCalorimetry/CastorCalib/interface/CastorDbASCIIIO.h"
#include "DataFormats/CastorReco/interface/CastorTower.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <HepMC/GenEvent.h>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Misc
#include <iostream>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
using namespace std;

//
// class declaration
//

#define PI 3.14159265359

double sectorToPhi(const int & sector) //not tested
{
  assert(1 <= sector && sector <= 16);

  double sec = double(sector) - 0.5;
  if (sec < 8)
    return sec/8.*PI;
  else
    return -PI + (sec-8)*PI/8;
}

int phiToSector(const double & phi)
{
double a=PI/8.;

int n = int(phi/a);
if(phi >= 0) n++;
else n = 16 + n;

return n;
}

class MyPhiAnalyzer : public edm::EDAnalyzer {
public:
  explicit MyPhiAnalyzer(const edm::ParameterSet&);
  ~MyPhiAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  TProfile *hEnergy;
  TProfile *hEnergy5;
  TGraph *hEnergy_res;
  TGraph *hEnergy_res5;
  TH1D *hEta_g_d;
  TH1D *hEta_g_d_N;
  TH1D *hEta_det;
  TH1D *hEta_det_N;
  TH1D *hPhi_det;
  TH1D *hPhi_det_N;
  TH1D *hPhi_det_1;
  TH1D *hPhi_det_2;
  TH1D *hPhi_det_3;
  TH1D *hPhi_det_5mod;
  TH1D *hPhi_gen;
  TH1D *hPhi_gen_N;
  TH2D *hPhi_gen_det;
  TH1D *hZ;
  TH1D *hZ_near;
  TH1D *hZ_far;
  TH1D *hE_det_vs_gen;
  TH1D *hE_det_vs_gen_helper;

  vector<TH1D*> hvZ;
  vector<TH1D*> hvZ_eta;
  vector<TH1D*> hvZ_eta_N;
  vector<TH1D*> hvZ_near;
  vector<TH1D*> hvZ_far;
  vector<TH1D*> hvPhi_det_eta;
  vector<TH1D*> hvPhi_det_eta_N;

  int nEvents;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MyPhiAnalyzer::MyPhiAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  TFileDirectory subDirHelper = fs->mkdir("helperHistos");
  TFileDirectory subDirEta = fs->mkdir("etaBinned");

  hEta_g_d     = fs->make<TH1D>("Efac_over_eta_gen","Eta boundaries of CASTOR; eta_{GEN}; E_{tot,gen}/E_{tot,det}",50,-6.8,-5);
  hEta_g_d_N   = subDirHelper.make<TH1D>("Efac_over_eta_gen_N","Eta boundaries of CASTOR; eta_{GEN}; N",50,-6.8,-5);
  hEta_det     = fs->make<TH1D>("E_over_eta_gen","Eta boundaries of CASTOR; eta_{GEN}; E_{tot,det}",50,-6.8,-5);
  hEta_det_N   = subDirHelper.make<TH1D>("E_over_eta_gen_N","Eta boundaries of CASTOR; eta_{GEN}; N",50,-6.8,-5);
  hPhi_det     = fs->make<TH1D>("phi_det","phi_detsector",16,0.5,16.5);
  hPhi_det_N   = subDirHelper.make<TH1D>("phi_det_N","phi_detsector",16,0.5,16.5);
  hPhi_det_1   = fs->make<TH1D>("phi_det_1","phi_detsector",16,0.5,16.5);
  hPhi_det_2   = fs->make<TH1D>("phi_det_2","phi_detsector",16,0.5,16.5);
  hPhi_det_3   = fs->make<TH1D>("phi_det_3","phi_detsector",16,0.5,16.5);
  hPhi_det_5mod     = fs->make<TH1D>("phi_det_5mod","phi_detsector",16,0.5,16.5);
  hPhi_gen     = fs->make<TH1D>("phi_gen","phi_gen;sector;#frac{1}{N} E /GeV",16,0.5,16.6);
  hPhi_gen_N   = subDirHelper.make<TH1D>("phi_gen_N","phi_gen;sector;N",16,0.5,16.6);
  hPhi_gen_det = fs->make<TH2D>("phi_gen_det","Hits in First Module (E>3GeV);#phi / degree;sector",16,-PI,PI,16,0.5,16.5);
  hZ           = fs->make<TH1D>("z","Profile;module;#frac{1}{N} #sumEnergy /GeV",14,0.5,14.5);
  hZ_near      = fs->make<TH1D>("z_near","Profile;module;#frac{1}{N} #sumEnergy /GeV",14,0.5,14.5);
  hZ_far       = fs->make<TH1D>("z_far","Profile;module;#frac{1}{N} #sumEnergy /GeV",14,0.5,14.5);
  hE_det_vs_gen = fs->make<TH1D>("E_det_vs_gen","Calo E over gen E;gen E in GeV; #sumEnergy /GeV",30,0,300);
  hE_det_vs_gen_helper = subDirHelper.make<TH1D>("E_det_vs_gen_helper","Calo E over gen E;gen E in GeV; #sumEnergy /GeV",30,0,300);
  hEnergy      = fs->make<TProfile>("energy","energy;energy (MC truth) in GeV;energy (reconstructed) in GeV",41,-5,405,"s");
  hEnergy5     = fs->make<TProfile>("energy5","energy;energy (MC truth) in GeV;energy (reconstructed) in GeV",41,-5,405,"s");
  hEnergy_res  = fs->make<TGraph>(41);
  hEnergy_res5 = fs->make<TGraph>(41);

  hEnergy_res->SetName("hEnergy_res");
  hEnergy_res5->SetName("hEnergy_res5");

  hEta_g_d->Sumw2();
  hEta_det->Sumw2();
  hPhi_gen->Sumw2();
  hPhi_det->Sumw2();
  hPhi_det_1->Sumw2();
  hPhi_det_2->Sumw2();
  hPhi_det_3->Sumw2();
  hPhi_det_5mod->Sumw2();
  hE_det_vs_gen->Sumw2();

  for (int i=1; i<=16; i++)
    {
      ostringstream name;
      name << "z_near" << i;
      hvZ_near.push_back(subDirHelper.make<TH1D>(name.str().c_str(),"Profile;#frac{1}{N} #sumEnergy /GeV",1000,-10,900));
      hvZ_near.back()->StatOverflows();
      name.str("");
      name << "z_far" << i;
      hvZ_far.push_back(subDirHelper.make<TH1D>(name.str().c_str(),"Profile;#frac{1}{N} #sumEnergy /GeV",1000,-10,900));
      hvZ_far.back()->StatOverflows();
      name.str("");
      name << "z" << i;
      hvZ.push_back(subDirHelper.make<TH1D>(name.str().c_str(),"Profile;#frac{1}{N} #sumEnergy /GeV",1000,-10,900));
      hvZ.back()->StatOverflows();
    }

  double eta;
  for (eta = -5.2; eta >= -6.6; eta-=0.1)
    {
      ostringstream name;
      name.str("");
      name << "z_eta_"; if(eta<0) name << "n"; name << fabs(eta);
      hvZ_eta.push_back(subDirEta.make<TH1D>(name.str().c_str(),";module no.; #sumEnergy /GeV",14,0.5,14.5));
      hvZ_eta.back()->Sumw2();

      name.str("");
      name << "z_eta_N_"; if(eta<0) name << "n"; name << fabs(eta);
      hvZ_eta_N.push_back(subDirEta.make<TH1D>(name.str().c_str(),";module no.; #sumEnergy /GeV",14,0.5,14.5));
 
      name.str("");
      name << "phi_det_eta_"; if(eta<0) name << "n"; name << fabs(eta);
      hvPhi_det_eta.push_back(subDirEta.make<TH1D>(name.str().c_str(),"sector; #sumEnergy /GeV",16,0.5,16.5));
      hvPhi_det_eta.back()->Sumw2();

      name.str("");
      name << "phi_det_eta_N_"; if(eta<0) name << "n"; name << fabs(eta);
      hvPhi_det_eta_N.push_back(subDirHelper.make<TH1D>(name.str().c_str(),";sector; N",16,0.5,16.5));
    }
}


MyPhiAnalyzer::~MyPhiAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyPhiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  ++nEvents; //event counter

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  ////CASTOR
  edm::Handle<CastorRecHitCollection> casrechits;
  try{ iEvent.getByLabel("castorreco",casrechits); }
  catch(...) { edm::LogWarning("CASTOR ") << " (*A*) Cannot get Castor RecHits " << std::endl; }

  //bool hasCastorHits = false;
  //double energyCastor = 0;

  if(casrechits.failedToGet()!=0 || !casrechits.isValid()) {

    edm::LogVerbatim("WARN") << " (*W*) Empty CastorRecHitCollection" << std::endl;
    return;

  }

  ////GENPARTICLES
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  try{ iEvent.getByLabel("genParticles", genParticleHandle); }
  catch(...) { edm::LogWarning("GEN ") << " Cannot get heavy ion GenParticles" << std::endl; }

  const reco::GenParticleCollection *genParticles = genParticleHandle.failedToGet()? 0 : &*genParticleHandle;

  if (!genParticles)
    {
      edm::LogVerbatim("WARN") << " (*W*) Empty CastorRecHitCollection" << std::endl;
      return;
    }



  vector<double> energyModuleMeanNear(14,0);
  vector<double> energyModuleMeanFar(14,0);
  vector<double> energySector(16,0);
  vector<double> energySector5Mod(16,0);
  vector<double> energySectorGen(16,0);
  double totalEInCastor = 0.;
  double totalEInCastor5Mod = 0.;
  double totalGenE = 0.;
  double primaryGenE = 0.;
  double particleGunEta = 1e99;
  int nFinalP = 0;

  for(size_t i1 = 0; i1 < casrechits->size(); ++i1)
    {

      const CastorRecHit & rh = (*casrechits)[i1];
      HcalCastorDetId castorid = rh.id();
      int secNb = castorid.sector();
      int modNb = castorid.module();

      double recHitE = rh.energy();

      totalEInCastor += recHitE;
      if (modNb <= 5)
        totalEInCastor5Mod += recHitE;
      //cout << i1 << " --- mod:" << modNb << " sec:" << secNb;
      //if (recHitE > 1) cout << " E: " << recHitE << endl;
      //else             cout << endl;
      if (secNb <= 4 || secNb >= 13)
        energyModuleMeanNear[modNb-1] += recHitE;
      else
        energyModuleMeanFar[modNb-1] += recHitE;

      energySector[secNb-1] += recHitE;
      if(modNb <= 5) energySector5Mod[secNb-1] += recHitE;

      if (modNb==1)
        hPhi_det_1->Fill(secNb,recHitE);
      if (modNb==2)
        hPhi_det_2->Fill(secNb,recHitE);
      if (modNb==3)
        hPhi_det_3->Fill(secNb,recHitE);
    }//rechits


  for(unsigned int j = 0; j < genParticles->size(); ++j)
    {

      const double energy = (*genParticles)[j].energy();
      const double eta = (*genParticles)[j].eta();
      const double phi = (*genParticles)[j].phi();

      if((*genParticles)[j].status() == 1)
	primaryGenE += energy;

      if ((*genParticles)[j].status() != 1)
	continue;

      if(eta<0) //remove antiparticle from particlegun
	{
	  particleGunEta = eta;
	  nFinalP++;
	}
      if (eta>-5.2 || eta<-6.6)
	continue;
      //only CASTOR

      totalGenE += energy;
      int genSec = phiToSector(phi);
      energySectorGen[genSec-1] += energy;
      
    }//genparticles

  //AFTER LOOPS
  //Getting histogram to corresponding particle gun eta
  int jeta=int((-5.2-particleGunEta)/0.1);
  TH1D* etaZHist = 0;
  TH1D* etaZHistN = 0;
  if (nFinalP != 1)
    { //remove
      cout << "more than 1 final particle " << particleGunEta << " " << jeta << endl;
    }
  if(nFinalP == 1 && jeta >= 0 && jeta < int(hvZ_eta.size()))
    {
      etaZHist = hvZ_eta[jeta];
      etaZHistN = hvZ_eta_N[jeta];
    }

  for (int i=1; i<=int(energySector.size()); i++)
    {
      hPhi_det->Fill(i,energySector[i-1]);
      hPhi_det_N->Fill(i);

      hPhi_det_5mod->Fill(i,energySector5Mod[i-1]);

      hPhi_gen->Fill(i,energySectorGen[i-1]);
      hPhi_gen_N->Fill(i);

      for (int j=1; j<=int(energySectorGen.size()); j++)
	{
	  if(energySectorGen[i-1] > 3 && energySector[j-1] > 3)
	    hPhi_gen_det->Fill(sectorToPhi(i),j);
	}

      if(nFinalP == 1 && jeta >= 0 && jeta < int(hvPhi_det_eta.size()))
	{
	  hvPhi_det_eta[jeta]->Fill(i,energySector[i-1]);
	  hvPhi_det_eta_N[jeta]->Fill(i);
	}
    }

  for (int i=1; i<=int(energyModuleMeanNear.size()); i++)
    { //fill every module
      hvZ_near[i-1]->Fill(energyModuleMeanNear[i-1]); //helper hists
      hvZ_far[i-1]->Fill(energyModuleMeanFar[i-1]);
      hvZ[i-1]->Fill(energyModuleMeanFar[i-1]+energyModuleMeanNear[i-1]);

      if (etaZHist && etaZHistN)
	{
	  etaZHist->Fill(i,energyModuleMeanNear[i-1]);
	  etaZHistN->Fill(i);
	}
    }

  if (nFinalP == 1)
    {
      hEta_g_d->Fill(particleGunEta,totalGenE/totalEInCastor);
      hEta_det->Fill(particleGunEta,totalEInCastor);
      hEta_g_d_N->Fill(particleGunEta);
      hEta_det_N->Fill(particleGunEta);
    }
  hE_det_vs_gen->Fill(totalGenE,totalEInCastor);
  hE_det_vs_gen_helper->Fill(totalGenE);
  hEnergy->Fill(primaryGenE,totalEInCastor);
  hEnergy5->Fill(primaryGenE,totalEInCastor5Mod);
}


// ------------ method called once each job just before starting event loop  ------------
void
MyPhiAnalyzer::beginJob()
{
  nEvents=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
MyPhiAnalyzer::endJob()
{
  //Fill near and far graphs from single histograms for each module
  for (int i=1; i<=16; i++)
    {
      hZ_far->SetBinContent(i,hvZ_far[i-1]->GetMean());
      hZ_far->SetBinError(i,hvZ_far[i-1]->GetRMS()/sqrt(hvZ_far[i-1]->GetEntries()));

      hZ_near->SetBinContent(i,hvZ_near[i-1]->GetMean());
      hZ_near->SetBinError(i,hvZ_near[i-1]->GetRMS()/sqrt(hvZ_near[i-1]->GetEntries()));

      hZ->SetBinContent(i,hvZ[i-1]->GetMean());
      hZ->SetBinError(i,hvZ[i-1]->GetRMS()/sqrt(hvZ[i-1]->GetEntries()));
    }

  //Scale total energy histos to get mean
  for (int i=1; i<=hE_det_vs_gen->GetNbinsX(); i++)
    {
      if(hE_det_vs_gen_helper->GetBinContent(i))
        {
          hE_det_vs_gen->SetBinContent(i,hE_det_vs_gen->GetBinContent(i)/ hE_det_vs_gen_helper->GetBinContent(i));
          hE_det_vs_gen->SetBinError  (i,hE_det_vs_gen->GetBinError  (i)/(hE_det_vs_gen_helper->GetBinContent(i)));
        }
    }

  //Phi histos
  for (int jeta = 0; jeta < int(hvZ_eta.size()); jeta++) //eta bins
    {
      for(int ibin=1; ibin<=hvZ_eta_N[jeta]->GetNbinsX(); ibin++) //zmodule bins
	{      
	  if(hvZ_eta_N[jeta]->GetBinContent(ibin)) hvZ_eta[jeta]->SetBinContent(ibin, hvZ_eta[jeta]->GetBinContent(ibin) / hvZ_eta_N[jeta]->GetBinContent(ibin));
	  if(hvZ_eta_N[jeta]->GetBinContent(ibin)) hvZ_eta[jeta]->SetBinError(ibin, hvZ_eta[jeta]->GetBinError(ibin) / hvZ_eta_N[jeta]->GetBinContent(ibin));
	}
    }

  for(int i=1; i<=hEta_g_d_N->GetNbinsX(); i++)
    {
      if(hEta_g_d_N->GetBinContent(i)) hEta_g_d->SetBinContent(i, hEta_g_d->GetBinContent(i) / hEta_g_d_N->GetBinContent(i));
      if(hEta_det_N->GetBinContent(i)) hEta_det->SetBinContent(i, hEta_det->GetBinContent(i) / hEta_det_N->GetBinContent(i));

      if(hEta_g_d_N->GetBinContent(i)) hEta_g_d->SetBinError(i, hEta_g_d->GetBinError(i) / hEta_g_d_N->GetBinContent(i));
      if(hEta_det_N->GetBinContent(i)) hEta_det->SetBinError(i, hEta_det->GetBinError(i) / hEta_det_N->GetBinContent(i));
    }

  for(int i=1; i<=hPhi_gen_N->GetNbinsX(); i++)
    {
      if(hPhi_gen_N->GetBinContent(i)) hPhi_gen->SetBinContent(i, hPhi_gen->GetBinContent(i) / hPhi_gen_N->GetBinContent(i));
      if(hPhi_det_N->GetBinContent(i)) hPhi_det->SetBinContent(i, hPhi_det->GetBinContent(i) / hPhi_det_N->GetBinContent(i));
      if(hPhi_det_N->GetBinContent(i)) hPhi_det_5mod->SetBinContent(i, hPhi_det_5mod->GetBinContent(i) / hPhi_det_N->GetBinContent(i));

      if(hPhi_gen_N->GetBinContent(i)) hPhi_gen->SetBinError(i, hPhi_gen->GetBinError(i) / hPhi_gen_N->GetBinContent(i));
      if(hPhi_det_N->GetBinContent(i)) hPhi_det->SetBinError(i, hPhi_det->GetBinError(i) / hPhi_det_N->GetBinContent(i));
      if(hPhi_det_N->GetBinContent(i)) hPhi_det_5mod->SetBinError(i, hPhi_det_5mod->GetBinError(i) / hPhi_det_N->GetBinContent(i));
 
      for (int jeta = 0; jeta < int(hvPhi_det_eta.size()); jeta++) //eta bins
	{
	  if(hvPhi_det_eta_N[jeta]->GetBinContent(i)) hvPhi_det_eta[jeta]->SetBinContent(i, hvPhi_det_eta[jeta]->GetBinContent(i) / hvPhi_det_eta_N[jeta]->GetBinContent(i));
	  if(hvPhi_det_eta_N[jeta]->GetBinContent(i)) hvPhi_det_eta[jeta]->SetBinError(i, hvPhi_det_eta[jeta]->GetBinError(i) / hvPhi_det_eta_N[jeta]->GetBinContent(i));
	}
    }

  //Energy resolution
  edm::Service<TFileService> fs;
  TGraph* curRes = hEnergy_res;
  TProfile* cur = hEnergy;
  for(int i=0; i<2; i++)
  {
    for (int j=1; j<cur->GetNbinsX()+1; j++)
      {
        if(cur->GetBinContent(j)<=0)
          continue;
        curRes->SetPoint(j,cur->GetBinCenter(j),cur->GetBinError(j)/cur->GetBinCenter(j));
      }
    curRes = hEnergy_res5;
    cur = hEnergy5;
  }

}

// ------------ method called when starting to processes a run  ------------
void
MyPhiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
MyPhiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
MyPhiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MyPhiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyPhiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyPhiAnalyzer);
