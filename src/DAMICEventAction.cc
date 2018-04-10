
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Muon Induced Background
//
//      For information related to this code contact: M. Settimo
//    e-mail: settimo@subatech.in2p3.fr
// --------------------------------------------------------------
// Comments
//
//                   Underground Advanced -  Modified by M.Settimo 2017
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// History/Additions:
// 16 Jan 2002  Added analysis
//
//
// EventAction program
// --------------------------------------------------------------

#include "DAMICEventAction.hh"

// pass parameters for messengers:
#include "DAMICRunAction.hh"
#include "DAMICPrimaryGeneratorAction.hh"

#include "DAMICEventActionMessenger.hh"
#include "DAMICAnalysisManager.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include <fstream>
#include <iomanip>

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DAMICEventAction::DAMICEventAction() 
  : runAct(0),genAction(0),hitsfile(0)
{

  // create messenger
  eventMessenger = new DAMICEventActionMessenger(this);

  // defaults for messenger
  drawColsFlag = "standard";
  drawTrksFlag = "all";
  drawHitsFlag = 1;
  saveHitsFlag = 1;
  
  printModulo = 1;

  // hits collections
  ccdCollID = -1;


  energy_pri=0;
  primX=primY=primZ=0;
  primType="";
  primLV="";
  seeds=NULL;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DAMICEventAction::~DAMICEventAction() {

  if (hitsfile)
    {
      hitsfile->close();
      delete hitsfile;
    }
 
  delete eventMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DAMICEventAction::BeginOfEventAction(const G4Event* evt) 
{

  //thread-local run action
  if (!runAct) 
    runAct = 
      dynamic_cast<const DAMICRunAction*>
      (G4RunManager::GetRunManager()->GetUserRunAction());
  
  if (!genAction)
    genAction = dynamic_cast<const DAMICPrimaryGeneratorAction*>
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());


  // grab seeds
  seeds = genAction->GetEventSeeds();

  // grab energy of primary
  energy_pri = genAction->GetEnergyPrimary();
  //SetEventEnergy(energy_pri);
  SetPrimaryEvent(evt);
  
  event_id = evt->GetEventID();
  primType = evt->GetPrimaryVertex(0)->GetPrimary()->GetParticleDefinition()->GetParticleType();
  primX = evt->GetPrimaryVertex()->GetX0(); //GetPosition()->x();
  primY = evt->GetPrimaryVertex()->GetY0(); //GetPosition()->x();
  primZ = evt->GetPrimaryVertex()->GetZ0(); //GetPosition()->x();

 
  // print this information event by event (modulo n)  	
  if (event_id%printModulo == 0)
    {
/*      G4cout << "\n---> Begin of event: " << event_id << G4endl;
      G4cout << "       Primary Energy: " << G4BestUnit(energy_pri,"Energy")  << G4endl;
      G4cout << "       Primary type: " << primType << " @ " << primX << " " << primY << " " << primZ   << G4endl;
*/    
      //      HepRandom::showEngineStatus(); 
    }


  // get ID for scintillator hits collection
  
  if (ccdCollID==-1) {
    // G4SDManager *SDman = G4SDManager::GetSDMpointer();
    // ccdCollID = SDman->GetCollectionID("ccdCollection");
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DAMICEventAction::EndOfEventAction(const G4Event* evt) {
  G4cout << " in endOfeventAction " << G4endl; 
  // check that both hits collections have been defined
  // if(ccdCollID<0) { G4cout << " ccdColID < 0. break" << G4endl; return;}

  // G4AnalysisManager* man = G4AnalysisManager::Instance();

  // // address hits collections
  // DAMICSiHitsCollection* SHC = NULL;
  // G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  // if(HCE) {
  //   SHC = (DAMICSiHitsCollection*)(HCE->GetHC(ccdCollID));
  // }

  // // event summary
  // totEnergy         = 0.;
  // totEnergyGammas   = 0.;
  // totEnergyNeutrons = 0.;
  // firstParticleE    = 0.;
  // particleEnergy    = 0.;
  // firstHitTime      = 0.;
  // firstTrackID      = 0;
  // firstStepID       = 0;
  // firstParticleName = "";
  // particleName      = "";


  // // particle flags
  // gamma_ev          = false;
  // neutron_ev        = false;
  // positron_ev       = false;
  // electron_ev       = false;
  // proton_ev         = false;
  // other_ev          = false;
  // start_gamma       = false;
  // start_neutron     = false;

  // // write out event summary
  // if(saveHitsFlag)  
  //   writeSiHitsToFile();
    
  // // draw trajectories
  // if(drawColsFlag=="standard" && drawTrksFlag!="none")
  //   drawTracks(evt);

  // // print this event by event (modulo n)  	
  // if (event_id%printModulo == 0) 
     G4cout << "\n---> End of event: " << event_id << G4endl << G4endl;	

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DAMICEventAction::writeSiHitsToFile() 
{
//   G4cout << " test " << G4endl << G4endl; 
//   G4String filename="hits.out";
//   if (runAct)
//     filename=runAct->GetsavehitsFile();

 
  
//   //First time it is inkoved
//   if (!hitsfile)
//     {
//       //check that we are in a worker: returns -1 in a master and -2 in sequential
//       //one file per thread is produced ending with ".N", with N= thread number
//       if (G4Threading::G4GetThreadId() >= 0)
// 	{
// 	  std::stringstream sss;
// 	  sss << filename.c_str() << "." << G4Threading::G4GetThreadId();	 
// 	  filename = sss.str();
// 	  //G4cout << "Filename is: " << filename << G4endl;
// 	}
      
//       hitsfile = new std::ofstream;
//       hitsfile->open(filename);
//       (*hitsfile) <<"Evt     Eprim   Etot    LXe     LXeTime   ddd  Time Seed1           Seed2           First   Flags" 
// 	       << G4endl;
//       (*hitsfile) <<"#       MeV     MeV     hits    ns      hits    ns                                      hit"
// 	       << G4endl
// 	       << G4endl;
//     }//end if !hist

// G4AnalysisManager* man = G4AnalysisManager::Instance();
//   if(S_hits) {

//     if(hitsfile->is_open()) {


//       (*hitsfile) << std::setiosflags(std::ios::fixed)
// 		  << std::setprecision(4)
// 		  << std::setiosflags(std::ios::left)
// 		  << std::setw(6)
// 		  << event_id << "\t"
// 		  << energy_pri/MeV << "\t" 
// 		  << totEnergy/MeV << "\t"
// 		  << S_hits  << "\t"
// 		  << std::setiosflags(std::ios::scientific) 
// 		  << std::setprecision(2)
// 		  << firstHitTime/nanosecond << "\t"
// 		  << P_hits << "\t"
// 		  << std::setiosflags(std::ios::fixed) 
// 		  << std::setprecision(4)
// 		  << *seeds     << "\t"
// 		  << *(seeds+1) << "\t"
// 		  << firstParticleName << "\t"
// 		  << (gamma_ev    ? "gamma " : "") 
// 		  << (neutron_ev  ? "neutron " : "") 
// 		  << (positron_ev ? "positron " : "") 
// 		  << (electron_ev ? "electron " : "") 
// 		  << (other_ev    ? "other " : "") 
// 		  << G4endl;

//       if (event_id%printModulo == 0)
// 	G4cout << "     Event summary in file " << filename << G4endl;  
//     }//end hitsfile open

//     // G4AnalysisManager* man = G4AnalysisManager::Instance();
//     G4int firstparticleIndex = 0;
//     if(firstParticleName == "gamma") firstparticleIndex = 1;
//     else if (firstParticleName == "neutron") firstparticleIndex = 2;
//     else if(firstParticleName == "electron") firstparticleIndex = 3;
//     else if(firstParticleName == "positron") firstparticleIndex = 4;
//     else{
//       firstparticleIndex = 5;
//       man->FillH1(3,totEnergy);
//     }

//     man->FillH1(4,P_hits,10); //weight
//     man->FillH1(5,P_hits);

//     man->FillH1(1,energy_pri/keV);
//     man->FillH1(2,totEnergy/keV);

//     long seed1 = *seeds;
//     long seed2 = *(seeds+1);

    
//     //Fill ntuple #2
//     man->FillNtupleDColumn(2,0,event_id);
//     man->FillNtupleDColumn(2,1,energy_pri/keV);
//     man->FillNtupleDColumn(2,2,totEnergy);
//     man->FillNtupleDColumn(2,3,S_hits);
//     man->FillNtupleDColumn(2,4,firstHitTime);
//     man->FillNtupleDColumn(2,5,P_hits);
//     man->FillNtupleDColumn(2,7,firstparticleIndex);
//     man->FillNtupleDColumn(2,8,firstParticleE);
//     man->FillNtupleDColumn(2,9,gamma_ev);
//     man->FillNtupleDColumn(2,10,neutron_ev);
//     man->FillNtupleDColumn(2,11,positron_ev);
//     man->FillNtupleDColumn(2,12,electron_ev);
//     man->FillNtupleDColumn(2,13,other_ev);
//     man->FillNtupleDColumn(2,14,seed1);
//     man->FillNtupleDColumn(2,15,seed2);
//     man->AddNtupleRow(2);
    
//   }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DAMICEventAction::drawTracks(const G4Event* evt) {

  if(G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");    
    G4TrajectoryContainer* trajContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;

    if(trajContainer) n_trajectories = trajContainer->entries();
    for (G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory*)(*trajContainer)[i];
      if (drawTrksFlag == "all") 
	trj->DrawTrajectory();
      else if ((drawTrksFlag == "charged") && (trj->GetCharge() != 0.))
	trj->DrawTrajectory();
      else if ((drawTrksFlag == "noscint") 
	       && (trj->GetParticleName() != "opticalphoton"))
	trj->DrawTrajectory();
    }
    
    // G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");    
  } 

}
