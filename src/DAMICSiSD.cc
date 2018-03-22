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
// SiSD (Si sensitive detector definition) program
// --------------------------------------------------------------

#include "DAMICSiSD.hh"

#include "DAMICSiHit.hh"
#include "DAMICDetectorConstruction.hh"
#include "DAMICAnalysisManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4Ions.hh"
#include "G4ios.hh"
#include <vector>
#include "G4OpticalPhoton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DAMICSiSD::DAMICSiSD(G4String name) 
  :G4VSensitiveDetector(name)
{
  G4String HCname="ccdCollection";
  collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DAMICSiSD::~DAMICSiSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICSiSD::Initialize(G4HCofThisEvent*)
{
    ccdCollection = new DAMICSiHitsCollection(SensitiveDetectorName,collectionName[0]);

    HitID = -1;

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool DAMICSiSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  //need to know if this is an optical photon and exclude it:
  //  if(aStep->GetTrack()->GetDefinition()  == G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  //  G4cout << " Process Hits for " <<  aStep->GetTrack()->GetTrackID() << G4endl;
  
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int stepID =  aStep->GetTrack()->GetCurrentStepNumber();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
  G4int motherID = aStep->GetTrack()->GetParentID();
  G4String particleName = particle->GetParticleName();
  G4int particlePDG = particle->GetPDGEncoding();

  G4StepPoint* point1 = aStep->GetPreStepPoint();
  G4StepPoint* point2 = aStep->GetPostStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4String ccdName = touch1->GetVolume()->GetName();
 
  if ( aStep->GetTrack()->GetCreatorProcess()!=0) {
  
     if (particlePDG==11 &&  aStep->GetTrack()->GetCreatorProcess()->GetProcessName()=="eIoni" &&  !aStep->IsFirstStepInVolume()) {
    //  if (particlePDG==11 &&  point1->GetProcessDefinedStep()->GetProcessName()=="eIoni" &&  !aStep->IsFirstStepInVolume()) {
   
    //I want to stop all e^- that have been produced by a ionization process in the CCD.
    //So I check if the creator process is a ionization and I have to ensure that this is not the 1st hit of the particle in SD
     // (because it creator would be somewhere else)
     aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
     //  aStep->UpdateTrack();
    
     // if (aStep->GetSecondary()!=0) aStep->DeleteSecondaryVector ();
    const  G4TrackVector *secondTr = aStep->GetSecondary();
    //  G4cout << " still status "<< aStep->GetTrack()->GetTrackStatus() << " Nsecondary: " << secondTr->size() << G4endl;
    
    for (int h=0; h<(int)secondTr->size(); ++h) {
         if (secondTr->at(h)->GetDefinition()->GetParticleName() == "e-")
               secondTr->at(h)->SetTrackStatus(fStopAndKill);
    }
    
    aStep->UpdateTrack(); 
   
   return false;
   }     
  }

  // if (particlePDG==11 && !aStep->IsFirstStepInVolume() && aStep->GetTrack()->GetCurrentStepNumber()>1) 
  //   return false;


  
  
  // --- get volume + mother  + proc  info ---
  
  //  G4int copyNumber = touch1->GetCopyNumber();
  // const G4VProcess* aProcess = point2->GetProcessDefinedStep();
  G4String processName  = point2->GetProcessDefinedStep()->GetProcessName();

  G4String creatorProcess = "";
  if (aStep->GetTrack()->GetCreatorProcess()!=0)
    creatorProcess = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
  else
     creatorProcess = "None";
  
  G4ProcessType processType =  point2->GetProcessDefinedStep()->GetProcessType();
    
  G4double stepl = 0.;

  stepl = aStep->GetStepLength();
  
  if ((edep==0.)&&(stepl==0.)) return false;      
  
  // fill in hit
  DAMICSiHit* newHit = new DAMICSiHit();
  newHit->SetTrackID(trackID);
  newHit->SetStepID(stepID);
  newHit->SetEdep(edep);
  newHit->SetStepLength(stepl);
  newHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
  newHit->SetGlobalTime(aStep->GetPreStepPoint()->GetGlobalTime());
  newHit->SetLocalTime(aStep->GetPreStepPoint()->GetLocalTime());
  // newHit->SetParticle(particle);
  newHit->SetParticleName(particleName);
  newHit->SetParticlePDG(particlePDG);
  newHit->SetParticleEnergy(aStep->GetPreStepPoint()->GetKineticEnergy());
  newHit->SetVolumeName(ccdName);
  newHit->SetProcessName(processName);
  newHit->SetCreatorProcess(processName);
  newHit->SetProcessType(processType);
  newHit->SetMotherID(motherID);
  
  HitID = ccdCollection->insert(newHit);
  
  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICSiSD::EndOfEvent(G4HCofThisEvent* HCE)
{
 G4AnalysisManager* man = G4AnalysisManager::Instance();
  G4String HCname = collectionName[0];
  
  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID,ccdCollection);

  G4int nHits = ccdCollection->entries();
  //  G4cout << G4endl << HCname << " " << nHits << G4endl;
  G4double firstParticleE;
  G4String firstParticleName;
  G4int firstTrackID;
  G4int firstStepID;
  G4double firstHitTime;
  G4double totEnergy=0;
  for (G4int i=0; i<nHits; i++) {
    
    if(i==0) {
      firstTrackID = (*ccdCollection)[0]->GetTrackID();
      firstStepID = (*ccdCollection)[0]->GetStepID();
      firstParticleName = (*ccdCollection)[0]->GetParticleName();
      firstHitTime   = (*ccdCollection)[0]->GetGlobalTime();
      firstParticleE = (*ccdCollection)[0]->GetParticleEnergy();
    }
    
    G4int ptrackID = (*ccdCollection)[i]->GetTrackID();
    G4int pstepID =  (*ccdCollection)[i]->GetStepID();
    G4double hitEnergy         = (*ccdCollection)[i]->GetEdep();
    totEnergy      += hitEnergy;
    G4double x =  (*ccdCollection)[i]->GetPos().x();
    G4double y =  (*ccdCollection)[i]->GetPos().y();
    G4double z =  (*ccdCollection)[i]->GetPos().z();
    
    G4String particleName      = (*ccdCollection)[i]->GetParticleName();
    G4int particlePDG          = (*ccdCollection)[i]->GetParticlePDG();
    G4double particleEnergy    = (*ccdCollection)[i]->GetParticleEnergy();
    G4double GlobalTime        = (*ccdCollection)[i]->GetGlobalTime();
    G4double LocalTime         = (*ccdCollection)[i]->GetGlobalTime();
    G4String volName           = (*ccdCollection)[i]->GetVolumeName();
    G4String proName           = (*ccdCollection)[i]->GetProcessName();
    G4String creatorName        = (*ccdCollection)[i]->GetCreatorProcess();
    G4int    proType           = (int)(*ccdCollection)[i]->GetProcessType();
    G4double stepL             = (*ccdCollection)[i]->GetStepLength();
    G4int    motherID          = (int)(*ccdCollection)[i]->GetMotherID();
    
    //Fill ntuple #3
    // man->FillNtupleDColumn(1,0,firstHitName);
    man->FillNtupleIColumn(1,0,ptrackID);
    man->FillNtupleIColumn(1,1,pstepID);
    man->FillNtupleIColumn(1,2,nHits);
    man->FillNtupleDColumn(1,3,hitEnergy/CLHEP::keV);
    man->FillNtupleDColumn(1,4,totEnergy/CLHEP::keV);
    man->FillNtupleDColumn(1,5,x/CLHEP::cm);
    man->FillNtupleDColumn(1,6,y/CLHEP::cm);
    man->FillNtupleDColumn(1,7,z/CLHEP::cm);
    man->FillNtupleDColumn(1,8,stepL/CLHEP::cm);
    man->FillNtupleDColumn(1,9,GlobalTime/CLHEP::s);
    man->FillNtupleDColumn(1,10,LocalTime/CLHEP::s);
    man->FillNtupleDColumn(1,11,particleEnergy/CLHEP::keV);
    man->FillNtupleSColumn(1,12,particleName);
    man->FillNtupleIColumn(1,13,particlePDG);
    man->FillNtupleSColumn(1,14,volName);
    man->FillNtupleSColumn(1,15,proName);
    man->FillNtupleIColumn(1,16,proType);
    man->FillNtupleSColumn(1,17,creatorName);
    man->FillNtupleIColumn(1,18,motherID);
    //man->FillNtupleDColumn(1,12,grandmotherID);
    man->AddNtupleRow(1);
  }
  
  if (verboseLevel>=1)
    G4cout << "     Si collection: " <<  nHits << " hits" << G4endl;
  if (verboseLevel>=2)
    ccdCollection->PrintAllHits();
  

       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICSiSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICSiSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICSiSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

