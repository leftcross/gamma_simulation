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
// History:
// 21 Feb 2002 AH: Added Analysis
//
// SteppingAction program
// --------------------------------------------------------------

#include "DAMICSteppingAction.hh"
#include "DAMICSteppingActionMessenger.hh"

#include "DAMICEventAction.hh"
#include "DAMICAnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh" 
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DAMICSteppingAction::DAMICSteppingAction()
  : evtAction(0)  {

  steppingMessenger = new DAMICSteppingActionMessenger(this);

  // defaults for messenger
  colourNeutronFlag      = "magenta";
  colourGammaFlag        = "cyan";
  colourOpticalFlag      = "white";
  colourChargedPlusFlag  = "red";
  colourChargedMinusFlag = "blue";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DAMICSteppingAction::~DAMICSteppingAction()
{

  delete steppingMessenger;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICSteppingAction::UserSteppingAction(const G4Step* fStep)
{
  if (!evtAction)
    evtAction = 
      dynamic_cast<const DAMICEventAction*>
      (G4RunManager::GetRunManager()->GetUserEventAction());

  
  const G4Event *evt1 = evtAction->GetPrimaryEvent();
  G4int event_id = evt1->GetEventID();
  
  G4String primType = evt1->GetPrimaryVertex(0)->GetPrimary()->GetParticleDefinition()->GetParticleType();
  G4double primX = evt1->GetPrimaryVertex()->GetX0(); 
  G4double primY = evt1->GetPrimaryVertex()->GetY0(); 
  G4double primZ = evt1->GetPrimaryVertex()->GetZ0(); 
  
  G4double energy_pri = evtAction->GetEventEnergy();
 
// G4int StepNo=fStep->GetTrack()->GetCurrentStepNumber(); 
 
  G4int stepNsec = 0;
  if (fStep->GetSecondary()!=0)
    stepNsec  = fStep->GetSecondary()->size();
  else stepNsec=-1;


  G4String CreatorProcess = "";
  if (fStep->GetTrack()->GetCreatorProcess()!=0)
     CreatorProcess = fStep->GetTrack()->GetCreatorProcess()->GetProcessName();
  else CreatorProcess = "None";

  G4String StepProcess = "";
  if (fStep->GetPreStepPoint()->GetProcessDefinedStep()!=0) 
    StepProcess = fStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
  else StepProcess= "N.D.";

  G4String Volume = "";
  G4String Material = "";
  if (fStep->GetPreStepPoint()) {
    Material = fStep->GetPreStepPoint()->GetMaterial()->GetName();
    Volume = fStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
  }
  
  G4String NextMaterial = "";  
  G4String NextVolume = "";
  if (!fStep->GetTrack()->GetNextVolume()) return;  
  else NextVolume = fStep->GetTrack()->GetNextVolume()->GetName();
 // std::cout<<fStep->GetTrack()->GetCurrentStepNumber()<<std::endl;
   //G4cout<<" Volume is "<<Volume<<std::endl; 
// G4cout<<"Present volume is "<<Volume<<" Next volume is "<<NextVolume<<std::endl;

  if (fStep->GetPostStepPoint()!=0) 
       NextMaterial = fStep->GetPostStepPoint()->GetMaterial()->GetName();
  //    NextMaterial = fStep->GetTrack()->GetMaterial()->GetName(); 

  G4int trackID = fStep->GetTrack()->GetTrackID();
  G4int StepNo = fStep->GetTrack()->GetCurrentStepNumber();
  G4int partPDG = fStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4String partType = fStep->GetTrack()->GetDefinition()->GetParticleName();
  G4double partEtot = fStep->GetTrack()->GetTotalEnergy();
  G4double partEkin = fStep->GetTrack()->GetKineticEnergy();

  //std::cout<<"Particle was "<<partType<<" with PDG Code "<<partPDG<<" Steprocess "<<StepProcess<<std::endl;

  G4double partX=fStep->GetTrack()->GetPosition().x();
  G4double partY=fStep->GetTrack()->GetPosition().y();
  G4double partZ=fStep->GetTrack()->GetPosition().z();

  /* 
  G4double partTau= fStep->GetTrack()->GetDefinition()->GetPDGLifeTime();
  G4double partGlobalTime = fStep->GetTrack()->GetGlobalTime();
  G4double partProperTime = fStep->GetTrack()->GetProperTime();
  */
 
  G4double Edep = fStep->GetTotalEnergyDeposit () ;
  G4double Edep_nonIon = fStep->GetNonIonizingEnergyDeposit ();
  G4double stepLength = fStep->GetStepLength ();
  G4int motherId = fStep->GetTrack()->GetParentID();
  G4double partTrackLength = fStep->GetTrack()->GetTrackLength();
  
 G4int StoreNeutron=0;
 G4double offset=20;
 G4double labHeight=800;
 G4double leadboxh=100.18;
 G4double leaboxlb=200;
 G4double leadboxposz=-(labHeight/2.0 - leadboxh/2.)+offset;

 if ((fStep->GetPreStepPoint()->GetProcessDefinedStep()==0)&&(StepNo==1))
	{
		if(partPDG==2112)
		{ 
	//	StoreNeutron=1;
 		if(((partZ/10.0)<(leadboxposz+leadboxh/4.0))&&((partZ/10.0)>(leadboxposz-leadboxh/4.0)))
			{
		StoreNeutron=1;
		//std::cout<<" Storing neutron for position "<<partZ<<" Track ID "<<trackID<<std::endl;
			}
		//else
		//std::cout<<" Skipping neutron for position "<<partZ<<std::endl;
		}
	}

//if(  (Volume == "extBonnerSphereLV") && ((NextVolume == "emptyLeadBoxPV")||(NextVolume =="extLeadBoxPV")) ) {
  if( NextVolume == "CCDPV") {  
   // if(StoreNeutron==1) {  
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn(0,0,event_id);
    man->FillNtupleDColumn(0,1,energy_pri/CLHEP::keV);
    man->FillNtupleDColumn(0,2,primX/CLHEP::cm); //3 
    man->FillNtupleDColumn(0,3,primY/CLHEP::cm); //4 
    man->FillNtupleDColumn(0,4,primZ/CLHEP::cm); //5
    man->FillNtupleIColumn(0, 5, trackID);   //7
    man->FillNtupleIColumn(0, 6, StepNo);   //7 
    man->FillNtupleIColumn(0, 7, partPDG);   //8
 
    man->FillNtupleDColumn(0, 8, partEtot/CLHEP::keV); //10 
    man->FillNtupleDColumn(0, 9, partEkin/CLHEP::keV); //11
    man->FillNtupleDColumn(0, 10, partX/CLHEP::cm);     //13
    man->FillNtupleDColumn(0, 11, partY/CLHEP::cm);     //14 
    man->FillNtupleDColumn(0, 12, partZ/CLHEP::cm);     //15 
    // man->FillNtupleDColumn(0, 16, partTau/CLHEP::s);   //12
    // man->FillNtupleDColumn(0, 17, partGlobalTime/CLHEP::s);   //12
    // man->FillNtupleDColumn(0, 18, partProperTime/CLHEP::s);   //12
    //  man->FillNtupleDColumn(0, 13, partTrackLength/CLHEP::cm); //16
    man->FillNtupleSColumn(0, 13, CreatorProcess);
    
    man->FillNtupleSColumn(0, 14, StepProcess);
    man->FillNtupleDColumn(0, 15, Edep/CLHEP::keV);            //18
    man->FillNtupleDColumn(0, 16, Edep_nonIon/CLHEP::keV);            //18
    man->FillNtupleDColumn(0, 17, stepLength/CLHEP::cm);            //18
    
    man->FillNtupleSColumn(0, 18, Volume);            //18
    man->FillNtupleSColumn(0, 19, NextVolume);        //19
    man->FillNtupleSColumn(0, 20, Material);          //20
    man->FillNtupleSColumn(0, 21, NextMaterial);      //21
    man->FillNtupleIColumn(0, 22, stepNsec); 
    man->FillNtupleIColumn(0, 23, motherId); //16
    man->AddNtupleRow(0);

    //fStep->GetTrack()->SetTrackStatus(fStopAndKill);
    //fStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

}
  // check what is to be drawn from EventAction/EventActionMessenger
  G4String drawColsFlag = evtAction->GetDrawColsFlag();
  G4String drawTrksFlag = evtAction->GetDrawTrksFlag();

  // draw by step (here) instead of by event (event action)
  if (drawColsFlag=="custom" && drawTrksFlag!="none") {

    // check that VisManager exists
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager) {

      // particle colour in a string
      G4String name = fStep->GetTrack()->GetDefinition()->GetParticleName();
      G4String strColour;
      if(name=="neutron") {
	if(drawTrksFlag=="charged") return;
	strColour = colourNeutronFlag;
      } else if (name=="gamma") {
	if(drawTrksFlag=="charged") return;
	strColour = colourGammaFlag;
      } else if (name=="opticalphoton") {
	if(drawTrksFlag!="all") return;
	strColour = colourOpticalFlag;
      }
      else if (name=="alpha" || name=="e+")
	strColour = colourChargedPlusFlag;
      else
	strColour = colourChargedMinusFlag;

      // convert string to G4Colour
      G4Colour colour;
      if     (strColour=="white")    colour=G4Colour(1.0, 1.0, 1.0);
      else if(strColour=="grey" )    colour=G4Colour(0.5, 0.5, 0.5);
      else if(strColour=="lgrey")    colour=G4Colour(.75, .75, .75);
      else if(strColour=="black")    colour=G4Colour(0.0, 0.0, 0.0);
      else if(strColour=="red")      colour=G4Colour(1.0, 0.0, 0.0);
      else if(strColour=="green")    colour=G4Colour(0.0, 1.0, 0.0);
      else if(strColour=="blue")     colour=G4Colour(0.0, 0.0, 1.0);
      else if(strColour=="cyan")     colour=G4Colour(0.0, 1.0, 1.0);
      else if(strColour=="magenta")  colour=G4Colour(1.0, 0.0, 1.0);
      else if(strColour=="yellow")   colour=G4Colour(1.0, 1.0, 0.0);
      else if(strColour=="lgreen")   colour=G4Colour(0.0, .75, 0.0);
      else if(strColour=="lblue")    colour=G4Colour(0.0, 0.0, .75);
      else                           colour=G4Colour(1.0, 1.0, 1.0);

      // create line with colour
      G4VisAttributes attribs(colour);
      G4Polyline polyline;
      polyline.SetVisAttributes(attribs);

      // draw line
      G4Point3D start(fStep->GetPreStepPoint()->GetPosition());
      G4Point3D end(fStep->GetPostStepPoint()->GetPosition());
      polyline.push_back(start);
      polyline.push_back(end);
      pVVisManager->Draw(polyline);
    }
    
  }
  
}


