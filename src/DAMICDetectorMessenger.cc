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
// SteppingActionMessenger program
// --------------------------------------------------------------

#include "DAMICDetectorMessenger.hh"

#include "DAMICDetectorConstruction.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4RunManager.hh"


DAMICDetectorMessenger::DAMICDetectorMessenger
   (DAMICDetectorConstruction* DC):detectorConstruction(DC) {
  
  EKineCutCmd = new G4UIcmdWithADoubleAndUnit("/damic/MinEnergyCut",this);
  EKineCutCmd->SetGuidance("Minimum Charged particle cut inside detector");
  EKineCutCmd->SetParameterName("ECut",false,false);
  EKineCutCmd->SetRange("ECut>=250.0*eV");
  EKineCutCmd->SetDefaultUnit("eV");
  EKineCutCmd->SetUnitCategory("Energy");
  EKineCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TimeCutCmd = new G4UIcmdWithADoubleAndUnit("/damic/TimeCut",this);
  TimeCutCmd->SetGuidance("Set Time Cut (for neutrons) inside detector");
  TimeCutCmd->SetParameterName("TCut",false,false);
  TimeCutCmd->SetRange("TCut>0.");
  TimeCutCmd->SetDefaultUnit("ns");
  TimeCutCmd->SetUnitCategory("Time");
  TimeCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}


//ooooooooooooooooooooooooooooooooooooooooo
DAMICDetectorMessenger::~DAMICDetectorMessenger() 
{
  delete EKineCutCmd;
  delete TimeCutCmd;
 }


//ooooooooooooooooooooooooooooooooooooooooo
void DAMICDetectorMessenger::SetNewValue(G4UIcommand* command, 
				       G4String newValue) 
{

  if(command == EKineCutCmd)
   detectorConstruction->
     SetEnergyCut(EKineCutCmd->GetNewDoubleValue(newValue));

  if(command == TimeCutCmd)
   detectorConstruction->
     SetTimeCut(TimeCutCmd->GetNewDoubleValue(newValue));

  //trigger a re-optimization of the geometry
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->GeometryHasBeenModified();


}



