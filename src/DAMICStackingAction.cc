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
// StackingAction program
// --------------------------------------------------------------

#include "DAMICStackingAction.hh"

#include "DAMICStackingActionMessenger.hh"
#include "DAMICDetectorConstruction.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"



DAMICStackingAction::DAMICStackingAction() {

  // new messenger
  theMessenger = new DAMICStackingActionMessenger(this);

  // messenger defaults
  killGammasFlag  = 0;

  // global geometry navigator
  gNavigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
}


DAMICStackingAction::~DAMICStackingAction() {
  
  delete theMessenger; 
}


G4ClassificationOfNewTrack DAMICStackingAction::ClassifyNewTrack 
(const G4Track* aTrack) {


  static G4int secondariesKilled = 0;

  G4ClassificationOfNewTrack classification = fWaiting;

  // Kill secondaries in the concrete wall since they will be probably stopped by the shielding
  if(killGammasFlag) {
    // Make sure that we are not killing muons
    G4ParticleDefinition* particleType = aTrack->GetDefinition();
    if(particleType!=G4MuonPlus::Definition() && particleType!=G4MuonMinus::Definition()) {
      // check if particle is in world_phys
      G4ThreeVector pos = aTrack->GetPosition();
      G4ThreeVector *ptr = NULL;
      G4VPhysicalVolume *theVolume;
      theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,false);
      if(theVolume->GetName() == "world_phys") {
	classification = fKill;
	G4cout << " Total for session particles killed in concrete wall (world_phys): "
	       << ++secondariesKilled << G4endl;
      }
    }
  }

  return classification;

}


void DAMICStackingAction::NewStage() {;}

    
void DAMICStackingAction::PrepareNewEvent() {;}




