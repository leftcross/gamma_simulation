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
// EventActionMessenger program
// --------------------------------------------------------------

#include "DAMICEventActionMessenger.hh"

#include <sstream>

#include "DAMICEventAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "globals.hh"


DAMICEventActionMessenger::DAMICEventActionMessenger(DAMICEventAction* EvAct)
:eventAction(EvAct){

  // saving event information
  damicDirectory = new G4UIdirectory("/damic/");
  damicDirectory->SetGuidance("DM Example commands.");

  SaveHitsCmd = new G4UIcmdWithABool("/damic/saveHits",this);
  SaveHitsCmd->SetGuidance("Set flag to save hits in each run");
  SaveHitsCmd->SetGuidance("into file 'hits.out'");
  SaveHitsCmd->SetGuidance("Default = true");
  SaveHitsCmd->SetParameterName("saveHitsFlag", false);


  // drawing event
  drawDirectory = new G4UIdirectory("/damic/draw/");
  drawDirectory->SetGuidance("DM Example draw commands.");

  DrawColsCmd = new G4UIcmdWithAString("/damic/draw/drawColours",this);
  DrawColsCmd->SetGuidance("Tracks drawn by Event (standard colours) or by Step (custom colours)");
  DrawColsCmd->SetGuidance("  Choice : custom, standard(default)");
  DrawColsCmd->SetParameterName("drawColsFlag", false);
  DrawColsCmd->SetCandidates("custom standard");
  DrawColsCmd->AvailableForStates(G4State_Idle);

  DrawTrksCmd = new G4UIcmdWithAString("/damic/draw/drawTracks",this);
  DrawTrksCmd->SetGuidance("Draw the tracks in the event");
  DrawTrksCmd->SetGuidance("  Choice : none, charged, noscint, all(default)");
  DrawTrksCmd->SetParameterName("drawTrksFlag", false);
  DrawTrksCmd->SetCandidates("none charged noscint all");
  DrawTrksCmd->AvailableForStates(G4State_Idle);

  DrawHitsCmd = new G4UIcmdWithABool("/damic/draw/drawHits",this);
  DrawHitsCmd->SetGuidance("Set flag to draw hits in PMT.");
  DrawHitsCmd->SetGuidance("Default = true");
  DrawHitsCmd->SetParameterName("drawHitsFlag", false);
  DrawHitsCmd->SetDefaultValue(true);

  PrintCmd = new G4UIcmdWithAnInteger("/damic/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);     

  
}


DAMICEventActionMessenger::~DAMICEventActionMessenger() {


  delete SaveHitsCmd;  
  delete damicDirectory;
  delete DrawColsCmd;
  delete DrawTrksCmd;
  delete DrawHitsCmd;  
  delete drawDirectory;

  delete PrintCmd;

}

void DAMICEventActionMessenger::SetNewValue
   (G4UIcommand* command, G4String newValue) { 

  
  if(command == DrawColsCmd)
    eventAction->SetDrawColsFlag(newValue);

  if(command == DrawTrksCmd)
    eventAction->SetDrawTrksFlag(newValue);

  if(command == DrawHitsCmd) {
    G4int vl;
    const char* t = newValue;
    std::istringstream is(t);
    is >> vl;
    eventAction->SetDrawHitsFlag(vl!=0);
  }


  if(command == SaveHitsCmd) {
    G4int vl;
    const char* t = newValue;
    std::istringstream is(t);
    is >> vl;
    eventAction->SetSaveHitsFlag(vl!=0);
  }

  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}


}

