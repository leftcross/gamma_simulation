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
// EventAction header
// --------------------------------------------------------------

#ifndef DAMICEventAction_h
#define DAMICEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ios.hh"

#include "DAMICSiHit.hh"


class DAMICRunAction;
class DAMICPrimaryGeneratorAction;
class DAMICEventActionMessenger;
class DAMICAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class DAMICEventAction : public G4UserEventAction {

  public:
    DAMICEventAction();
    virtual ~DAMICEventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    void writeSiHitsToFile();
    void drawTracks(const G4Event*);

  public:
    void SetDrawTrksFlag (G4String val)     {drawTrksFlag    = val;};
    G4String GetDrawTrksFlag() const         {return drawTrksFlag;};

    void SetDrawColsFlag (G4String val)     {drawColsFlag    = val;};
    G4String GetDrawColsFlag() const  {return drawColsFlag;};

    void SetDrawHitsFlag (G4int val)        {drawHitsFlag    = val;};
    void SetSaveHitsFlag (G4int val)        {saveHitsFlag    = val;};
    void SetPrintModulo  (G4int val)        {printModulo     = val;};

  
  void SetPrimaryEvent(const G4Event *evt) {evt_pri = evt;};
    const G4Event* GetPrimaryEvent() const {return evt_pri;};
    
  void SetEventEnergy(G4double val) {energy_pri = val;};
    G4double GetEventEnergy() const {return energy_pri;};
  
  private:

    const G4Event *evt_pri;
    const long* seeds;
    G4int event_id;
    G4double energy_pri;
    G4double primX;
    G4double primY;
    G4double primZ;
    G4String primType;
    G4String primLV;
    
    // hits collections
    G4int ccdCollID;                
    G4int S_hits;
    G4int P_hits;

    // event summary
    G4double totEnergy;
    G4double totEnergyGammas;
    G4double totEnergyNeutrons;
    G4double hitEnergy;
    G4double firstHitTime;
    G4int trackID;
    G4int stepID;
  
    G4double firstParticleE;
    G4double particleEnergy;
    G4String particleName;
    G4String firstParticleName;
    G4int firstTrackID;
    G4int firstStepID;
  
    G4bool gamma_ev;
    G4bool neutron_ev;
    G4bool positron_ev;
    G4bool electron_ev;
    G4bool proton_ev;
    G4bool other_ev;
    G4bool start_gamma;
    G4bool start_neutron;

    // messenger
    G4String drawTrksFlag;
    G4String drawColsFlag;
    G4int drawHitsFlag;         
    G4int saveHitsFlag;         
    G4int printModulo;                         
    DAMICEventActionMessenger*  eventMessenger;

    const DAMICRunAction*    runAct;  //pointer to run action
    const DAMICPrimaryGeneratorAction* genAction; // pointer to particle generator
  std::ofstream *hitsfile;
};

#endif

