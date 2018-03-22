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
// PrimaryGeneratorAction header
// --------------------------------------------------------------

#ifndef DAMICPrimaryGeneratorMessenger_h
#define DAMICPrimaryGeneratorMessenger_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#ifdef DMXENV_GPS_USE
class G4GeneralParticleSource;
#else
class DAMICParticleSource;
#endif
class DAMICAnalysisManager;

class G4Event;

class DAMICPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

   public:
      DAMICPrimaryGeneratorAction();
      ~DAMICPrimaryGeneratorAction();
  
   public:
      void GeneratePrimaries(G4Event* anEvent);

   private:
#ifdef DMXENV_GPS_USE
  G4GeneralParticleSource* particleGun;
#else
  DAMICParticleSource* particleGun;
#endif

   private:
  //      const long* seeds;
  long seeds[2];
  G4double energy_pri;

   public:
  const long* GetEventSeeds() const       {return seeds;};
  G4double GetEnergyPrimary() const  {return energy_pri;};


};

#endif

