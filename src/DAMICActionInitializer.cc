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
// $Id: DAMICActionInitializer.cc 66241 2012-12-13 18:34:42Z gunter $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DAMICActionInitializer.hh"
#include "DAMICPrimaryGeneratorAction.hh"
#include "DAMICRunAction.hh"
#include "DAMICEventAction.hh"
#include "DAMICSteppingAction.hh"
#include "DAMICStackingAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DAMICActionInitializer::DAMICActionInitializer() : 
  G4VUserActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DAMICActionInitializer::Build() const 
{
  SetUserAction(new DAMICPrimaryGeneratorAction());
  SetUserAction(new DAMICRunAction());
  SetUserAction(new DAMICEventAction());
  SetUserAction(new DAMICSteppingAction());
  SetUserAction(new DAMICStackingAction());
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DAMICActionInitializer::BuildForMaster() const
{
  //Run action in the master
  SetUserAction(new DAMICRunAction());
}

