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
// 17 Jan 2002 Alex Howard Added Analysis
// 23 Oct 2009 Luciano Pandola Removed un-necessary calls from EndOfRun()
//
// RunAction program
// --------------------------------------------------------------

#include "DAMICRunActionMessenger.hh"
#include "DAMICRunAction.hh"


#include "G4Run.hh"
#include "G4ios.hh"

#include <fstream>

#include "DAMICAnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DAMICRunAction::DAMICRunAction()
{
  runMessenger = new DAMICRunActionMessenger(this);
  savehitsFile = "hits.out";
  savehistFile = "damic";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DAMICRunAction::~DAMICRunAction()
{
  delete runMessenger;
  runMessenger = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICRunAction::BeginOfRunAction(const G4Run* aRun)
{
  //Master mode or sequential
  if (IsMaster())    
    G4cout << "### Run " << aRun->GetRunID() << " starts (master)." << G4endl;
  else
    G4cout << "### Run " << aRun->GetRunID() << " starts (worker)." << G4endl;
  
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  
  // Open an output file
  man->OpenFile(savehistFile);

  // Book histograms and ntuples
  Book();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICRunAction::EndOfRunAction(const G4Run*)
{
  // Merge parameters
  //G4ParameterManager* parameterManager = G4ParameterManager::Instance();
  //parameterManager->Merge();

  //G4AnalysisManager* man = G4AnalysisManager::Instance();
  //man->Write();
  //man->CloseFile();  
  G4cout << " done writing!  " ;
  G4cout << G4endl;
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



void DAMICRunAction::Book()
{  
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  
  // Open an output file
  //man->OpenFile(savehistFile);
  man->SetFirstHistoId(1);
  man->SetFirstNtupleId(0);


  // ---- primary ntuple ------
  // id==1
    man->CreateNtuple("tree1", "Particle Source Energy");
    //primary infos
    man->CreateNtupleIColumn("EventID"); //0
    //  man->CreateNtupleSColumn("primType"); //1
    man->CreateNtupleDColumn("primEnergy"); //2
    man->CreateNtupleDColumn("primX"); //3 
    man->CreateNtupleDColumn("primY"); //4 
    man->CreateNtupleDColumn("primZ"); //5 
    //   man->CreateNtupleSColumn("primLV"); //6 Logical Volume
    //current particle infos
    man->CreateNtupleIColumn("trackID"); //7
    man->CreateNtupleIColumn("stepID"); //8
    man->CreateNtupleIColumn("partPDG"); //9 
    //    man->CreateNtupleSColumn("partType"); //10
    man->CreateNtupleDColumn("partEtot"); //11 
    man->CreateNtupleDColumn("partEkin"); //12
    man->CreateNtupleDColumn("partX");  //13
    man->CreateNtupleDColumn("partY");  //14 
    man->CreateNtupleDColumn("partZ");  //15
    // man->CreateNtupleDColumn("partTau");//16
    // man->CreateNtupleDColumn("partGlobalTime");//17
    // man->CreateNtupleDColumn("partProperTime");//18
    // man->CreateNtupleDColumn("partTrackLength"); //19
    man->CreateNtupleSColumn("partCreatorProcess");  //20 //at the current step

    //step-level
    man->CreateNtupleSColumn("stepProcess");  //21 //at the current step
    man->CreateNtupleDColumn("stepEdep"); //22
    man->CreateNtupleDColumn("stepEdep_nonIon"); //23
    man->CreateNtupleDColumn("stepLength"); //24

    man->CreateNtupleSColumn("stepVolume");  //24 //from G4start
    man->CreateNtupleSColumn("stepNextVolume");  //25 //from G4start
    man->CreateNtupleSColumn("stepMaterial");  //26 //from G4start
    man->CreateNtupleSColumn("stepNextMaterial");  //27 //from G4start
    man->CreateNtupleIColumn("stepNsecondaries");  //28 //from G4start
        //mother info - particle Change
    man->CreateNtupleIColumn("partMotherID"); //29  : 0,1 if there is a part change
//    man->CreateNtupleDColumn("partTheta");   // 30, particle theta 
//    man->CreateNtupleDColumn("partPhi");     // 31, particle phi
    man->CreateNtupleSColumn("prodVolume");
    man->CreateNtupleDColumn("prodX");
    man->CreateNtupleDColumn("prodY");
    man->CreateNtupleDColumn("prodZ");
    man->CreateNtupleIColumn("isEntering");

    man->FinishNtuple();

  // ---- tertiary ntuple ------   
  //id==3
  man->CreateNtuple("tree2", " Hits Info");
  man->CreateNtupleIColumn("trackID");
  man->CreateNtupleIColumn("stepID");
  man->CreateNtupleIColumn("Nhits");
  man->CreateNtupleDColumn("Edep");
  man->CreateNtupleDColumn("TotEdep");
  man->CreateNtupleDColumn("xpos");
  man->CreateNtupleDColumn("ypos");
  man->CreateNtupleDColumn("zpos");
  man->CreateNtupleDColumn("stepLength");
  man->CreateNtupleDColumn("globalTime");
  man->CreateNtupleDColumn("localTime");
  man->CreateNtupleDColumn("partEnergy");
  man->CreateNtupleSColumn("partName");
  man->CreateNtupleIColumn("partPDG");
  man->CreateNtupleSColumn("volName");
  man->CreateNtupleSColumn("procName");
  man->CreateNtupleIColumn("procType");
  man->CreateNtupleSColumn("creatorProcess");
  man->CreateNtupleIColumn("motherID");  
  man->FinishNtuple();


    //id==3
  man->CreateNtuple("tree3", " Hits Info - evt level");
  man->CreateNtupleDColumn("event");
  man->CreateNtupleDColumn("hits");
  man->CreateNtupleDColumn("xpos");
  man->CreateNtupleDColumn("ypos");
  man->CreateNtupleDColumn("zpos");
  man->FinishNtuple();
  
  return;

}


