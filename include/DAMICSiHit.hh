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
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef DAMICSiHit_h
#define DAMICSiHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ProcessType.hh"
#include "G4ParticleDefinition.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DAMICSiHit : public G4VHit
{
  public:

      DAMICSiHit();
      ~DAMICSiHit();
      DAMICSiHit(const DAMICSiHit&);
      const DAMICSiHit& operator=(const DAMICSiHit&);
      int operator==(const DAMICSiHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
   
      void SetTrackID        (G4int de)          { trackID = de; };
      void SetStepID         (G4int de)          { stepID = de; };
      void SetEdep           (G4double de)       { edep = de; };
      void SetPos            (G4ThreeVector xyz) { pos = xyz; };
      void SetParticleName   (G4String name)     { particleName = name; };
      void SetParticlePDG    (G4int pdg)         { particlePDG = pdg; };
      void SetParticleEnergy (G4double e1)       { particleEnergy = e1; };
      void SetGlobalTime     (G4double t2)       { time = t2; };
      void SetLocalTime     (G4double t2)        { loctime = t2; };
  //      void SetParticle(G4ParticleDefinition* pa) {particle = pa;};
      void SetVolumeName(G4String name)          {volName = name;};
      void SetProcessName(G4String name)         {procName = name;};
      void SetCreatorProcess(G4String name)      {creatorProc = name;};
      void SetProcessType(G4ProcessType pa)      {procType = pa;};
      void SetStepLength(G4double pa)            {stepLength = pa;};
      void SetMotherID(G4int ma)                 {motherID = ma;};
  

      G4double GetTrackID()                      { return trackID; };
      G4double GetStepID()                       { return stepID; };
      G4double GetEdep()                         { return edep; };
      G4int GetParticlePDG()                     { return particlePDG;};
      G4ThreeVector GetPos()                     { return pos; };
      G4String GetParticleName()                 { return particleName;};
  //      G4ParticleDefinition* GetParticle()        { return particle;};
      G4double GetParticleEnergy()               { return particleEnergy;};
      G4double GetGlobalTime()                   { return time; };
      G4double GetLocalTime()                    { return loctime; };      
      G4ProcessType GetProcessType()             { return procType; } ;
      G4String GetProcessName()                  { return procName; } ;
      G4String GetCreatorProcess()               { return creatorProc; } ;
      G4String GetVolumeName()                   { return volName; } ;
      G4double GetStepLength()                   { return stepLength; };
      G4double GetMotherID()                     { return motherID; };      
  
  private:
  //      G4ParticleDefinition* particle; 
      G4int         trackID;
      G4int         stepID;
      G4double      edep;
      G4ThreeVector pos;
      G4double      loctime;
      G4double      time;
      G4int         motherID;
      G4String      particleName;
      G4String      creatorProc;
      G4double      particleEnergy;
      G4int         particlePDG;
      G4String volName;
      G4double stepLength;
      G4ProcessType procType;
      G4String procName;
  
      G4double hitEnergy;
      G4double firstHitTime;
  
      G4double firstParticleE;
      G4String firstParticleName;
      G4int firstTrackID;
      G4int firstStepID;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4THitsCollection<DAMICSiHit> DAMICSiHitsCollection;

extern G4ThreadLocal G4Allocator<DAMICSiHit> *DAMICSiHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* DAMICSiHit::operator new(size_t)
{
  if (!DAMICSiHitAllocator)
    DAMICSiHitAllocator = new G4Allocator<DAMICSiHit>;
  return (void*) DAMICSiHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void DAMICSiHit::operator delete(void* aHit)
{
  DAMICSiHitAllocator->FreeSingle((DAMICSiHit*) aHit);
}

#endif

