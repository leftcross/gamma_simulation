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
//    For information related to this code contact: M. Settimo
//    e-mail: settimo@subatech.in2p3.fr
// --------------------------------------------------------------
// Comments
//
//                   Underground Advanced -  Modified by M.Settimo 2017
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//               Modified by M.Settimo 2017
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DAMICDetectorConstruction.hh"
#include "DAMICDetectorMessenger.hh"

#include "DAMICSiSD.hh"
#include "DAMICDetectorConstruction.hh"
#include "DAMICDetectorModules.hh"
#include "G4NistManager.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSolid.hh" 

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DAMICDetectorConstruction::DAMICDetectorConstruction()  
{
  // create commands for interactive definition of time cuts:
  detectorMessenger = new DAMICDetectorMessenger(this);

  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theMinEkine         = 0.0*eV; // minimum kinetic energy required in volume

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DAMICDetectorConstruction::~DAMICDetectorConstruction() 
{

  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DAMICDetectorConstruction::DefineMaterials() 
{

#include "DAMICDetectorMaterial.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DAMICDetectorConstruction::Construct() {

  DefineMaterials();

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;


  // Universe - room wall - CONCRETE ****************************************
  //NB: measured INSIDE of lab, therefore have to add twice wall thickness

  G4double Offset_DAMICGeoZ = 0 ;  
  G4double wallThick   = 5.0*cm;
  G4double worldWidth  = 800.0*cm + 2.*wallThick; // "x"
  G4double worldLength = 800.0*cm + 2.*wallThick; // "y"
  G4double worldHeight = 800.0*cm + 2.*wallThick; // "z"


  G4Box* WorldBox = new G4Box("world", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight);
  WorldLV  = new G4LogicalVolume(WorldBox, WorldMat, "WorldLV");
  WorldPV = new G4PVPlacement(0, G4ThreeVector(0.,0.,0), "WorldPV", WorldLV, NULL, false,0);

  G4VisAttributes* world_vat= new G4VisAttributes(orange);
  // WorldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  world_vat->SetVisibility(true);
  //world_vat->SetVisibility(false);
  //world_log->SetVisAttributes(world_vat);
  
  
  // Lab Space - AIR ********************************************************

  G4double labWidth  = worldWidth  - 2.*wallThick; //X
  G4double labLength = worldLength - 2.*wallThick; //Y
  G4double labHeight = worldHeight - 2.*wallThick; //Z

  G4Box* LabBox = new G4Box("Lab_Box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight);
  LabLV = new G4LogicalVolume(LabBox, LabMat, "LabLV");

  G4VisAttributes* lab_vat= new G4VisAttributes(gray);
  //  lab_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  lab_vat->SetVisibility(true);
  lab_vat->SetVisibility(true);
  LabLV->SetVisAttributes(lab_vat);
  
//  ----- Lead shielding  ---- 
  G4double OutLeadBoxX = (2+6+2)*2.54*cm;
  G4double OutLeadBoxY = (4+8+4)*2.54*cm;
  G4double OutLeadBoxZ = (8+6+6)*2.54*cm;

  G4double InLeadBoxX = (8+0.1)*2.54*cm; //defining only 1 type of lead 
  G4double InLeadBoxY = (8+0.1)*2.54*cm;
  G4double InLeadBoxZ = (8+0.1)*2.54*cm;
  G4double LeadBoxThickness = (6+2)*2.54*cm; //6" of Fermilab lead + 2 ancient lead =  thickness from each side
  
  
   G4Box *outerLeadBox = new G4Box("Outer Box", OutLeadBoxX/2.,OutLeadBoxY/2.,OutLeadBoxZ/2.);   
   extLeadBoxLV = new G4LogicalVolume(outerLeadBox, LeadMat, "extLeadBoxLV");
   G4VisAttributes* lead_vat= new G4VisAttributes(blue);
   lead_vat->SetVisibility(true);
   extLeadBoxLV->SetVisAttributes(lead_vat); 


   //--- empty space in the hollow lead box --
    G4Box *emptyLeadBox = new G4Box("emptyLeadBox", InLeadBoxX/2., InLeadBoxY/2., InLeadBoxZ/2.);
    G4LogicalVolume * emptyLeadBoxLV = new G4LogicalVolume(emptyLeadBox, VacuumMat, "emptyLeadBoxLV");
   G4VisAttributes* vacuumlead_vat= new G4VisAttributes(blue);
   vacuumlead_vat->SetVisibility(true);
   emptyLeadBoxLV->SetVisAttributes(vacuumlead_vat); 

  G4double OutSteelBoxX= 6*2.54*cm;
  G4double OutSteelBoxY= 6*2.54*cm;
  G4double OutSteelBoxZ= 6*2.54*cm;
  
 // G4double InSteelBoxX= 40*cm;
 // G4double InSteelBoxY= 40*cm;
 // G4double InSteelBoxZ= 40*cm;

  G4Box* outerSteelBox= new G4Box("Outer chamber",OutSteelBoxX/2., OutSteelBoxY/2., OutSteelBoxZ/2. );
  
  G4Sphere* emptySteelSphere=new G4Sphere("empty steel sphere",0.,(3*2.54*cm), 0*degree, 360*degree, 0*degree, 180*degree);

//  G4VSolid* boxminussphere=new G4SubtractionSolid("cube-sphere",outerSteelBox,emptySteelSphere); 

  G4LogicalVolume * extSteelBoxLV = new G4LogicalVolume(outerSteelBox, StainSteelMat, "extSteelBoxLV");
  G4VisAttributes* steel_vat=new G4VisAttributes(green);
  steel_vat->SetVisibility(true);
  extSteelBoxLV->SetVisAttributes(steel_vat);

 

 // G4Box *emptySteelBox = new G4Box("emptySteelBox",InSteelBoxX/2.,InSteelBoxY/2.,InSteelBoxZ/2.);
 // G4LogicalVolume * emptySteelBoxLV = new G4LogicalVolume(emptySteelBox,VacuumMat, "emptySteelBoxLV");
  
  G4double rad_cyl_cavity=7*2.54*cm;
  G4double thick_cyl_cavity=1*2.54*cm;
  G4VSolid* cylinder1=new G4Tubs("Cylinder1",0.,rad_cyl_cavity,thick_cyl_cavity,0.,360*degree);
  
 // G4VSolid* union1=new G4UnionSolid("Cylinder1+Sphere",emptySteelSphere,cylinder1); 

  G4LogicalVolume* emptySteelSphereLV= new G4LogicalVolume(emptySteelSphere, VacuumMat, "emptySteelSphereLV");
  G4VisAttributes* vacuumsteel_vat=new G4VisAttributes(green);
  vacuumsteel_vat->SetVisibility(true);
  emptySteelSphereLV->SetVisAttributes(vacuumsteel_vat);

  LabPV = new G4PVPlacement(0, G4ThreeVector(0.,0.,0), "LabPV", LabLV, WorldPV, false,0);
  G4PVPlacement* extLeadBoxPV = new G4PVPlacement(0, G4ThreeVector(0.,0., 0.),
						     "extLeadBoxPV", extLeadBoxLV, LabPV, false,true);
  G4PVPlacement* emptyLeadBoxPV = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
						  "emptyLeadBoxPV", emptyLeadBoxLV, extLeadBoxPV, false,true);
 
  G4double distance_chamber=150*cm;
 
 G4PVPlacement* extSteelBoxPV = new G4PVPlacement(0, G4ThreeVector(0.,0.,distance_chamber), "extSteelBoxPV", extSteelBoxLV, LabPV, false, true);	
 G4PVPlacement* emptySteelSpherePV= new G4PVPlacement(0, G4ThreeVector(0*cm, 0., 0*cm), "emptySteelSpherePV", emptySteelSphereLV, extSteelBoxPV, false, true);   
 

//G4double r_min=2*2.54*cm;
G4double r_max=4.*2.54*cm;

G4Sphere* extbonner= new G4Sphere("bonner sphere",0., r_max, 0*degree, 360*degree, 0*degree, 180*degree); 
G4LogicalVolume* extBonnerSphereLV= new G4LogicalVolume(extbonner,PolyMat,"extBonnerSphereLV");
G4VisAttributes* extbonnerPoly_vat=new G4VisAttributes(red);
extbonnerPoly_vat->SetVisibility(true);
extBonnerSphereLV->SetVisAttributes(extbonnerPoly_vat);


//G4Sphere* emptybonner= new G4Sphere("inner bonner sphere", 0. , r_min, 0*degree, 360*degree, 0*degree, 180*degree);
//G4LogicalVolume* emptyBonnerSphereLV= new G4LogicalVolume(emptybonner, VacuumMat, "emptyBonnerSphereLV");
//G4VisAttributes* emptybonnerVacuum_vat=new G4VisAttributes(lblue);
//emptybonnerVacuum_vat->SetVisibility(true);
//emptyBonnerSphereLV->SetVisAttributes(emptybonnerVacuum_vat);


G4PVPlacement* extBonnerSpherePV= new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),"extBonnerSpherePV",extBonnerSphereLV,emptyLeadBoxPV, false, true);

G4double tableX=2.54*cm;
G4double tableY=95*2.54*cm;
G4double tableZ=36*2.54*cm;

G4Box* tablePolyBox= new G4Box("Upper table poly",tableX/2., tableY/2., tableZ/2.);
G4LogicalVolume* polyTableLV=new G4LogicalVolume(tablePolyBox,PolyMat,"polyTableLV");
G4VisAttributes* polyTable_vat=new G4VisAttributes(red);
polyTable_vat->SetVisibility(true);
polyTableLV->SetVisAttributes(polyTable_vat);

G4double distance_table=distance_chamber+25*cm;
G4PVPlacement* polyTablePV=new G4PVPlacement(0, G4ThreeVector(-1*OutSteelBoxX/2.-tableX/2.,0.,distance_table),"polyTablePV",polyTableLV,LabPV,false,true);

G4double thickness_plastic=0.4*cm;
G4Box* tableWoodBox= new G4Box("upper table wood", (tableX/2.-thickness_plastic),(tableY/2.0 - thickness_plastic), (tableZ/2.0 - thickness_plastic));  
G4LogicalVolume* woodTableLV= new G4LogicalVolume(tableWoodBox,TableMat,"woodTableLV");
G4VisAttributes* woodTable_vat=new G4VisAttributes(orange);
woodTable_vat->SetVisibility(true);
woodTableLV->SetVisAttributes(woodTable_vat);

G4PVPlacement* woodTablePV= new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"woodTablePV",woodTableLV,polyTablePV,false,true);

//G4PVPlacement* emptyBonnerSpherePV= new G4PVPlacement(0,G4ThreeVector(0.,0.,0.), "emptyBonnerSpherePV", emptyBonnerSphereLV, extBonnerSpherePV, false, true);


// -------------------------------------------MODULES---------------------------------------

  // CCD01_44LV = GetConstructionCCDSensor44();
  // CCD02_44LV = GetConstructionCCDSensor44();
  // CCD03_44LV = GetConstructionCCDSensor44();
  // CCD04_44LV = GetConstructionCCDSensor44();
  // CCD05_44LV = GetConstructionCCDSensor44();
  // CCD06_44LV = GetConstructionCCDSensor44();
  // //CCD07_44LV = GetConstructionCCDSensor44();
  
  // G4double PosZModule44 = 6.35/2*mm + 5.994/2*mm+3*mm; //PosZBottomPlate + 6.35/2*mm + 5.994/2*mm+3*mm;
  // G4double PosYModule44 =  6.35/2*mm - 116.332/2*mm -1.800*mm; //PosYRearPlate - 6.35/2*mm - 116.332/2*mm -1.800*mm;

  // G4ThreeVector Module44Vect = G4ThreeVector(0,PosYModule44,PosZModule44);
  // G4ThreeVector uModule44 = G4ThreeVector(0,-1,0);
  // G4ThreeVector vModule44 = G4ThreeVector(1,0,0);
  // G4ThreeVector wModule44 = G4ThreeVector(0,0,1);
  // G4RotationMatrix* Module44Rot = new G4RotationMatrix(uModule44, vModule44, wModule44);

  // G4PVPlacement* CCD01PV = new G4PVPlacement(Module44Rot, G4ThreeVector(0,PosYModule44,PosZModule44+ 00*mm), CCD01_44LV, "CCD01PV", emptyLeadBoxLV, false, 0, false);
  // G4PVPlacement* CCD02PV = new G4PVPlacement(Module44Rot, G4ThreeVector(0,PosYModule44,PosZModule44+ 20*mm), CCD02_44LV, "CCD02PV", emptyLeadBoxLV, false, 0, false);
  // G4PVPlacement* CCD03PV = new G4PVPlacement(Module44Rot, G4ThreeVector(0,PosYModule44,PosZModule44+ 40*mm), CCD03_44LV, "CCD03PV", emptyLeadBoxLV, false, 0, false);
  // G4PVPlacement* CCD04PV = new G4PVPlacement(Module44Rot, G4ThreeVector(0,PosYModule44,PosZModule44+ 60*mm), CCD04_44LV, "CCD04PV", emptyLeadBoxLV, false, 0, false);
  // G4PVPlacement* CCD05PV = new G4PVPlacement(Module44Rot, G4ThreeVector(0,PosYModule44,PosZModule44+ 80*mm), CCD05_44LV, "CCD05PV", emptyLeadBoxLV, false, 0, false);
  // G4PVPlacement* CCD06PV = new G4PVPlacement(Module44Rot, G4ThreeVector(0,PosYModule44,PosZModule44+100*mm), CCD06_44LV, "CCD06PV", emptyLeadBoxLV, false, 0, false);
  
 
  return WorldPV;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICDetectorConstruction::ConstructSDandField()
{
  // ......................................................................
  // sensitive detectors ..................................................
  // ......................................................................
  
  // if (SiSD.Get() == 0) 
  //  {    
  //     G4String name="/DMXDet/SiSD";
  //     DAMICSiSD* aSD = new DAMICSiSD(name);
  //     SiSD.Put(aSD);
  //   }
  
  // G4SDManager::GetSDMpointer()->AddNewDetector(SiSD.Get());  
  // if ( CCD01_44LV)
  //   SetSensitiveDetector(CCD01_44LV,SiSD.Get());
  // if ( CCD02_44LV)
  //   SetSensitiveDetector(CCD02_44LV,SiSD.Get());
  // if ( CCD03_44LV)
  //   SetSensitiveDetector(CCD03_44LV,SiSD.Get());
  // if ( CCD04_44LV)
  //   SetSensitiveDetector(CCD04_44LV,SiSD.Get());
  // if ( CCD05_44LV)
  //   SetSensitiveDetector(CCD05_44LV,SiSD.Get());
  // if ( CCD06_44LV)
  //   SetSensitiveDetector(CCD06_44LV,SiSD.Get());  
 
  // return;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void DAMICDetectorConstruction::SetEnergyCut(G4double val)
{

}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DAMICDetectorConstruction::SetTimeCut(G4double val)
{
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



