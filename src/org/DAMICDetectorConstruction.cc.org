#include "DAMICDetectorConstruction.hh"
#include "DAMICDetectorCopperBox.hh"
#include "DAMICDetectorInnerLead.hh"
#include "DAMICDetectorVessel.hh"
#include "DAMICDetectorOutterLead.hh"
#include "DAMICDetectorModules.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4MultiUnion.hh"
#include "G4Transform3D.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>


DAMICDetectorConstruction::DAMICDetectorConstruction():G4VUserDetectorConstruction(),fCheckOverlaps(true)
{
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* WorldMat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* Copper = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* Lead = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* Aluminium = nist->FindOrBuildMaterial("G4_Al");
  G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Galac = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* StainSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;


}

DAMICDetectorConstruction::~DAMICDetectorConstruction(){}

G4VPhysicalVolume* DAMICDetectorConstruction::Construct()
{

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


  G4Material* WorldMat = G4Material::GetMaterial("G4_Galactic");

  G4Box* WorldBox = new G4Box( "WorldBox", 1.1*m, 1.1*m, 1.1*m);

  G4LogicalVolume* WorldLV = new G4LogicalVolume(WorldBox, WorldMat, "WorldLV");


  G4VPhysicalVolume* WorldPV = new G4PVPlacement (0, G4ThreeVector(), WorldLV, "World",0, false, 0, fCheckOverlaps);


  G4LogicalVolume* VesselLV = GetConstructionVesselandOutter();


  G4PVPlacement* VesselPV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), VesselLV, "VesselPV", WorldLV, false, 0, fCheckOverlaps);

  return WorldPV;

}
