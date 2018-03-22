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


G4LogicalVolume* GetConstructionModule44(){
  G4double ModuleX = 116.332*mm;
  G4double ModuleY = 116.84*mm;
  G4double ModuleZ = 5.994*mm;

  G4double erasePartX = 111.252*mm;
  G4double erasePartY = (ModuleY-102.616)/2;
  G4double erasePartZ = ModuleZ;

  G4Box* ModuleBox = new G4Box("ModuleBox", ModuleX/2, ModuleY/2, ModuleZ/2);
  G4Box* Box1 = new G4Box("Box1", erasePartX/2, erasePartY/2, erasePartZ);

  G4RotationMatrix* rot0Erase = new G4RotationMatrix;
  G4double posYErase = 102.616/2+erasePartY/2;
  G4double posXErase = ModuleX/2 - erasePartX/2+0.005;
  G4ThreeVector vectErase1 = G4ThreeVector(posXErase,posYErase,0);
  G4ThreeVector vectErase2 = G4ThreeVector(posXErase,-posYErase,0);
  G4Transform3D TrErase1 = G4Transform3D(*rot0Erase, vectErase1);
  G4Transform3D TrErase2 = G4Transform3D(*rot0Erase, vectErase2);

  G4SubtractionSolid* sub1 = new G4SubtractionSolid("Sub1",ModuleBox, Box1, TrErase1);
  G4SubtractionSolid* sub2 = new G4SubtractionSolid("Sub2",sub1, Box1, TrErase2);

  G4Material* Galac = G4Material::GetMaterial("G4_Galactic");
  G4LogicalVolume* Module44LV = new G4LogicalVolume(sub2, Galac,"Module44LV");

  G4LogicalVolume* CopperBasePlateLV = GetConstructionCopperBasePlate();
  G4LogicalVolume* CopperMountBarLV = GetConstructionCopperMountBar();
  //G4LogicalVolume* ProtectiveCoverLV = GetConstructionProtectiveCover();
  G4LogicalVolume* SubAssyLV = GetConstructionCCDSubAssy44();
  G4LogicalVolume* CopperTopPlate44LV = GetConstructionCopperTopPlate44();


  /*--------PosCopperMountBar-----*/
  G4double PosBarZ = -ModuleZ/2+5.334/2;
  G4double PosBarY = 0;
  G4double PosBarX = -ModuleX/2+11.43/2*mm;
  G4ThreeVector VectPosBar = G4ThreeVector(PosBarX, PosBarY, PosBarZ);
  G4ThreeVector uPosBar = G4ThreeVector(0,0,1);
  G4ThreeVector vPosBar = G4ThreeVector(1,0,0);
  G4ThreeVector wPosBar = G4ThreeVector(0,1,0);
  G4RotationMatrix* RotPosBar = new G4RotationMatrix(uPosBar, vPosBar, wPosBar);
   /*-----------PosCopperBasePlate-------*/
  G4double PosZCopperBasePlate = PosBarZ-2.006/2;
  G4double PosYCopperBasePlate = 0;
  G4double PosXCopperBasePlate = PosBarX - 0.635*mm +111.252/2*mm;
  G4ThreeVector VectCopperBasePlate = G4ThreeVector(PosXCopperBasePlate,PosYCopperBasePlate,PosZCopperBasePlate);

  /*---------PosSubAssy-----*/
  G4double PosSubAssyZ = PosZCopperBasePlate + 1.206/2*mm+1.118*mm ;
  G4double PosSubAssyY =0;
  G4double PosSubAssyX = PosXCopperBasePlate + 12.891*mm - 9.505*mm;
  G4ThreeVector VectSubAssy = G4ThreeVector(PosSubAssyX,PosSubAssyY,PosSubAssyZ);
  G4ThreeVector uAssy = G4ThreeVector(-1,0,0);
  G4ThreeVector vAssy = G4ThreeVector(0,-1,0);
  G4ThreeVector wAssy = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotSubAssy = new G4RotationMatrix(uAssy, vAssy, wAssy);

  /*-------PosCopperTopPlate-----*/
  G4double PosZTop = PosSubAssyZ + 1.206/2*mm+2.006/2*mm;
  G4double PosYTop = 0;
  G4double PosXTop = PosXCopperBasePlate+ 111.252/2*mm-85.852/2*mm;
  G4ThreeVector VectTopPlate = G4ThreeVector(PosXTop,PosYTop,PosZTop);
  /*------PosProtectiveCover-----*/
  G4double PosProtectiveZ = PosSubAssyZ + 0.6*mm + 0.25/2*mm+1*mm;
  G4double PosProtectiveY = 0;
  G4double PosProtectiveX = 11.1125*mm;
  G4ThreeVector VectProtective = G4ThreeVector(PosProtectiveX,PosProtectiveY,PosProtectiveZ);

  G4PVPlacement* CopperBasePlatePV = new G4PVPlacement(0, VectCopperBasePlate, CopperBasePlateLV, "CopperBasePlatePV", Module44LV, false, 0, false);
  G4PVPlacement* CopperMountBarPV = new G4PVPlacement(RotPosBar, VectPosBar, CopperMountBarLV, "CopperMountBarPV", Module44LV, false, 0, false);
  //G4PVPlacement* ProtectiveCoverPV = new G4PVPlacement(0, VectProtective, ProtectiveCoverLV, "ProtectiveCoverPV", Module44LV, false, 0, false);
  G4PVPlacement* SubAssyPV = new G4PVPlacement(RotSubAssy, VectSubAssy, SubAssyLV, "SubAssyPV", Module44LV, false, 0, false);
  G4PVPlacement* CopperTopPlate44PV = new G4PVPlacement(0, VectTopPlate, CopperTopPlate44LV, "CopperTopPlate44PV", Module44LV, false, 0, false);



  return Module44LV;
}

G4LogicalVolume* GetConstructionCCDSubAssy44()
{
  G4double SubX = 95.21*mm;
  G4double SubY = 101.6*mm;
  G4double SubZ = 1.206*mm;
  G4Box* SubBox = new G4Box("SubBox", SubX/2, SubY/2, SubZ/2);

  G4Material* Galac = G4Material::GetMaterial("G4_Galactic");
  G4LogicalVolume* SubLV = new G4LogicalVolume(SubBox, Galac,"Sub44LV");


  G4LogicalVolume* CCDSensor44LV = GetConstructionCCDSensor44();
  G4LogicalVolume* SiliconSubstrate44LV = GetConstructionSiliconSubstrate44();
  G4LogicalVolume* FlexCable44LV = GetConstructionFlexCable44();


  G4double CopyNeeded = (SubX-90.358)/2-0.02;
  /*---------SensorPos----------*/
  G4double PosZSensor = -SubZ/2+0.25/2;
  G4double PosYSensor = 0;
  G4double PosXSensor = -SubX/2+76.2/2*mm;
  G4ThreeVector VectSensor = G4ThreeVector(PosXSensor, PosYSensor, PosZSensor);
  /*----------FlexCablePos-----------*/
  G4double PosZFlexCable = PosZSensor+ 0.25/2*mm+0.305/2*mm;
  G4double PosYFlexCable = 0;
  G4double PosXFlexCable = CopyNeeded;
  G4ThreeVector VectFlexCable = G4ThreeVector(PosXFlexCable,PosYFlexCable , PosZFlexCable);

  /*-----------SiliconSubstratePos----------*/
  G4double PosZSiliconSubstrate = PosZFlexCable+0.305/2*mm+0.65/2*mm;
  G4double PosYSiliconSubstrate = 0;
  G4double PosXSiliconSubstrate = -SubX/2+76.2/2*mm;
  G4ThreeVector VectSiliconSubstrate = G4ThreeVector(PosXSiliconSubstrate, PosYSiliconSubstrate, PosZSiliconSubstrate);


  G4PVPlacement* CCDSensor44PV = new G4PVPlacement(0, VectSensor, CCDSensor44LV, "CCDSensor44PV", SubLV, false, 0, false);
  G4PVPlacement* SiliconSubstrate44PV = new G4PVPlacement(0, VectSiliconSubstrate, SiliconSubstrate44LV, "SiliconSubstrate44PV", SubLV, false, 0, false);
  G4PVPlacement* FlexCable44PV = new G4PVPlacement(0, VectFlexCable, FlexCable44LV, "FlexCable44PV", SubLV, false, 0, false);

  return SubLV;
}

G4LogicalVolume* GetConstructionCopperBasePlate()
{
  G4double boxXBase = 111.252*mm;
  G4double boxYBase = 102.616*mm;
  G4double boxZBase = 2.006*mm;
  G4double rmBoxX = 77.216*mm;
  G4double rmBoxY = 86.766*mm;
  G4double rmBoxZ = 0.889*mm;
  G4double rm2BoxX = 56.769*mm;
  G4double rm2BoxY = 54.864*mm;
  G4double rm2BoxZ = 2.0066*mm;

  G4Box* boxBase = new G4Box("BoxBase", boxXBase/2, boxYBase/2, boxZBase/2);
  G4Box* rmBox = new G4Box("rmBox", rmBoxX/2, rmBoxY/2, rmBoxZ);
  G4Box* rm2Box = new G4Box("rm2Box", rm2BoxX/2, rm2BoxY/2, rm2BoxZ);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4double posRM2BoxX = 11.938*mm;
  G4double posRM2BoxY = 0*mm;
  G4ThreeVector posRM2BoxVect = G4ThreeVector(posRM2BoxX, posRM2BoxY, 0);
  G4Transform3D posRM2BoxTr = G4Transform3D(*rot0, posRM2BoxVect);

  G4double posRMBoxX = 11.938*mm;
  G4double posRMBoxZ = 2.006/2*mm;
  G4ThreeVector posRMBoxVect = G4ThreeVector(posRMBoxX, 0, posRMBoxZ);
  G4Transform3D posRMBoxTr = G4Transform3D(*rot0, posRMBoxVect);

  G4SubtractionSolid* copperBase1 = new G4SubtractionSolid("CopperBase1", boxBase, rmBox, posRMBoxTr);
  G4SubtractionSolid* copperBaseFinal = new G4SubtractionSolid("CopperBaseFinal", copperBase1, rm2Box, posRM2BoxTr);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *copperBasePlate= new G4LogicalVolume(copperBaseFinal, Copper,"CopperBasePlate");

  return copperBasePlate;

}


G4LogicalVolume* GetConstructionCopperMountBar(){

G4double boxBaseX = 116.84*mm;
G4double boxBaseY = 5.334*mm;
G4double boxBaseZ = 11.43*mm;

G4double rm1BoxX = 28.448*mm;
G4double rm1BoxY = 1.524*mm;
G4double rm1BoxZ = 11.43*mm;

G4double rm2BoxX = 116.84*mm;
G4double rm2BoxY = 2.794*mm;
G4double rm2BoxZ = 6.35*mm;

G4double rm3BoxX = 7.62*mm;
G4double rm3BoxY = 5.334*mm;
G4double rm3BoxZ = 6.35*mm;


G4Box* boxBase = new G4Box("BoxBase", boxBaseX/2, boxBaseY/2, boxBaseZ/2);
G4Box* rm1Box = new G4Box("Rm1Box", rm1BoxX/2, rm1BoxY, rm1BoxZ);
G4Box* rm2Box = new G4Box("Rm2Box", rm2BoxX, rm2BoxY, rm2BoxZ);
G4Box* rm3Box = new G4Box("Rm3Box", rm3BoxX, rm3BoxY, rm3BoxZ);

G4RotationMatrix* rot0 = new G4RotationMatrix;
G4double posRmBox1Y = 5.334/2*mm;
G4ThreeVector posRmBox1Vect = G4ThreeVector(0,posRmBox1Y,0);
G4Transform3D posRmBox1Tr = G4Transform3D(*rot0, posRmBox1Vect);

G4double posRmBox2Y = -5.334/2*mm;
G4double posRmBox2Z = 11.43/2*mm;
G4ThreeVector posRmBox2Vect = G4ThreeVector(0,posRmBox2Y, posRmBox2Z);
G4Transform3D posRmBox2Tr = G4Transform3D (*rot0, posRmBox2Vect);

G4double posRmBox3X = 116.84/2*mm;
G4double posRmBox3Z = 11.43/2*mm;
G4ThreeVector posRmBox3Vect = G4ThreeVector(posRmBox3X,0,posRmBox3Z);
G4ThreeVector posRmBox4Vect = G4ThreeVector(-posRmBox3X,0,posRmBox3Z);

G4Transform3D posRmBox3Tr = G4Transform3D(*rot0, posRmBox3Vect);
G4Transform3D posRmBox4Tr = G4Transform3D(*rot0, posRmBox4Vect);


G4SubtractionSolid* copperMount1 = new G4SubtractionSolid("CopperMount1", boxBase, rm1Box, posRmBox1Tr);
G4SubtractionSolid* copperMount2 = new G4SubtractionSolid("CopperMount2", copperMount1, rm2Box, posRmBox2Tr);
G4SubtractionSolid* copperMount3 = new G4SubtractionSolid("CopperMount3", copperMount2, rm3Box, posRmBox3Tr);
G4SubtractionSolid* copperMountFinal = new G4SubtractionSolid("CopperMountFinal", copperMount3, rm3Box, posRmBox4Tr);

G4Material* Copper = G4Material::GetMaterial("G4_Cu");
G4LogicalVolume *copperMountBar= new G4LogicalVolume(copperMountFinal, Copper,"CopperMountBar");

return copperMountBar;

}

G4LogicalVolume* GetConstructionProtectiveCover()
{
  G4double mainBoxX = 88.9*mm;
  G4double mainBoxY = 85.852*mm;
  G4double mainBoxZ = 1.524*mm;
  G4double rmBoxXY = 7.62;
  G4double rmBoxZ = 1.524*mm;


  G4Box* mainBox = new G4Box( "MainBox", mainBoxX/2, mainBoxY/2, mainBoxZ/2);
  G4Box* rmBox = new G4Box("RmBox", rmBoxXY/2, rmBoxXY/2, rmBoxZ);

  G4RotationMatrix* rot45 = new G4RotationMatrix;
  rot45->rotateZ(M_PI/4.*rad);
  G4double posX = mainBoxX/2;
  G4double posY = mainBoxY/2;
  G4ThreeVector rmBoxVect1 = G4ThreeVector(posX, posY, 0);
  G4ThreeVector rmBoxVect2 = G4ThreeVector(posX, -posY, 0);
  G4ThreeVector rmBoxVect3 = G4ThreeVector(-posX, -posY, 0);
  G4ThreeVector rmBoxVect4 = G4ThreeVector(-posX, posY, 0);

  G4Transform3D rmBoxTr1 = G4Transform3D(*rot45, rmBoxVect1);
  G4Transform3D rmBoxTr2 = G4Transform3D(*rot45, rmBoxVect2);
  G4Transform3D rmBoxTr3 = G4Transform3D(*rot45, rmBoxVect3);
  G4Transform3D rmBoxTr4 = G4Transform3D(*rot45, rmBoxVect4);


  G4SubtractionSolid* protective1 = new G4SubtractionSolid("Protective1", mainBox, rmBox, rmBoxTr1);
  G4SubtractionSolid* protective2 = new G4SubtractionSolid("Protective2", protective1, rmBox, rmBoxTr2);
  G4SubtractionSolid* protective3 = new G4SubtractionSolid("Protective3", protective2, rmBox, rmBoxTr3);
  G4SubtractionSolid* protectiveCoverFinal = new G4SubtractionSolid("ProtectiveCoverFinal", protective3, rmBox, rmBoxTr4);

  G4Material* Aluminium = G4Material::GetMaterial("G4_Al");
  G4LogicalVolume *protectiveCover= new G4LogicalVolume(protectiveCoverFinal, Aluminium,"ProtectiveCover");

  return protectiveCover;
}

G4LogicalVolume* GetConstructionCCDSensor44(){

  G4double mainBoxX = 63.99*mm;
  G4double mainBoxY = 64.36*mm;
  G4double mainBoxZ = 0.25*mm;

  G4Box* mainBox  = new G4Box("mainBox", mainBoxX/2, mainBoxY/2, mainBoxZ/2);

  G4Material* Silicon = G4Material::GetMaterial("G4_Si");
  G4LogicalVolume *CCDSensor= new G4LogicalVolume(mainBox, Silicon,"CCDSensor");

  return CCDSensor;

}


G4LogicalVolume* GetConstructionSiliconSubstrate44(){

  G4double mainBoxX = 76.20*mm;
  G4double mainBoxY = 88.90*mm;
  G4double mainBoxZ = 0.65*mm;
  G4double rmBoxX = 56.39*mm;
  G4double rmBoxY = 56.64*mm;
  G4double rmBoxZ = 0.65*mm;

  G4Box* mainBox = new G4Box("MainBox", mainBoxX/2, mainBoxY/2, mainBoxZ/2);
  G4Box* rmBox = new G4Box("RmBox", rmBoxX/2, rmBoxY/2, rmBoxZ);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4ThreeVector Vect0 = G4ThreeVector(0,0,0);
  G4Transform3D Tr0 = G4Transform3D(*rot0, Vect0);


  G4SubtractionSolid* siliconSubstrateFinal = new G4SubtractionSolid("SiliconSubstrateFinal", mainBox, rmBox, Tr0 );

  G4Material* Silicon = G4Material::GetMaterial("G4_Si");
  G4LogicalVolume* siliconSubstrate = new G4LogicalVolume(siliconSubstrateFinal, Silicon, "SiliconSubstrate");

  return siliconSubstrate;
}


G4LogicalVolume* GetConstructionCopperTopPlate44(){

  G4double mainBoxX = 85.852*mm;
  G4double mainBoxY = 101.6*mm;
  G4double mainBoxZ = 2.006*mm;
  G4double rmBoxX1 = 67.31*mm;
  G4double rmBoxY1 = 66.700*mm;
  G4double rmBoxZ1 = 2.006*mm;
  G4double rmBoxX2 = 10.033*mm;
  G4double rmBoxY2 = 10.033*mm;
  G4double rmBoxZ2 = 2.006*mm;

  G4Box* mainBox = new G4Box("MainBox", mainBoxX/2, mainBoxY/2, mainBoxZ/2);
  G4Box* rmBox1 = new G4Box("RmBox1", rmBoxX1/2, rmBoxY1/2, rmBoxZ1);
  G4Box* rmBox2 = new G4Box("RmBox2", rmBoxX2/2, rmBoxY2/2, rmBoxZ2);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4RotationMatrix* rot45 = new G4RotationMatrix;
  rot45->rotateZ(M_PI/4.*rad);
  G4double posX = mainBoxX/2;
  G4double posY = mainBoxY/2;
  G4ThreeVector Vect0 = G4ThreeVector(0,0,0);
  G4ThreeVector rmBoxVect1 = G4ThreeVector( posX, posY,0);
  G4ThreeVector rmBoxVect2 = G4ThreeVector( posX, -posY,0);
  G4Transform3D rmBoxTr0 = G4Transform3D(*rot0, Vect0);
  G4Transform3D rmBoxTr1 = G4Transform3D(*rot45, rmBoxVect1);
  G4Transform3D rmBoxTr2 = G4Transform3D(*rot45, rmBoxVect2);

  G4SubtractionSolid * topPlate1 = new G4SubtractionSolid("TopPlate1", mainBox, rmBox1, rmBoxTr0);
  G4SubtractionSolid * topPlate2 = new G4SubtractionSolid("TopPlate2", topPlate1, rmBox2, rmBoxTr1);
  G4SubtractionSolid * copperTopPlateFinal = new G4SubtractionSolid("CopperTopPlateFinal", topPlate2, rmBox2, rmBoxTr2);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume* copperTopPlate = new G4LogicalVolume(copperTopPlateFinal, Copper, "CopperTopPlate");

  return copperTopPlate;
}

G4LogicalVolume* GetConstructionFlexCable44()
{
  G4double longPartFlexX = 556.768*mm;
  G4double longPartFlexY = 25.4*mm;
  G4double longPartFlexZ = 0.305*mm;
  G4double boxMainX = 90.398*mm;
  G4double boxMainY = 101.6*mm;
  G4double boxMainZ = 0.305*mm;
  G4double rmBox1X = 77.698*mm;
  G4double rmBox1Y = 67.411*mm;
  G4double rmBox1Z = 0.305*mm;
  G4double rmBox2X = 12.7*mm;
  G4double rmBox2Y = 82.55*mm;
  G4double rmBox2Z = 0.305*mm;

  G4Box* longPartFlex = new G4Box("LongPartFlex", longPartFlexX/2, longPartFlexY/2, longPartFlexZ/2);
  G4Box* boxMain = new G4Box("BoxMain", boxMainX/2, boxMainY/2, boxMainZ/2);
  G4Box* rmBox1 = new G4Box("RmBox1", rmBox1X, rmBox1Y/2, rmBox1Z);
  G4Box* rmBox2 = new G4Box("RmBox2", rmBox2X/2, rmBox2Y/2, rmBox2Z);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4double posRmBoxX = -boxMainX/2;
  G4ThreeVector posRmVect = G4ThreeVector(posRmBoxX, 0, 0);
  G4Transform3D posRmTr = G4Transform3D(*rot0, posRmVect);

  G4double posRmBox2X = posRmBoxX+rmBox1X-rmBox2X/2;
  G4ThreeVector posRmBox2Vect = G4ThreeVector(posRmBox2X, 0,0);
  G4Transform3D posRmBox2Tr = G4Transform3D(*rot0, posRmBox2Vect);

  G4double posLongPartX = boxMainX/2 + longPartFlexX/2;
  G4ThreeVector posLongPartVect = G4ThreeVector(posLongPartX,0,0);
  G4Transform3D posLongPartTr = G4Transform3D(*rot0, posLongPartVect);

  G4SubtractionSolid* flex1 = new G4SubtractionSolid("Flex1", boxMain, rmBox1, posRmTr);
  G4SubtractionSolid* flex2 = new G4SubtractionSolid("Flex2", flex1, rmBox2, posRmBox2Tr);
  //G4UnionSolid* flexCableFinal = new G4UnionSolid("FlexCableFinal", flex2, longPartFlex, posLongPartTr);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume* flexCable = new G4LogicalVolume(flex2, Copper, "FlexCable");

  return flexCable;
}
