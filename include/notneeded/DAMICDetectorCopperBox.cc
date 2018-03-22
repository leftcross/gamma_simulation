#include "DAMICDetectorCopperBox.hh"
#include "DAMICDetectorModules.hh"
#include "DAMICDetectorInnerLead.hh"
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


G4LogicalVolume* GetConstructionCopperAndInner()
{
  G4double prinCylR = 191.287/2*mm;
  G4double prinCylH = 413.512*mm;
  G4double zero = 0;
  G4double a2pi = M_PI*2+1;

  G4Tubs* prinTubGeo = new G4Tubs("PrinTubGeo",zero, prinCylR, prinCylH/2, zero, a2pi);
  G4Material* Galac = G4Material::GetMaterial("G4_Galactic");

  G4LogicalVolume* PrinTubLV = new G4LogicalVolume(prinTubGeo, Galac, "PrinTubLV");

  G4LogicalVolume* ColdFingerLV = GetConstructionColdFinger();
  G4double PosZColdFinger = prinCylH/2-248.158/2*mm;
  G4ThreeVector VectColdFinger = G4ThreeVector(0,0,PosZColdFinger);


  G4LogicalVolume* TopPlateLV = GetConstructionTopPlate();
  G4double PosZTopPlate = PosZColdFinger-269.748/2*mm-9.652/2*mm;;
  G4ThreeVector VectTopPlate = G4ThreeVector(0,0,PosZTopPlate);
  G4RotationMatrix* rotTopPlate = new G4RotationMatrix;
  rotTopPlate->rotateZ(M_PI/2);

  G4LogicalVolume* RearPlateLV = GetConstructionRearPlate();
  G4double PosZRearPlate = PosZTopPlate - 148.59/2+9.652/2*mm;
  G4double PosYRearPlate = 131.826/2+6.35/2*mm;
  G4double PosXRearPlate = 0;
  G4ThreeVector VectRearPlate = G4ThreeVector(PosXRearPlate,PosYRearPlate,PosZRearPlate);
  G4ThreeVector uRear = G4ThreeVector(0,1,0);
  G4ThreeVector vRear = G4ThreeVector(0,0,1);
  G4ThreeVector wRear = G4ThreeVector(1,0,0);
  G4RotationMatrix* rotRearPlate = new G4RotationMatrix(uRear,vRear,wRear);

  G4LogicalVolume* EndCoverPlateLV = GetConstructionEndCoverPlate();
  G4double PosZEndCoverPlate = PosZRearPlate;
  G4double PosYEndCoverPlate = -131.826/2-3.048/2*mm;
  G4double PosXEndCoverPlate = 0;
  G4ThreeVector VectEndPlate = G4ThreeVector(PosXEndCoverPlate, PosYEndCoverPlate, PosZEndCoverPlate );
  G4ThreeVector uEndP = G4ThreeVector(0,-1,0);
  G4ThreeVector vEndP = G4ThreeVector(0,0,-1);
  G4ThreeVector wEndP = G4ThreeVector(1,0,0);
  G4RotationMatrix* rotEndPlate = new G4RotationMatrix(uEndP,vEndP,wEndP);


  G4LogicalVolume* SidePlate1LV = GetConstructionSidePlate();
  G4double PosZSidePlate1 = PosZRearPlate;
  G4double PosYSidePlate1 = PosYRearPlate-6.35/2*mm-113.03/2*mm;
  G4double PosXSidePlate1 = -56.606*mm;
  G4ThreeVector VectSidePlate1 = G4ThreeVector(PosXSidePlate1, PosYSidePlate1, PosZSidePlate1);
  G4ThreeVector uSide1 = G4ThreeVector(0,0,1);
  G4ThreeVector vSide1 = G4ThreeVector(0,1,0);
  G4ThreeVector wSide1 = G4ThreeVector(-1,0,0);
  G4RotationMatrix* rotSidePlate1 = new G4RotationMatrix(uSide1, vSide1, wSide1);

  G4LogicalVolume* SidePlate2LV = GetConstructionSidePlate();
  G4double PosZSidePlate2 = PosZRearPlate;
  G4double PosYSidePlate2 = PosYSidePlate1;
  G4double PosXSidePlate2 = -PosXSidePlate1;
  G4ThreeVector VectSidePlate2 = G4ThreeVector(PosXSidePlate2, PosYSidePlate2, PosZSidePlate2);
  G4ThreeVector uSide2 = G4ThreeVector(0,0,-1);
  G4ThreeVector vSide2 = G4ThreeVector(0,1,0);
  G4ThreeVector wSide2 = G4ThreeVector(1,0,0);
  G4RotationMatrix* rotSidePlate2 = new G4RotationMatrix(uSide2, vSide2, wSide2);

  G4LogicalVolume* EndCoverBar1LV = GetConstructionEndCoverBar();
  G4double PosZEndCoverBar1 = PosZRearPlate;
  G4double PosYEndCoverBar1 = -56.7055*mm;
  G4double PosXEndCoverBar1 = -61.570*mm;
  G4ThreeVector VectEndCoverBar1 = G4ThreeVector(PosXEndCoverBar1, PosYEndCoverBar1, PosZEndCoverBar1);
  G4ThreeVector uBar1 = G4ThreeVector(0,1,0);
  G4ThreeVector vBar1 = G4ThreeVector(0,0,1);
  G4ThreeVector wBar1 = G4ThreeVector(1,0,0);
  G4RotationMatrix* rotEndCoverBar1 = new G4RotationMatrix(uBar1, vBar1, wBar1);

  G4LogicalVolume* EndCoverBar2LV = GetConstructionEndCoverBar();
  G4double PosZEndCoverBar2 = PosZRearPlate;
  G4double PosYEndCoverBar2 = PosYEndCoverBar1;
  G4double PosXEndCoverBar2 = -PosXEndCoverBar1;
  G4ThreeVector VectEndCoverBar2 = G4ThreeVector(PosXEndCoverBar2, PosYEndCoverBar2, PosZEndCoverBar2);
  G4ThreeVector uBar2 = G4ThreeVector(0,-1,0);
  G4ThreeVector vBar2 = G4ThreeVector(0,0,1);
  G4ThreeVector wBar2 = G4ThreeVector(-1,0,0);
  G4RotationMatrix* rotEndCoverBar2 = new G4RotationMatrix(uBar2, vBar2, wBar2);

  G4LogicalVolume* BottomPlateLV = GetConstructionBottomPlate();
  G4double PosZBottomPlate = PosZRearPlate-148.59/2+6.35*mm;
  G4double PosYBottomPlate = 0;
  G4double PosXBottomPlate = 0;
  G4ThreeVector VectBottomPlate = G4ThreeVector(PosXBottomPlate, PosYBottomPlate, PosZBottomPlate);
  G4ThreeVector uBottom = G4ThreeVector(0,1,0);
  G4ThreeVector vBottom = G4ThreeVector(-1,0,0);
  G4ThreeVector wBottom = G4ThreeVector(0,0,1);
  G4RotationMatrix* rotBottomPlate = new G4RotationMatrix(uBottom, vBottom, wBottom);

  G4LogicalVolume* Insulator1LV = GetConstructionInsulator();
  G4double PosZInsulator1 = PosZSidePlate1 + 89.916*mm;
  G4double PosYInsulator1 = 72.517*mm;
  G4double PosXInsulator1 = 0;
  G4ThreeVector VectInsulator1 = G4ThreeVector(PosXInsulator1, PosYInsulator1, PosZInsulator1);
  G4ThreeVector uInsu1 = G4ThreeVector(1,0,0);
  G4ThreeVector vInsu1 = G4ThreeVector(0,0,1);
  G4ThreeVector wInsu1 = G4ThreeVector(0,-1,0);
  G4RotationMatrix* rotInsulator1 = new G4RotationMatrix(uInsu1, vInsu1, wInsu1);

  G4LogicalVolume* Insulator2LV = GetConstructionInsulator();
  G4double PosZInsulator2 = PosZInsulator1;
  G4double PosYInsulator2 = PosYSidePlate1-24.765*mm;
  G4double PosXInsulator2 = -60.734*mm;
  G4ThreeVector VectInsulator2 = G4ThreeVector(PosXInsulator2, PosYInsulator2, PosZInsulator2);
  G4ThreeVector uInsu2 = G4ThreeVector(0,0,-1);
  G4ThreeVector vInsu2 = G4ThreeVector(1,0,0);
  G4ThreeVector wInsu2 = G4ThreeVector(0,-1,0);
  G4RotationMatrix* rotInsulator2 = new G4RotationMatrix(uInsu2, vInsu2, wInsu2);


  G4LogicalVolume* Insulator3LV = GetConstructionInsulator();
  G4double PosZInsulator3 = PosZInsulator1;
  G4double PosYInsulator3 = PosYInsulator2;
  G4double PosXInsulator3 = -PosXInsulator2;
  G4ThreeVector VectInsulator3 = G4ThreeVector(PosXInsulator3, PosYInsulator3, PosZInsulator3);
  G4ThreeVector uInsu3 = G4ThreeVector(0,0,1);
  G4ThreeVector vInsu3 = G4ThreeVector(-1,0,0);
  G4ThreeVector wInsu3 = G4ThreeVector(0,-1,0);
  G4RotationMatrix* rotInsulator3 = new G4RotationMatrix(uInsu3, vInsu3, wInsu3);

  G4LogicalVolume* UpperPlateLV = GetConstructionUpperPlate();
  G4double PosZUpperPlate = PosZInsulator1 + 18.923*mm+ 0.127*mm;
  G4double PosYUpperPlate = PosYInsulator1-0.508/2*mm-139.7/2*mm;
  G4double PosXUpperPlate = 0 ;
  G4ThreeVector VectUpperPlate = G4ThreeVector(PosXUpperPlate, PosYUpperPlate, PosZUpperPlate);
  G4ThreeVector uUpper = G4ThreeVector(0,1,0);
  G4ThreeVector vUpper = G4ThreeVector(-1,0,0);
  G4ThreeVector wUpper = G4ThreeVector(0,0,1);
  G4RotationMatrix* rotUpperPlate = new G4RotationMatrix(uUpper, vUpper, wUpper);


  G4LogicalVolume* ColdFingerSpacer1LV = GetConstructionColdFingerSpacer();
  G4double PosZColdFingerS1 = PosZColdFinger + 177.8/2+1.016/2;
  G4double PosYColdFingerS1 = 0;
  G4double PosXColdFingerS1 = 0;
  G4ThreeVector VectColdFingerS1 = G4ThreeVector(PosXColdFingerS1, PosYColdFingerS1, PosZColdFingerS1);

  G4LogicalVolume* ColdFingerSpacer2LV = GetConstructionColdFingerSpacer();
  G4double PosZColdFingerS2 = PosZColdFinger - (177.8/2+1.016/2);
  G4double PosYColdFingerS2 = 0;
  G4double PosXColdFingerS2 = 0;
  G4ThreeVector VectColdFingerS2 = G4ThreeVector(PosXColdFingerS2, PosYColdFingerS2, PosZColdFingerS2);

  /*-------------------------------------------MODULES---------------------------------------*/
  G4LogicalVolume* ModuleLV1 = GetConstructionModule44();
  G4double PosZModule44 = PosZBottomPlate + 6.35/2*mm + 5.994/2*mm+3*mm;
  G4double PosYModule44 = PosYRearPlate - 6.35/2*mm - 116.332/2*mm -1.800*mm;
  G4ThreeVector Module44Vect = G4ThreeVector(0,PosYModule44,PosZModule44);
  G4ThreeVector uModule44 = G4ThreeVector(0,-1,0);
  G4ThreeVector vModule44 = G4ThreeVector(1,0,0);
  G4ThreeVector wModule44 = G4ThreeVector(0,0,1);
  G4RotationMatrix* Module44Rot = new G4RotationMatrix(uModule44, vModule44, wModule44);

  /*-------------------------------------------INNERLEAD------------------------------------------*/

  G4LogicalVolume* BoxMountingPlateLV = GetConstructionBoxMountingPlate();
  G4double BoxMountingPlateZ = PosZUpperPlate + 6.35/2*mm +9.652*mm/2;
  G4double BoxMountingPlateY = 0;
  G4double BoxMountingPlateX = 0;
  G4ThreeVector VectBoxMountingPlate = G4ThreeVector(BoxMountingPlateX, BoxMountingPlateY, BoxMountingPlateZ);
  G4ThreeVector uPlate = G4ThreeVector(0,1,0);
  G4ThreeVector vPlate = G4ThreeVector(-1,0,0);
  G4ThreeVector wPlate = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotBoxMountingPlate = new G4RotationMatrix(uPlate,vPlate,wPlate);

  G4LogicalVolume* AncientLead2LV = GetConstructionAncientLead2();
  G4double AncientLead2Z = BoxMountingPlateZ+9.652*mm/2*mm+25.4/2*mm;
  G4double AncientLead2Y = 0;
  G4double AncientLead2X = 0;
  G4double COSSIN3P4 = sqrt(2)/2;
  G4ThreeVector VectAncientLead2 = G4ThreeVector(AncientLead2X, AncientLead2Y, AncientLead2Z);
  G4ThreeVector uAncient2 = G4ThreeVector(-COSSIN3P4,COSSIN3P4,0);
  G4ThreeVector vAncient2 = G4ThreeVector(-COSSIN3P4,-COSSIN3P4,0);
  G4ThreeVector wAncient2 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotAncientLead2 = new G4RotationMatrix(uAncient2,vAncient2,wAncient2);

  G4LogicalVolume* AncientLead1LV = GetConstructionAncientLead1();
  G4double AncientLead1Z = AncientLead2Z;
  G4double AncientLead1Y = 0;
  G4double AncientLead1X = 0;
  G4ThreeVector VectAncientLead1 = G4ThreeVector(AncientLead1X, AncientLead1Y, AncientLead1Z);
  G4ThreeVector uAncient1 = G4ThreeVector(COSSIN3P4,-COSSIN3P4,0);
  G4ThreeVector vAncient1 = G4ThreeVector(COSSIN3P4,COSSIN3P4,0);
  G4ThreeVector wAncient1 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotAncientLead1 = new G4RotationMatrix(uAncient1,vAncient1,wAncient1);

  G4LogicalVolume* GoodLead3LV = GetConstructionGoodLead();
  G4double PosGoodLead3Z = AncientLead2Z+25.4/2*mm+50.8/2*mm;
  G4double PosGoodLead3Y = 0;
  G4double PosGoodLead3X = 0;
  G4ThreeVector VectGoodLead3 = G4ThreeVector(PosGoodLead3X, PosGoodLead3Y, PosGoodLead3Z);
  G4ThreeVector uGoodLead3 = G4ThreeVector(0,-1,0);
  G4ThreeVector vGoodLead3 = G4ThreeVector(1,0,0);
  G4ThreeVector wGoodLead3 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotGoodLead3 = new G4RotationMatrix(uGoodLead3, vGoodLead3, wGoodLead3);

  G4LogicalVolume* SpacerPlate2LV = GetConstructionSpacerPlate();
  G4double PosSpacerPlate2Z = PosGoodLead3Z+3.048/2*mm+50.8/2*mm;
  G4double PosSpacerPlate2Y = 0;
  G4double PosSpacerPlate2X = 0;
  G4ThreeVector VectSpacerPlate2 = G4ThreeVector(PosSpacerPlate2X, PosSpacerPlate2Y, PosSpacerPlate2Z);
  G4ThreeVector uSpacerPlate2 = G4ThreeVector(0,-1,0);
  G4ThreeVector vSpacerPlate2 = G4ThreeVector(1,0,0);
  G4ThreeVector wSpacerPlate2 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotSpacerPlate2 = new G4RotationMatrix(uSpacerPlate2, vSpacerPlate2, wSpacerPlate2);

  G4LogicalVolume* GoodLead2LV = GetConstructionGoodLead();
  G4double PosGoodLead2Z = PosSpacerPlate2Z+50.8/2*mm+3.048/2*mm;
  G4double PosGoodLead2Y = 0;
  G4double PosGoodLead2X = 0;
  G4ThreeVector VectGoodLead2 = G4ThreeVector(PosGoodLead2X, PosGoodLead2Y, PosGoodLead2Z);
  G4ThreeVector uGoodLead2 = G4ThreeVector(0,-1,0);
  G4ThreeVector vGoodLead2 = G4ThreeVector(1,0,0);
  G4ThreeVector wGoodLead2 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotGoodLead2 = new G4RotationMatrix(uGoodLead2, vGoodLead2, wGoodLead2);

  G4LogicalVolume* GoodLead1LV = GetConstructionGoodLead();
  G4double PosGoodLead1Z = PosGoodLead2Z+50.8*mm;
  G4double PosGoodLead1Y = 0;
  G4double PosGoodLead1X = 0;
  G4ThreeVector VectGoodLead1 = G4ThreeVector(PosGoodLead1X, PosGoodLead1Y, PosGoodLead1Z);
  G4ThreeVector uGoodLead1 = G4ThreeVector(0,-1,0);
  G4ThreeVector vGoodLead1 = G4ThreeVector(1,0,0);
  G4ThreeVector wGoodLead1 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotGoodLead1 = new G4RotationMatrix(uGoodLead1, vGoodLead1, wGoodLead1);

  G4LogicalVolume* SpacerPlate1LV = GetConstructionSpacerPlate();
  G4double PosSpacerPlate1Z = PosGoodLead1Z+50.8/2*mm+3.048/2*mm;
  G4double PosSpacerPlate1Y = 0;
  G4double PosSpacerPlate1X = 0;
  G4ThreeVector VectSpacerPlate1 = G4ThreeVector(PosSpacerPlate1X, PosSpacerPlate1Y, PosSpacerPlate1Z);
  G4ThreeVector uSpacerPlate1 = G4ThreeVector(0,-1,0);
  G4ThreeVector vSpacerPlate1 = G4ThreeVector(1,0,0);
  G4ThreeVector wSpacerPlate1 = G4ThreeVector(0,0,1);
  G4RotationMatrix* RotSpacerPlate1 = new G4RotationMatrix(uSpacerPlate1, vSpacerPlate1, wSpacerPlate1);



  G4PVPlacement* ColdFingerPV = new G4PVPlacement(0, VectColdFinger, ColdFingerLV, "ColdFingerPV", PrinTubLV, false, 0, false);
  G4PVPlacement* TopPlatePV = new G4PVPlacement(rotTopPlate, VectTopPlate, TopPlateLV, "TopPlatePV", PrinTubLV, false, 0, false);
  G4PVPlacement* RearPlatePV = new G4PVPlacement(rotRearPlate, VectRearPlate, RearPlateLV, "RearPlatePV",PrinTubLV, false, 0, false);
  G4PVPlacement* EndCoverPlatePV = new G4PVPlacement(rotEndPlate,VectEndPlate, EndCoverPlateLV, "EndCoverPlatePV", PrinTubLV, false, 0, false);
  G4PVPlacement* SidePlate1PV = new G4PVPlacement(rotSidePlate1, VectSidePlate1, SidePlate1LV, "SidePlate1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* SidePlate2PV = new G4PVPlacement(rotSidePlate2, VectSidePlate2, SidePlate2LV, "SidePlate2PV", PrinTubLV, false, 0, false);
  G4PVPlacement* EndCoverBar1PV = new G4PVPlacement(rotEndCoverBar1, VectEndCoverBar1, EndCoverBar1LV, "EndCoverBar1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* EndCoverBar2PV = new G4PVPlacement(rotEndCoverBar2, VectEndCoverBar2, EndCoverBar2LV, "EndCoverBar2PV", PrinTubLV, false, 0, false);
  G4PVPlacement* BottomPlatePV = new G4PVPlacement(rotBottomPlate, VectBottomPlate, BottomPlateLV, "BottomPlatePV", PrinTubLV, false, 0, false);
  G4PVPlacement* Insulator1PV = new G4PVPlacement(rotInsulator1, VectInsulator1 , Insulator1LV, "Insulator1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* Insulator2PV = new G4PVPlacement(rotInsulator2, VectInsulator2, Insulator2LV, "Insulator2PV", PrinTubLV, false, 0, false);
  G4PVPlacement* Insulator3PV = new G4PVPlacement(rotInsulator3, VectInsulator3, Insulator3LV, "Insulator3PV", PrinTubLV, false, 0, false);
  G4PVPlacement* UpperPlatePV = new G4PVPlacement(rotUpperPlate, VectUpperPlate, UpperPlateLV, "UpperPlatePV", PrinTubLV, false, 0, false);
  G4PVPlacement* ColdFingerSpacer1PV = new G4PVPlacement(0, VectColdFingerS1, ColdFingerSpacer1LV, "ColdFingerSpacer1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* ColdFingerSpacer2PV = new G4PVPlacement(0, VectColdFingerS2, ColdFingerSpacer2LV, "ColdFingerSpacer2PV", PrinTubLV, false, 0, false);


  G4PVPlacement* ModulePV1 = new G4PVPlacement(Module44Rot, Module44Vect, ModuleLV1, "ModulePV1", PrinTubLV, false, 0, false);

  G4PVPlacement* SpacerPlate1PV = new G4PVPlacement(RotSpacerPlate1, VectSpacerPlate1, SpacerPlate1LV, "SpacerPlate1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* GoodLead1PV = new G4PVPlacement(RotGoodLead1, VectGoodLead1, GoodLead1LV, "GoodLead1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* GoodLead2PV = new G4PVPlacement(RotGoodLead2, VectGoodLead2, GoodLead2LV, "GoodLead2PV", PrinTubLV, false, 0, false);
  G4PVPlacement* SpacerPlate2PV = new G4PVPlacement(RotSpacerPlate2, VectSpacerPlate2, SpacerPlate2LV, "SpacerPlate2PV", PrinTubLV, false, 0, false);
  G4PVPlacement* GoodLead3PV = new G4PVPlacement(RotGoodLead3, VectGoodLead3, GoodLead3LV, "GoodLead3PV", PrinTubLV, false, 0, false);
  G4PVPlacement* AncientLead1PV = new G4PVPlacement(RotAncientLead1, VectAncientLead1, AncientLead1LV, "AncientLead1PV", PrinTubLV, false, 0, false);
  G4PVPlacement* AncientLead2PV = new G4PVPlacement(RotAncientLead2, VectAncientLead2, AncientLead2LV, "AncientLead2PV", PrinTubLV, false, 0, false);
  G4PVPlacement* BoxMountingPlatePV = new G4PVPlacement(RotBoxMountingPlate, VectBoxMountingPlate, BoxMountingPlateLV, "BoxMountingPlatePV", PrinTubLV, false, 0, false);

  return PrinTubLV;
}

G4LogicalVolume* GetConstructionBottomPlate()
{
  G4double longueurX = 131.826*mm;
  G4double longueurY = 100.0252*mm;
  G4double largeurZ = 6.35*mm;
  G4double vis1Rayon2 = 0;
  G4double vis2Rayon6 = 0;
  G4double vis1Longueur2 = 0;
  G4double vis2Longueur6 = 0;
  G4double espaceEntreVisY = 14.605*mm;
  G4double espaceEntreVisX1 = 28.224*mm;
  G4double espaceEntreVisX2 = 9.6012*mm;
  G4Box* bottomBox = new G4Box("BottomPlate", longueurX/2, longueurY/2, largeurZ/2);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume* BottomPlateLV = new G4LogicalVolume(bottomBox, Copper,"BottomPlateLV");

  return BottomPlateLV;

}

G4LogicalVolume* GetConstructionColdFinger(){
  G4double petitsCylindresR = 22.86/2*mm;
  G4double petitsCylindresH = 1.016*mm;
  G4double grandCylindreR = 25.4/2*mm;
  G4double grandCylindreH = 177.8*mm;
  G4double boutCylindre1H = 34.163*mm;
  G4double boutCylindre2H = 44.958*mm;
  G4double boutCoupe1LX = 3.175*mm;
  G4double boutCoupe1LY = 100*mm;
  G4double boutCoupe1LZ = 25.4*mm;
  G4double boutCoupe2LX = 2.54*mm;
  G4double boutCoupe2LY = 100*mm;
  G4double boutCoupe2LZ = 22.86*mm;
  G4double trouDiam = 7.1374*mm;
  G4double vis1Longueur2 =0;
  G4double vis1Rayon2 =0;
  G4double vis2Longueur2 = 0;
  G4double vis2Rayon2 = 0;

  G4double zeros = 0;
  G4double angle360 = 360;
  G4Tubs* middlePart = new G4Tubs( "MiddlePart",zeros, grandCylindreR, grandCylindreH/2, zeros, angle360);
  G4Tubs* littleTubs1 = new G4Tubs ( "LittleTub1", zeros, petitsCylindresR, petitsCylindresH/2,zeros, angle360);
  G4Tubs* littleTubs2 = new G4Tubs ("LittleTub2", zeros, petitsCylindresR, petitsCylindresH/2,zeros,angle360);

  G4Tubs* boutTubs1 = new G4Tubs ("boutTubs1", zeros, grandCylindreR, boutCylindre1H/2, zeros, angle360);
  G4Tubs* boutTubs2 = new G4Tubs ("boutTubs2", zeros, grandCylindreR, boutCylindre2H/2, zeros, angle360);
  G4Tubs* trouTubs = new G4Tubs ("trouTubs", zeros, trouDiam/2, boutCoupe1LZ, zeros, angle360);

  G4Box* boutBox1 = new G4Box ("boutBox1", boutCoupe1LX, boutCoupe1LY/2, boutCoupe1LZ);
  G4Box* boutBox2 = new G4Box ("boutBox2",boutCoupe2LX, boutCoupe2LY/2, boutCoupe2LZ);

  G4RotationMatrix * zRot = new G4RotationMatrix;
  zRot->rotateZ(M_PI/4.*rad);
  G4RotationMatrix * zzRot = new G4RotationMatrix;
  zzRot->rotateZ(0*rad);
  G4RotationMatrix *zRotbis = new G4RotationMatrix;
  zRotbis->rotateZ(3*M_PI/4*rad);
  G4RotationMatrix * rotHole = new G4RotationMatrix;
  rotHole->rotateY(M_PI/2*rad);
  rotHole->rotateZ(M_PI/4*rad);

  G4double bout1Trans1X = grandCylindreR;
  G4double bout1Trans1Z = boutCylindre1H/2;
  G4double bout1PosHole = (boutCylindre1H-boutCoupe1LZ)/2;
  G4ThreeVector bout1Trans1 = G4ThreeVector(bout1Trans1X*sqrt(2)/2, sqrt(2)*bout1Trans1X/2 ,bout1Trans1Z);
  G4ThreeVector bout1Trans2 = G4ThreeVector(-bout1Trans1X*sqrt(2)/2, -sqrt(2)*bout1Trans1X/2 ,bout1Trans1Z);
  G4ThreeVector bout1TransHole = G4ThreeVector(0,0,bout1PosHole);

  G4Transform3D boutTransform1 = G4Transform3D (*zRot, bout1Trans1);
  G4Transform3D boutTransform2 = G4Transform3D (*zRot, bout1Trans2);
  G4Transform3D boutHole = G4Transform3D(*rotHole, bout1TransHole);

  G4SubtractionSolid * bout1Sub = new G4SubtractionSolid("Bout1Sub", boutTubs1, boutBox1, boutTransform1);
  G4SubtractionSolid * bout1Sub2 = new G4SubtractionSolid("bout1Sub2", bout1Sub, boutBox1, boutTransform2);
  G4SubtractionSolid * bout1Final = new G4SubtractionSolid("bout1Final", bout1Sub2, trouTubs, boutHole);

  G4RotationMatrix * zRot2 = new G4RotationMatrix;
  zRot2->rotateZ(M_PI/2*rad);
  G4RotationMatrix* zRot0 = new G4RotationMatrix;
  zRot0->rotateZ(0*rad);
  G4double bout2TransX = grandCylindreR;
  G4double bout2TransZ = boutCylindre2H/2;
  G4ThreeVector bout2Trans = G4ThreeVector(bout2TransX, 0, bout2TransZ);
  G4ThreeVector bout2Trans2 = G4ThreeVector(0, bout2TransX, bout2TransZ);

  G4Transform3D bout2SubTr = G4Transform3D(*zRot0,bout2Trans);
  G4Transform3D bout2SubTr1 = G4Transform3D(*zRot2, bout2Trans2);
  G4SubtractionSolid * bout2Sub = new G4SubtractionSolid("Bout2Sub", boutTubs2, boutBox2, bout2SubTr);
  G4SubtractionSolid * bout2Final = new G4SubtractionSolid("bout2Final", bout2Sub, boutBox2, bout2SubTr1 );

  G4double positionLittleCylZ = grandCylindreH/2+petitsCylindresH/2;
  G4double positionBout1Z = grandCylindreH/2+petitsCylindresH+boutCylindre1H/2;
  G4double positionBout2Z = -grandCylindreH/2-petitsCylindresH-boutCylindre2H/2;
  G4RotationMatrix rotm = G4RotationMatrix();
  G4RotationMatrix* rotY = new G4RotationMatrix;
  rotY->rotateY(M_PI*rad);
  rotY->rotateZ(M_PI*rad);
  G4ThreeVector positionMiddlePart = G4ThreeVector(0,0,0);
  G4ThreeVector positionLittleCyl = G4ThreeVector(0,0,positionLittleCylZ);
  G4ThreeVector positionBout1 = G4ThreeVector(0,0,positionBout1Z);
  G4ThreeVector positionBout2 = G4ThreeVector(0,0,positionBout2Z);

  G4Transform3D trMiddle = G4Transform3D(rotm, positionMiddlePart);
  G4Transform3D trLittleCyl1 = G4Transform3D(rotm, positionLittleCyl);
  G4Transform3D trLittleCyl2 = G4Transform3D(rotm, -positionLittleCyl);
  G4Transform3D trBout1 = G4Transform3D(rotm, positionBout1);
  G4Transform3D trBout2 = G4Transform3D(*rotY, positionBout2);

  G4UnionSolid* coldFinger1 = new G4UnionSolid("first", middlePart, littleTubs1, trLittleCyl1);
  G4UnionSolid* coldFinger2 = new G4UnionSolid("second", coldFinger1, littleTubs2, trLittleCyl2);
  G4UnionSolid* coldFinger3 = new G4UnionSolid("third", coldFinger2, bout1Final, trBout1);
  G4UnionSolid* coldFingerFinal = new G4UnionSolid("final", coldFinger3, bout2Final, trBout2);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *ColdFingerLV = new G4LogicalVolume(coldFingerFinal, Copper,"ColdFingerLV");

  return ColdFingerLV;
}
G4LogicalVolume* GetConstructionColdFingerSpacer()
{
  G4double bigRadius = 25.4/2*mm;
  G4double inRadius = 22.86/2*mm;
  G4double bout1Lenght = 2.032*mm;
  G4double allRadius = 27.53/2*mm;
  G4double CFSH = 0.381*mm;
  G4double angle360 = 360;
  G4double angle0 = 0;

  G4double spaceBtwRadius = allRadius - bigRadius;
  G4Tubs * SpacerTub = new G4Tubs("SpaceTub", inRadius, bigRadius, CFSH/2, angle0, angle360);
  G4Box * SpacerBoxX = new G4Box("SpacerBox", spaceBtwRadius, bout1Lenght/2, CFSH/2);
  G4Box * SpacerBoxY = new G4Box("SpacerBox", bout1Lenght/2, spaceBtwRadius, CFSH/2);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4double Trans = allRadius-bout1Lenght/2;
  G4ThreeVector XPosition = G4ThreeVector(Trans,0,0);
  G4ThreeVector XPositionM = G4ThreeVector(-Trans,0,0);
  G4ThreeVector YPosition = G4ThreeVector(0,Trans,0);
  G4ThreeVector YPositionM = G4ThreeVector(0,-Trans,0);

  G4Transform3D trans1 = G4Transform3D(*rot0, XPosition);
  G4Transform3D trans2 = G4Transform3D(*rot0, XPositionM);
  G4Transform3D trans3 = G4Transform3D(*rot0, YPosition);
  G4Transform3D trans4 = G4Transform3D(*rot0, YPositionM);


  G4UnionSolid* spacer1 = new G4UnionSolid("Spacer1", SpacerTub, SpacerBoxX, trans1 );
  G4UnionSolid* spacer2 = new G4UnionSolid("Spacer2", spacer1, SpacerBoxX, trans2 );
  G4UnionSolid* spacer3 = new G4UnionSolid("Spacer3", spacer2, SpacerBoxY, trans3 );
  G4UnionSolid* spacerFinal = new G4UnionSolid("SpacerFinal", spacer3, SpacerBoxY, trans4 );
  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *coldFingerSpacer = new G4LogicalVolume(spacerFinal, Copper,"SpacerFinal");

  return coldFingerSpacer;
}

G4LogicalVolume* GetConstructionEndCoverBar()
{
  G4double longSideBig = 148.59*mm;
  G4double longSideLittle = 4.826*mm;
  G4double littleSideBig = 13.97*mm;
  G4double littleSideLittle = 6.35*mm;
  G4double depth = 18.415*mm;

  G4Box* longSide = new G4Box("longSide", longSideBig/2, longSideLittle/2, depth/2);
  G4Box* littleSide = new G4Box("littleSide", littleSideLittle/2, littleSideBig/2, depth/2);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4double calcTrX = longSideBig/2-littleSideLittle/2;
  G4double calcTrY = (littleSideBig-longSideLittle)/2;
  G4ThreeVector rightBox = G4ThreeVector(calcTrX, calcTrY, 0);
  G4ThreeVector leftBox = G4ThreeVector(-calcTrX, calcTrY, 0);

  G4Transform3D boxRTr = G4Transform3D(*rot0, rightBox);
  G4Transform3D boxLTr = G4Transform3D(*rot0, leftBox);

  G4UnionSolid* endCover1 = new G4UnionSolid("EndCover1", longSide, littleSide, boxRTr);
  G4UnionSolid* endCoverFinal = new G4UnionSolid("EndCoverFinal", endCover1, littleSide, boxLTr);
  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *EndCoverBar = new G4LogicalVolume(endCoverFinal, Copper,"endCoverBar");

  return EndCoverBar;

}

G4LogicalVolume*  GetConstructionEndCoverPlate()
{
  G4double angleTrap = 45; /*in degrees*/
  G4double angle0 = 0;
  G4double depth = 3.048;
  G4double XSide = 148.59;
  G4double YSideDown = 129.032;
  G4double YSideUp = YSideDown-depth*2;
  G4double XBox = 12.7;
  G4double YBox = 31.75;

  G4Trd* EndCoverPlateTrp = new G4Trd("EndCoverPlateTrp",
             XSide/2,
             XSide/2,
             YSideDown/2,
             YSideUp/2,
             depth/2);
  G4Box* boxPlateRM = new G4Box("BoxPlateRM", XBox, YBox/2, depth+1);

  G4RotationMatrix* rot0 = new G4RotationMatrix;

  G4ThreeVector fBox = G4ThreeVector(XSide/2, 0, 0);

  G4Transform3D boxRTr = G4Transform3D(*rot0, fBox);
  G4SubtractionSolid* EndCoverPlateFinal = new G4SubtractionSolid("EndCoverPlateFinal", EndCoverPlateTrp, boxPlateRM, boxRTr);
/*  G4Trap* EndCoverPlateTrap = new G4Trap("EndCoverPlateTrap",
             depth, angle0,
             angle0, YSideDown/2,
             XSide/2, XSide/2,
             angle0, YSideDown/2,
             XSide/2, XSide/2,
             angle0);*/

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *EndCoverPlate = new G4LogicalVolume(EndCoverPlateFinal, Copper,"endCoverPlate");

  return EndCoverPlate;

}

G4LogicalVolume* GetConstructionInsulator()
{
  G4double XBigPlate = 63.5*mm;
  G4double YBigPlate = 44.45*mm;
  G4double XLittleHole = 57.15*mm;
  G4double YLittleHole = 32.766*mm;
  G4double inDepth = 0.508*mm;
  G4double sideDepth = 1.5748*mm;

  G4double YSideBox = (YBigPlate - YLittleHole)/2;
  G4Box * sideBox = new G4Box("SideBox", XBigPlate/2,YSideBox/2, sideDepth/2);
  G4Box * inBox = new G4Box ("InBox", XBigPlate/2, YBigPlate/2, inDepth/2);
  G4Box * holeBox = new G4Box ("HoleBox", XLittleHole/2, YLittleHole/2, sideDepth/2);

  G4RotationMatrix * rot0 = new G4RotationMatrix;
  G4ThreeVector sideBoxV1 = G4ThreeVector(0, YLittleHole/2+YSideBox/2, sideDepth/2-inDepth/2);
  G4ThreeVector sideBoxV2 = G4ThreeVector(0, -YLittleHole/2-YSideBox/2, sideDepth/2-inDepth/2);
  G4ThreeVector fullInsulator = G4ThreeVector(0,0,0);

  G4Transform3D sideBoxTr = G4Transform3D(*rot0, sideBoxV1);
  G4Transform3D sideBoxTr2 = G4Transform3D(*rot0, sideBoxV2);
  G4Transform3D fullInsulatorTr = G4Transform3D(*rot0, fullInsulator);

  G4UnionSolid * insulator1 = new G4UnionSolid("Insulator1",inBox, sideBox, sideBoxTr);
  G4UnionSolid * insulator2 = new G4UnionSolid("Insulator2", insulator1, sideBox, sideBoxTr2);
  G4SubtractionSolid* insulatorFinal = new G4SubtractionSolid("InsulatorFinal", insulator2, holeBox, fullInsulatorTr);

  G4Material* StainSteel = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  G4LogicalVolume * Insulator = new G4LogicalVolume(insulatorFinal, StainSteel,"Insulator");

  return Insulator;
}

G4LogicalVolume*  GetConstructionRearPlate()
{
  G4double angleTrap = 45; /*in degrees*/
  G4double angle0 = 0;
  G4double depth = 6.35;
  G4double XSide = 148.59;
  G4double YSideDown = 119.126;
  G4double YSideUp = YSideDown-depth*2;

  G4Trd* rearPlateFinal = new G4Trd("RearPlateTrd",
             XSide/2,
             XSide/2,
             YSideDown/2,
             YSideUp/2,
             depth/2);

   G4Material* Copper = G4Material::GetMaterial("G4_Cu");
   G4LogicalVolume *rearPlate = new G4LogicalVolume(rearPlateFinal, Copper,"RearPlate");

   return rearPlate;



}

G4LogicalVolume* GetConstructionUpperPlate()
{
  G4double XSide = 139.7*mm;
  G4double YSide = 119.126*mm;
  G4double XHole = 101.6*mm;
  G4double YHole = 81.026*mm;
  G4double depth = 6.35*mm;

  G4Box * upperPlateBox = new G4Box("upperPlateBox", XSide/2, YSide/2, depth/2);
  G4Box * holeBox = new G4Box("holeBox", XHole/2, YHole/2, depth*2);

  G4RotationMatrix * rot0 = new G4RotationMatrix;
  G4ThreeVector V0 = G4ThreeVector(0,0,0);
  G4Transform3D Tr0 = G4Transform3D(*rot0,V0);

  G4SubtractionSolid* upperPlateFinal = new G4SubtractionSolid("UpperPlateFinal", upperPlateBox, holeBox,Tr0 );

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *upperPlate = new G4LogicalVolume(upperPlateFinal, Copper,"UpperPlate");

  return upperPlate;
}

G4LogicalVolume* GetConstructionSidePlate(){

  G4double bigBoxSideX = 148.59*mm;
  G4double bigBoxSideY = 113.03*mm;
  G4double bigBoxSideZ = 7.747*mm;
  G4double middleBoxZ = 9.525*mm;
  G4double middleBoxX = 10.668*mm;
  G4double middleBoxY = 113.03*mm;
  G4double littleBoxX = 1.016*mm;
  G4double littleBoxY = 113.03*mm;
  G4double littleBoxZ = 9.525*mm;

  G4Box * bigBox  = new G4Box(" BigBox", bigBoxSideX/2, bigBoxSideY/2, bigBoxSideZ/2);
  G4Box * middleBox = new G4Box("MiddleBox", middleBoxX/2, middleBoxY/2, middleBoxZ/2);
  G4Box * littleBox = new G4Box("LittleBox", littleBoxX/2, littleBoxY/2, littleBoxZ/2);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4double UnionX = bigBoxSideX/2-middleBoxX/2;
  G4double UnionZ = middleBoxZ - bigBoxSideZ;
  G4ThreeVector Vec1 = G4ThreeVector(UnionX, 0, UnionZ);
  G4ThreeVector Vec2 = G4ThreeVector(-UnionX, 0, UnionZ);
  G4Transform3D trUn = G4Transform3D(*rot0, Vec1);
  G4Transform3D trUn2 = G4Transform3D(*rot0, Vec2);

  G4UnionSolid * bigMidUn1 = new G4UnionSolid("BigMidUn1", bigBox, middleBox, trUn );
  G4UnionSolid * sidePlateFinal = new G4UnionSolid("SidePlateFinal", bigMidUn1, middleBox, trUn2);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *sidePlate = new G4LogicalVolume(sidePlateFinal, Copper,"SidePlate");

  return sidePlate;

}

G4LogicalVolume* GetConstructionTopPlate()
{
  G4double mainBoxX = 131.826*mm;
  G4double mainBoxY = 100.0252*mm;
  G4double mainBoxZ = 9.652*mm;
  G4double rmBoxX = 9.652*mm;
  G4double rmBoxY = 31.75*mm;
  G4double box1X = 45.72*mm;
  G4double box1Y = 11.43*mm;
  G4double boxZ = 22.352*mm;
  G4double box2X = 11.43*mm;
  G4double box2Y = 43.18*mm;

  G4Box * mainBox = new G4Box(" MainBox", mainBoxX/2, mainBoxY/2, mainBoxZ/2);
  G4Box * rmBox = new G4Box("RmBox", rmBoxX, rmBoxY/2, mainBoxZ);
  G4Box * box1 = new G4Box("Box1", box1X/2, box1Y/2, boxZ/2);
  G4Box * box2 = new G4Box("Box2", box2X/2, box2Y/2, boxZ/2);

  G4RotationMatrix* rot0 = new G4RotationMatrix;
  G4double posRMX = mainBoxX/2;
  G4double posBoxZ = boxZ/2 - mainBoxZ/2;
  G4double posBox1Y = 60.198- mainBoxY/2 + box1Y/2;
  G4double posBox2X = 77.724 - mainBoxX/2 + box2X/2;

  G4ThreeVector rmVec = G4ThreeVector(posRMX,0,0);
  G4ThreeVector box1Vec = G4ThreeVector(0,posBox1Y,posBoxZ);
  G4ThreeVector box2Vec = G4ThreeVector(posBox2X, 0, posBoxZ);

  G4Transform3D rmTr = G4Transform3D(*rot0, rmVec);
  G4Transform3D box1Tr = G4Transform3D(*rot0, box1Vec);
  G4Transform3D box2Tr = G4Transform3D(*rot0, box2Vec);

  G4SubtractionSolid * topPlate1 = new G4SubtractionSolid("TopPlate1", mainBox, rmBox, rmTr);
  G4UnionSolid * topPlate2 = new G4UnionSolid("TopPlate2", topPlate1, box1, box1Tr);
  G4UnionSolid* topPlateFinal = new G4UnionSolid("TopPlateFinal", topPlate2, box2, box2Tr);

  G4Material* Copper = G4Material::GetMaterial("G4_Cu");
  G4LogicalVolume *topPlate = new G4LogicalVolume(topPlateFinal, Copper,"TopPlate");

  return topPlate;

}
