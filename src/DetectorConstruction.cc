//This will be a simple goemetry just to check and see if I can get it working

#include "DetectorConstruction.hh"
//include header files for all classes used in the methods
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Isotope.hh"
#include "G4SystemOfUnits.hh"

//for sensitive detector definition
#include "SensitiveDetector.hh"
#include "SensitiveHDPE.hh"
#include "G4SDManager.hh"

#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "DetectorMessenger.hh"

//constructor / destructor do nothing
DetectorConstruction::DetectorConstruction(){
  AbsThickness = 2*223*mm;//subtract 10 mm for flattener. Doubled as no rear.
  DetMess = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction(){ 

}

void DetectorConstruction::SetAbsThickness(G4double value){
  AbsThickness = value;
}

G4double DetectorConstruction::GetAbsThickness(){
  return AbsThickness;
}

void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct(), true);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
/* materials definition */

/*define the elements that will be used in our materials*/
//define hydrogen 
  G4double A = 1.01 * g/mole;
  G4double Z = 1;
  G4Element* elH = new G4Element ("Hydrogen", "H", Z, A);

  //define oxygen
  A = 16.0 * g/mole;
   Z = 8;
  G4Element* elO = new G4Element ("Oxygen", "O", Z, A);

  //define nitrogen
  A = 14.0 * g/mole;
  Z = 7;
  G4Element* elN = new G4Element("Nitrogen", "N", Z, A);

  //define carbon
  A = 12.0107 * g/mole;
  Z = 6;
  G4Element* elC = new G4Element("Carbon", "C", Z, A);

  //define Silicon
  A = 28.086 * g/mole;
  Z = 14;
  G4Element* elSi = new G4Element("Silicon", "Si", Z, A);

  //define sodium
  A = 22.990 * g/mole;
  Z = 11;
  G4Element* elNa = new G4Element("Sodium", "Na", Z, A);

  //define Calcium
  A = 40.08 * g/mole;
  Z = 20;
  G4Element* elCa = new G4Element("Calcium", "Ca", Z, A);

  //define Phosphorus
  A = 30.973761*g/mole;
  Z = 15;
  G4Element* elP = new G4Element("Phosphorus", "P", Z, A);

  //define Sulphur
  A = 32.065*g/mole;
  Z = 16;
  G4Element* elS = new G4Element("Sulphur", "S", Z, A);

  //define Aluminium
  A = 26.961538*g/mole;
  Z = 13;
  G4Element* elAl = new G4Element("Aluminium", "Al",Z,A);

  //define Iron
  A = 55.845*g/mole;
  Z = 26;
  G4Element* elFe = new G4Element("Iron", "Fe",Z,A);
 
  //Define Copper
  A = 63.546*g/mole;
  Z = 29;
  G4Element* elCu = new G4Element("Copper", "Cu",Z,A);

  //Define Manganese
  A = 54.938045*g/mole;
  Z = 25;
  G4Element* elMn = new G4Element("Manganese", "Mn", Z, A);

  //Define Magnesium
  A = 24.3050*g/mole;
  Z = 12;
  G4Element* elMg = new G4Element("Magnesium", "Mg", Z, A);

  //Define Titanium
  A = 47.867*g/mole;
  Z = 22;
  G4Element* elTi = new G4Element("Titanium", "Ti", Z, A);

  //Define Chromium
  A = 51.9961*g/mole;
  Z = 24;
  G4Element* elCr = new G4Element("Chromium", "Cr",Z,A);

  //Define Zinc
  A = 65.409*g/mole;
  Z = 30;
  G4Element* elZn = new G4Element("Zinc", "Zn",Z,A);

  //Define Boron
  A = 10.811*g/mole;
  Z = 5;
  G4Element* elB = new G4Element("Boron", "B",Z,A);

  //Define Potassium
  A = 39.0983*g/mole;
  Z = 19;
  G4Element* elK = new G4Element("Potassium", "K",Z,A);

  //constructor of the G4Material class requires arguments: string conG4OpticalSurfacetaining name of material, density, number of elements
  G4Material* water = new G4Material("water", 1.0 * g/cm3, 2);
  water->AddElement(elH,2);
  water->AddElement(elO,1);

/*now we define air for the world volume*/
  G4Material* air = new G4Material("dry air", 0.01*mg/cm3, 2, kStateGas, 293*kelvin, 1*atmosphere);
//we can also specify the percentage (by mass) composition
  air->AddElement(elN, 75*perCent);
  air->AddElement(elO, 25*perCent);
  
  //Define Glass (Silicon Dioxide mixed with other stuff.Source - http://www.physics.nist.gov/cgi-bin/Star/compos.pl?matno=171)
  //G4Material* glass = new G4Material("Glass", 2.400 * g/cm3, 4, kStateSolid, 293*kelvin, 1*atmosphere);
  //now do fraction by weight of each element
  //glass->AddElement(elO,  0.459801);
  //glass->AddElement(elNa, 0.096441);
  //glass->AddElement(elSi, 0.336553);
  //glass->AddElement(elCa, 0.107205);

  //High Density Polyethylene
  //G4Material* HDPE = new G4Material("HD Poly-Ethylene", 0.941*g/cm3, 2);
  //G4Material* HDPE = new G4Material("HD Poly-Ethylene", 0.93*g/cm3, 2);
  G4Material* HDPE = new G4Material("HD Poly-Ethylene", 0.97*g/cm3, 2);
  HDPE->AddElement(elC, 33.33333*perCent);
  HDPE->AddElement(elH, 66.66667*perCent);

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* water_nist = nist->FindOrBuildMaterial("G4_WATER");
  G4Material* air_nist = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* HDPE_nist = nist->FindOrBuildMaterial("G4_POLYETHYLENE");


  //Ultima Gold Liquid Scintillator (Perkin-Elmer's advice)
  //G4Material* UltimaGold = new G4Material("Ultima Gold", 0.978*g/cm3, 7, kStateLiquid, 293*kelvin, 1*atmosphere);
  //UltimaGold->AddElement(elC, 0.39191371);
  //UltimaGold->AddElement(elH, 0.56822235);
  //UltimaGold->AddElement(elN, 0.00085188);
  //UltimaGold->AddElement(elO, 0.03542433);
  //UltimaGold->AddElement(elP, 0.00269660);
  //UltimaGold->AddElement(elS, 0.00037212);
  //UltimaGold->AddElement(elNa,0.00051901);

  //Insta Gel Plus Liquid Scintillator (from Perkin-Elmer's site)
  //G4Material* InstaGel = new G4Material("Insta-Gel", 0.960*g/cm3, 4, kStateLiquid, 293*kelvin, 1*atmosphere);
  //InstaGel->AddElement(elC, 0.348865649);
  //InstaGel->AddElement(elH, 0.577868215);
  //InstaGel->AddElement(elN, 0.000084866);
  //InstaGel->AddElement(elO, 0.073181270);

  //using another variable 'Cocktail' to allow inclusion of other scintillants.
  //G4Material* Cocktail = new G4Material("Scint Cocktail", 1*g/cm3, 1, kStateLiquid, 293*kelvin, 1*atmosphere);
  //Cocktail->AddMaterial(InstaGel, 1);

  //Aluminium Alloy used in the walls Alloy #6063
  G4Material* AlAlloy = new G4Material("Al Alloy", 2680.*kg/m3, 9, kStateSolid);
  //Adding elements in this way seems pointless now but is useful for later
  AlAlloy->AddElement(elAl, 98.5*perCent);
  AlAlloy->AddElement(elSi, 0.40*perCent);
  AlAlloy->AddElement(elFe, 0.175*perCent);
  AlAlloy->AddElement(elCu, 0.05*perCent);
  AlAlloy->AddElement(elMn, 0.05*perCent);
  AlAlloy->AddElement(elMg, 0.675*perCent);
  AlAlloy->AddElement(elCr, 0.05*perCent);
  AlAlloy->AddElement(elZn, 0.05*perCent);
  AlAlloy->AddElement(elTi, 0.05*perCent);

  //assumed that cork has the same composition as cellulose.
  //Nitrogen makes up the rounding error.
  G4Material* Cork = new G4Material("Cork", 200*kg/m3,4);
  Cork->AddElement(elC, 0.28571);
  Cork->AddElement(elH, 0.47619);
  Cork->AddElement(elO, 0.23809);
  Cork->AddElement(elN, 0.00001);


/*the volumes: */

  G4double worldx = 1 * m;
  G4double worldy = 1 * m;
  G4double worldz = 1 * m;

//the whole simulation must be contained within a "world" volume
//defining a volume requires definition of solid, logical, and physical volume
//solid volume is the shape, has dimensions
  G4Box* world = new G4Box("world_box", worldx, worldy, worldz);

//logical volume: has a solid volume as a member, a material, last 3???
  G4LogicalVolume* logical_world = new G4LogicalVolume(world, air_nist, "world_log", 0,0,0);

  //make the world invisible!
  logical_world->SetVisAttributes(G4VisAttributes::Invisible);

  //physical volume: G4PVPlacement class, represents a logical volume that is placed somewhere in the mother volumeG4OpticalSurface
  G4VPhysicalVolume* physical_world = new G4PVPlacement(0,G4ThreeVector(),logical_world, "world_phys", 0, false, 0);



  //now to make the curve of the vial in the middle (volume representing and
  //the contents inside will be included as daughter).
  //Arguments: Name (string), Inner Radius (double), Outer radius (Double), Half length in z axis of tube, starting angle (double), spanning angle of tube - ie for a full circle make it 360 degrees

    G4double R_inner = 12.72*mm;
    G4double R_outer = 13.77 * mm;
    G4double R_outerV  = R_outer;
    G4double HalfLength = 24.00 * mm;
    G4double StartPhi = 0*deg;
    G4double SpanningPhi = 360*deg;

    G4Tubs* vialCurve = new G4Tubs("vialCurve", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_vialCurve = new G4LogicalVolume(vialCurve, HDPE_nist, "vialCurve_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_vialCurve = new G4OpticalSurface("vialCurve optical surface", unified, ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_vialCurve = new G4LogicalSkinSurface("Vial Curved Surface",logical_vialCurve,opt_surface_vialCurve);

    G4VisAttributes* VA_vialCurve = new G4VisAttributes();
    VA_vialCurve->SetColour(1,0,1);
    logical_vialCurve->SetVisAttributes(VA_vialCurve);

    G4VPhysicalVolume* physical_vialCurve = new G4PVPlacement(0,G4ThreeVector(0,0,0),logical_vialCurve,"vialCurve_phys", logical_world, false, 0);


    //now for the base of the vial. A flat Cylinder that sits snug inside vialCurve.
    R_inner = 0*mm;
    R_outer = 12.72 * mm;
    HalfLength = 0.425 * mm;

    G4Tubs* vialBase = new G4Tubs("vialBase", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_vialBase = new G4LogicalVolume(vialBase, HDPE_nist, "vialBase_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_vialBase = new G4OpticalSurface("vialBase optical surface", unified, ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_vialBase = new G4LogicalSkinSurface("Vial Base Surface",logical_vialBase,opt_surface_vialBase);

    G4VisAttributes* VA_vialBase = new G4VisAttributes();
    VA_vialBase->SetColour(1,0,1);
    logical_vialBase->SetVisAttributes(VA_vialBase);

    //Placement is -24mm+0.425 in z axis (level with vialCurve as it sits snug.
    G4VPhysicalVolume* physical_vialBase = new G4PVPlacement(0,G4ThreeVector(0,0,-23.575*mm),logical_vialBase,"vialBase_phys", logical_world, false, 0);


    //ok, now for the vial top. An annulus at the top of the vial. The neck will sit level with (ie not on top of) this, and
    //this sits level and fits snugly with the top of vialCurve
    R_inner = 9.68*mm;
    R_outer = 12.72 * mm;
    HalfLength = 0.525 * mm;

    G4Tubs* vialTop = new G4Tubs("vialTop", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_vialTop = new G4LogicalVolume(vialTop, HDPE_nist, "vialTop_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_vialTop = new G4OpticalSurface("vialTop opt surf", unified, ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_vialTop = new G4LogicalSkinSurface("Vial Top Surface",logical_vialTop,opt_surface_vialTop);

    G4VisAttributes* VA_vialTop = new G4VisAttributes();
    VA_vialTop->SetColour(1,0,1);
    logical_vialTop->SetVisAttributes(VA_vialTop);

    //Placement is +24mm - 0.525mm in z axis (so the top is level with the top of vialCurve).
    G4VPhysicalVolume* physical_vialTop = new G4PVPlacement(0,G4ThreeVector(0,0,23.475*mm),logical_vialTop,"vialTop_phys", logical_world, false, 0);


    //The vial Neck. An annulus at the top of the vial. This will sit level with (ie not on top of) the vial top.
    R_inner = 8.00*mm;
    R_outer = 9.68 * mm;
    //halflegth is (10+1.05)/2
    HalfLength = 5.525 * mm;

    G4Tubs* vialNeck = new G4Tubs("vialNeck", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_vialNeck = new G4LogicalVolume(vialNeck, HDPE_nist, "vialNeck_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_vialNeck = new G4OpticalSurface("vialNeck opt surf", unified,ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_vialNeck = new G4LogicalSkinSurface("Vial Neck Surface",logical_vialNeck,opt_surface_vialNeck);

    G4VisAttributes* VA_vialNeck = new G4VisAttributes();
    VA_vialNeck->SetColour(1,0,1);
    //VA_vialNeck->SetForceSolid(true);
    logical_vialNeck->SetVisAttributes(VA_vialNeck);

    //Placement is +24mm - 1.05 + 5.525mm in z axis (so the bottom is level with the inside of the vial - or bottom of vialTop).
    G4VPhysicalVolume* physical_vialNeck = new G4PVPlacement(0,G4ThreeVector(0,0,28.475*mm),logical_vialNeck,"vialNeck_phys", logical_world, false, 0);


    //The Alminium reflector. A flat cylinder sitting on top of the vial neck.
    R_inner = 0.00*mm;
    R_outer = 11.025 * mm;
    HalfLength = 0.065 * mm;

    G4Tubs* AlReflect = new G4Tubs("Al Reflector", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_AlRefl = new G4LogicalVolume(AlReflect, HDPE_nist, "AlRefl_log", 0,0,0);
   
    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_AlRefl = new G4OpticalSurface("AlRefl optical surface", unified, groundfrontpainted, dielectric_metal,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_AlRefl = new G4LogicalSkinSurface("Al Reflector Surface",logical_AlRefl,opt_surface_AlRefl);

    G4VisAttributes* VA_AlRefl = new G4VisAttributes();
    VA_AlRefl->SetColour(1,0.5,1);
    logical_AlRefl->SetVisAttributes(VA_AlRefl);

    //Placement is +24mm + 10mm + 0.065mm in z axis (so the top is level with the top of the vial neck).
    G4VPhysicalVolume* physical_AlRefl = new G4PVPlacement(0,G4ThreeVector(0,0,34.065*mm),logical_AlRefl,"AlRefl_phys", logical_world, false, 0);

    //The Cork backing. A flat cylinder sitting on top of the Al Reflector.
    R_inner = 0.00*mm;
    R_outer = 11.025 * mm;
    HalfLength = 0.77 * mm;

    G4Tubs* CorkBack = new G4Tubs("Cork Backing", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_CorkBack = new G4LogicalVolume(CorkBack, HDPE_nist, "CorkBack_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_CorkBack = new G4OpticalSurface("CorkBack opt surf", unified,ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_CorkBack = new G4LogicalSkinSurface("CorkBack Surface",logical_CorkBack,opt_surface_CorkBack);

    G4VisAttributes* VA_CorkBack = new G4VisAttributes();
    VA_CorkBack->SetColour(1,0.75,1);
    VA_CorkBack->SetForceSolid(true);
    logical_CorkBack->SetVisAttributes(VA_CorkBack);

    //Placement is +24mm + 10mm + 0.13mm + 0.77mm in z axis (so the bottom is level with the top of the Al Reflector).
    G4VPhysicalVolume* physical_CorkBack = new G4PVPlacement(0,G4ThreeVector(0,0,34.9*mm),logical_CorkBack,"CorkBack_phys", logical_world, false, 0);

    //The lid top. A flat cylinder sitting on top of the cork backing.
    R_inner = 0.00*mm;
    R_outer = 11.025 * mm;
    HalfLength = 1.3 * mm;

    G4Tubs*  LidTop= new G4Tubs("Lid Top", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_LidTop = new G4LogicalVolume(LidTop, HDPE_nist, "LidTop_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_LidTop = new G4OpticalSurface("LidTop opt surf", unified,ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_LidTop = new G4LogicalSkinSurface("LidTop Surface",logical_LidTop,opt_surface_LidTop);

    G4VisAttributes* VA_LidTop = new G4VisAttributes();
    VA_LidTop->SetColour(1,1,0);
    VA_LidTop->SetForceSolid(true);
    logical_LidTop->SetVisAttributes(VA_LidTop);

    //Placement is +24mm + 10mm + 0.13mm + 1.54mm + 1.3mm in z axis (sitting on cork backing).
    G4VPhysicalVolume* physical_LidTop = new G4PVPlacement(0,G4ThreeVector(0,0,36.97*mm),logical_LidTop,"LidTop_phys", logical_world, false, 0);


    //The lid curve. A long annular cylinder sitting so its top is level with the top of LidTop.
    R_inner = 11.025*mm;
    R_outer = 12.25 * mm;
    HalfLength = 6.6025 * mm;

    G4Tubs*  LidCurve = new G4Tubs("Lid Curve", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_LidCurve = new G4LogicalVolume(LidCurve, HDPE_nist, "LidCurve_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_LidCurve = new G4OpticalSurface("LidCurve opt surf", unified,ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_LidCurve = new G4LogicalSkinSurface("LidCurve Surface",logical_LidCurve,opt_surface_LidCurve);

    G4VisAttributes* VA_LidCurve = new G4VisAttributes();
    VA_LidCurve->SetColour(1,1,0);
    VA_LidCurve->SetForceSolid(true);
    logical_LidCurve->SetVisAttributes(VA_LidCurve);

    //Placement is +24mm + 2.145mm(gap between bottom of the lid and top of vial) + 6.6025mm in z axis.
    //alternatively try 36.97mm (middle of lidTop)+ 1.3mm (lidTop half lenght) - 6.0625mm
    G4VPhysicalVolume* physical_LidCurve = new G4PVPlacement(0,G4ThreeVector(0,0,31.6675*mm),logical_LidCurve,"LidCurve_phys", logical_world, false, 0);


    //define the liquid bit inside the vial (walls are 1.05mm thick) the liquid depth is correct for 10mL
    //of scintillant. Liquid must sit flush with the top of VialBottom.
    //Modifying liquid volume to 15 mL.
    R_inner = 0;
    R_outer = 12.72 * mm;
    HalfLength = (9.83662*1.5)*mm;

    G4Tubs* liquid = new G4Tubs("Liquid", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_liquid = new G4LogicalVolume(liquid,water_nist, "liquid_log", 0,0,0);
        
    //create SensitiveDetector object
    SensitiveDetector* SDscintillator = new SensitiveDetector("liquid_log");
    //get pointer to the sensitive detector manager: this class is used by G4RunManager to see which sensitive detectors there are
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //pass to manager
    SDman->AddNewDetector(SDscintillator);

    //now we pass the sensitive detector pointer to the logical volume of our scoring volume
    logical_liquid->SetSensitiveDetector(SDscintillator);

    G4VisAttributes* VA_liq = new G4VisAttributes();
    VA_liq->SetColour(0,1,1);
    VA_liq->SetForceSolid(true);
    logical_liquid->SetVisAttributes(VA_liq);

    //Liquid centre is -24mm + 0.85mm + 9.83662mm
    G4double liquidCentre = -24*mm+0.85*mm+HalfLength;
    G4VPhysicalVolume* physical_liquid = new G4PVPlacement(0,G4ThreeVector(0,0,liquidCentre),logical_liquid,"liquid_phys", logical_world, false, 0);

    //now for the finishing touches, just to show off. Add the screw bits, one ring of HDPE_nist protruding from the neck
    //and one ring of plastic (lower down) protruding from the inside of the lid.
    //neck ring will sit 2.5mm lower than the Al Reflector
    R_inner = 9.68*mm;
    R_outer = 10.205* mm;
    HalfLength = 0.525 * mm;

    G4Tubs* NeckRing = new G4Tubs("Neck Ring", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_NeckRing = new G4LogicalVolume(NeckRing, HDPE_nist, "NeckRing_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_NeckRing = new G4OpticalSurface("NeckRing opt surf", unified,ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_NeckRing = new G4LogicalSkinSurface("NeckRing Surface",logical_NeckRing,opt_surface_NeckRing);

    G4VisAttributes* VA_NeckRing = new G4VisAttributes();
    VA_NeckRing->SetColour(1,0,1);
    //VA_NeckRing->SetForceSolid(true);
    logical_NeckRing->SetVisAttributes(VA_NeckRing);

    //NeckRing centre is 34.0mm (top of neck) - 2.5mm.
    G4VPhysicalVolume* physical_NeckRing = new G4PVPlacement(0,G4ThreeVector(0,0,31.5*mm),logical_NeckRing,"NeckRing_phys", logical_world, false, 0);

    //lid ring will sit 5mm lower than the Al Reflector
    R_inner = 9.975*mm;
    R_outer = 11.025* mm;
    HalfLength = 0.525 * mm;

    G4Tubs* LidRing = new G4Tubs("Lid Ring", R_inner, R_outer, HalfLength, StartPhi, SpanningPhi);

    G4LogicalVolume* logical_LidRing = new G4LogicalVolume(LidRing, HDPE_nist, "LidRing_log", 0,0,0);

    //arguments of OpticalSurface:name,model,finish,type of interface,?(double)
    G4OpticalSurface* opt_surface_LidRing = new G4OpticalSurface("LidRing opt surf", unified,ground, dielectric_dielectric,1.);

    //arguments of LogicalSkinSurface: name, logical vol pointer, op properties pointer
    G4LogicalSkinSurface* skin_LidRing = new G4LogicalSkinSurface("LidRing Surface",logical_LidRing,opt_surface_LidRing);

    G4VisAttributes* VA_LidRing = new G4VisAttributes();
    VA_LidRing->SetColour(1,1,0);
    //VA_LidRing->SetForceSolid(true);
    logical_LidRing->SetVisAttributes(VA_LidRing);

    //LidRing centre is 34.0mm (top of neck) - 5mm.
    G4VPhysicalVolume* physical_LidRing = new G4PVPlacement(0,G4ThreeVector(0,0,29*mm),logical_LidRing,"LidRing_phys", logical_world, false, 0);



    //OK, now I need to add the HDPE_nist absorber.
    //G4Box arguments: name, 0.5*(x dimension, y dimension, z dimension).
    G4double HDPE_nistthickness = AbsThickness;//Defaults to 233*mm;
    G4Box* HDPE_nistabs = new G4Box("HDPE_nist absorber", 100*mm,
				    (HDPE_nistthickness/2), 100*mm);

    G4LogicalVolume* logical_HDPE_nistabs =
      new G4LogicalVolume(HDPE_nistabs, HDPE_nist, "HDPE_nistabs_log", 0, 0, 0);


    //create SensitiveDetector object
    SensitiveHDPE* SDHDPE = new SensitiveHDPE("HDPE_nistabs_log");
    //pass new sensitive detector to manager
    SDman->AddNewDetector(SDHDPE);

    //now we pass the sensitive detector pointer to the logical volume of our
    //scoring volume
    logical_HDPE_nistabs->SetSensitiveDetector(SDHDPE);


    G4VisAttributes* VA_HDPE_nistabs = new G4VisAttributes();
    VA_HDPE_nistabs->SetColour(0.8,0.8,0.8);
    //VA_HDPE_nistabs->SetForceSolid(true);
    logical_HDPE_nistabs->SetVisAttributes(VA_HDPE_nistabs);

    //HDPE_nistabs centre is
    //vial outer diameter + 5mm (flattener) + 1/2*HDPE_nistthickness.
    G4double HDPE_nistcentreY = R_outerV + 5*mm + HDPE_nistthickness/2;
    G4VPhysicalVolume* physical_HDPE_nistabs =
      new G4PVPlacement(0, G4ThreeVector(0,HDPE_nistcentreY+1,0), logical_HDPE_nistabs,
			"HDPE_nistabs_phys", logical_world, false, 0);

    //G4double Vial2Centre =
    //HDPE_nistcentreY + R_outerV + HDPE_nistthickness/2 + 5;//made 5mm off.

    //Place vial2 components
    //G4VPhysicalVolume* physical_vial2Curve =
    //new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,0),logical_vialCurve,
    //			"vial2Curve_phys", logical_world, false, 0);

//G4VPhysicalVolume* physical_vial2Base =
//    new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,-23.575*mm),
    //logical_vialBase,"vial2Base_phys", logical_world,
			  //			false, 0);

//  G4VPhysicalVolume* physical_vial2Top = 
//    new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,23.475*mm),
//			logical_vialTop,"vial2Top_phys", logical_world,
//			false, 0);

//  G4VPhysicalVolume* physical_vial2Neck =
//    new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,28.475*mm),
    //logical_vialNeck,"vial2Neck_phys", logical_world,
			  //			false, 0);
    
    //Not including the lid liner in the second vial.
//  G4VPhysicalVolume* physical_Lid2Top =
//    new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,36.97*mm),logical_LidTop,
//			"Lid2Top_phys", logical_world, false, 0);

//  G4VPhysicalVolume* physical_Lid2Curve =
    //    new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,31.6675*mm),
    //			logical_LidCurve,"Lid2Curve_phys", logical_world,
    //			false, 0);


    //Now add the liquid scintillator for the second vial.
    //G4LogicalVolume* logical_liquid2 =
    //new G4LogicalVolume(liquid,water_nist, "liquid2_log", 0,0,0);
        
    //create SensitiveDetector object
    //SensitiveDetectorRear* SDscintillator2 = 
    //new SensitiveDetectorRear("liquid2_log");
    //pass new sensitive detector to manager
    //SDman->AddNewDetector(SDscintillator2);

    //now we pass the sensitive detector pointer to the logical volume of our
    //scoring volume
    //logical_liquid2->SetSensitiveDetector(SDscintillator2);

    //logical_liquid2->SetVisAttributes(VA_liq);

    //Place the liquid
    //G4VPhysicalVolume* physical_liquid2 =
    //new G4PVPlacement(0,G4ThreeVector(0,Vial2Centre,liquidCentre),
    //			logical_liquid2,"liquid2_phys",logical_world, false, 0);


    //Now add another vial at the front, this one below the existing front vial.
    //There is a 10 mm gap between the top of this and the bottom of the next
    //vial.
    //Therefore the centre of this vial is at -72.27 mm along z.
    G4double Vial3Centre = -72.27*mm;

    //Place vial3 components
    G4VPhysicalVolume* physical_vial3Curve =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial3Centre),logical_vialCurve,
			"vial3Curve_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_vial3Base =
      new G4PVPlacement(0,G4ThreeVector(0,0,-23.575*mm+Vial3Centre),
			logical_vialBase,"vial3Base_phys", logical_world,
			false, 0);

    G4VPhysicalVolume* physical_vial3Top = 
      new G4PVPlacement(0,G4ThreeVector(0,0,23.475*mm+Vial3Centre),
    			logical_vialTop,"vial3Top_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_vial3Neck =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial3Centre+28.475*mm),
    			logical_vialNeck,"vial3Neck_phys", logical_world,
    			false, 0);
    
    //Not including the lid liner in the second vial.
    G4VPhysicalVolume* physical_Lid3Top =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial3Centre+36.97*mm),
    			logical_LidTop,"Lid3Top_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_Lid3Curve =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial3Centre+31.6675*mm),
    			logical_LidCurve,"Lid3Curve_phys", logical_world,
    			false, 0);

    //Now add the liquid scintillator for the non-SD vials.
    G4LogicalVolume* logical_liquid_nonSD =
      new G4LogicalVolume(liquid,water_nist, "liquidNonSD_log", 0,0,0);

    G4VPhysicalVolume* physical_liquid_vial3 =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial3Centre+liquidCentre),
    			logical_liquid_nonSD,"liquid3_phys",
    			logical_world, false, 0);

    //The vial next to the centre vial has a 5 mm gap, so the centre of the vial
    //is at 5 + R_outerV*2 mm along X.
    G4double Vial4Centre = R_outerV*2 + 5*mm;
    
    //Place vial4 components
    G4VPhysicalVolume* physical_vial4Curve =
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,0),logical_vialCurve,
    			"vial4Curve_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_vial4Base =
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,-23.575*mm),
			logical_vialBase,"vial4Base_phys", logical_world,
			false, 0);

    G4VPhysicalVolume* physical_vial4Top = 
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,23.475*mm),
    			logical_vialTop,"vial4Top_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_vial4Neck =
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,28.475*mm),
    			logical_vialNeck,"vial4Neck_phys", logical_world,
    			false, 0);
    
    //Not including the lid liner in the second vial.
    G4VPhysicalVolume* physical_Lid4Top =
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,36.97*mm),
    			logical_LidTop,"Lid4Top_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_Lid4Curve =
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,31.6675*mm),
    			logical_LidCurve,"Lid4Curve_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_liquid_vial4 =
      new G4PVPlacement(0,G4ThreeVector(Vial4Centre,0,liquidCentre),
    			logical_liquid_nonSD,"liquid4_phys",
    			logical_world, false, 0);

    //The vial next to the centre vial has a 5 mm gap, so the centre of the vial
    //is at 5 + R_outerV*2 mm along X.
    G4double Vial5Centre = -(R_outerV*2 + 5*mm);
    
    //Place vial5 components
    G4VPhysicalVolume* physical_vial5Curve =
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,0),logical_vialCurve,
    			"vial5Curve_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_vial5Base =
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,-23.575*mm),
    			logical_vialBase,"vial5Base_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_vial5Top = 
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,23.475*mm),
    			logical_vialTop,"vial5Top_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_vial5Neck =
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,28.475*mm),
    			logical_vialNeck,"vial5Neck_phys", logical_world,
    			false, 0);
    
    //Not including the lid liner in the second vial.
    G4VPhysicalVolume* physical_Lid5Top =
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,36.97*mm),
    			logical_LidTop,"Lid5Top_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_Lid5Curve =
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,31.6675*mm),
    			logical_LidCurve,"Lid5Curve_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_liquid_vial5 =
      new G4PVPlacement(0,G4ThreeVector(Vial5Centre,0,liquidCentre),
    			logical_liquid_nonSD,"liquid5_phys",
    			logical_world, false, 0);

    //Now add a vial above the centre of the front vial.
    G4double Vial6Centre = 72.27*mm;

    //Place vial6 components
    G4VPhysicalVolume* physical_vial6Curve =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial6Centre),logical_vialCurve,
    			"vial6Curve_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_vial6Base =
      new G4PVPlacement(0,G4ThreeVector(0,0,-23.575*mm+Vial6Centre),
    			logical_vialBase,"vial6Base_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_vial6Top = 
      new G4PVPlacement(0,G4ThreeVector(0,0,23.475*mm+Vial6Centre),
    			logical_vialTop,"vial6Top_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_vial6Neck =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial6Centre+28.475*mm),
    			logical_vialNeck,"vial6Neck_phys", logical_world,
			false, 0);
    
    //Not including the lid liner in the second vial.
    G4VPhysicalVolume* physical_Lid6Top =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial6Centre+36.97*mm),
    			logical_LidTop,"Lid6Top_phys", logical_world, false, 0);

    G4VPhysicalVolume* physical_Lid6Curve =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial6Centre+31.6675*mm),
    			logical_LidCurve,"Lid6Curve_phys", logical_world,
    			false, 0);

    G4VPhysicalVolume* physical_liquid_vial6 =
      new G4PVPlacement(0,G4ThreeVector(0,0,Vial6Centre+liquidCentre),
    			logical_liquid_nonSD,"liquid6_phys",
    			logical_world, false, 0);

    //Now make a HDPE flattener, 20x20 cm aligned with beam centre, and with
    //G4double flatHalfThickness = R_outerV + 5*mm;
    //Now that I've made the mount, adding more realistic dimensions.
    G4double flatHalfThickness = 19*mm;

    G4Box* flattenerBox = new G4Box("FlattenerBox",
				    10*cm, flatHalfThickness, 10*cm);

    //Now make the piece to make holes for the subtraction solid.
    //Radius of holes is vial outer radius + 1 mm. length just needs to be >20cm
    //G4Tubs* thruHole =
    //new G4Tubs("Through hole", 0, (R_outerV+1*mm), 21*cm, 0, 360*deg);
    G4Tubs* thruHole =
      new G4Tubs("Through hole", 0, 30*mm, 21*cm, 0, 360*deg);

    //Subtract central hole
    G4VSolid *Flat1Hole =
      new G4SubtractionSolid("Flat1Hole", flattenerBox, thruHole,
			     0, G4ThreeVector(0,0,0));

    //Subtract holes on either side, offset by Vial4Centre and Vial5Centre along
    //X.
    G4VSolid* Flat2Holes =
      new G4SubtractionSolid("Flat2Holes", Flat1Hole, thruHole,
			     0, G4ThreeVector(Vial4Centre, 0, 0));

    G4VSolid* Flat3Holes = 
      new G4SubtractionSolid("Flat3Holes", Flat2Holes, thruHole,
			     0, G4ThreeVector(Vial5Centre, 0, 0));

    G4LogicalVolume* Flattener_log = 
      new G4LogicalVolume(Flat3Holes, HDPE_nist, "Flattener_log", 0,0,0);

    G4VPhysicalVolume* Flattener_phys = 
      new G4PVPlacement(0, G4ThreeVector(0,0,0), Flattener_log, "Flattener_phys"
			, logical_world, false, 0);

    //Place another flattener centred on Vial2Centre.
//    G4VPhysicalVolume* Flattener_phys2 = 
//    new G4PVPlacement(0, G4ThreeVector(0,Vial2Centre,0), Flattener_log,
//			"Flattener_phys", logical_world, false, 0);

    return physical_world;
}
