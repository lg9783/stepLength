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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "globals.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include <G4SystemOfUnits.hh>
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "Randomize.hh"

#include <fstream>
#include "CommandLineParser.hh"


using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::angstrom;
using CLHEP::degree;
using CLHEP::micrometer;
using CLHEP::mm;
using CLHEP::nanometer;

static G4VisAttributes visInvBlue(false, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visInvWhite(false, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visInvPink(false, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visInvCyan(false, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visInvRed(false, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visInvGreen(false, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visWhite(true, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visPink(true, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visCyan(true, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visRed(true, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visGreen(true, G4Colour(0.0, 1.0, 0.0));

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),
                                               fpGun(new G4MoleculeGun())
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  /***************************************************************************/
  //                               World
  /***************************************************************************/

  G4NistManager *man = G4NistManager::Instance();
  G4Material *waterMaterial = man->FindOrBuildMaterial("G4_WATER");
  G4Material *air = man->FindOrBuildMaterial("G4_AIR");

  G4Box *solidWorld{nullptr};
  G4Orb *solidWaterBox{nullptr};

  if (CommandLineParser::GetParser()->GetCommandIfActive("-ref") == 0)
  {

    solidWorld = new G4Box("world", 1 * mm, 1 * mm, 1 * mm);

    solidWaterBox = new G4Orb("waterBox", 10 * um);
  }
  else if (CommandLineParser::GetParser()->GetCommandIfActive("-ref") != 0)
  {
    // larger world volume to provide build up for CPE for photon beam
    solidWorld = new G4Box("world", 100 * mm, 100 * mm, 100 * mm);

    solidWaterBox = new G4Orb("waterBox", 2.5 * mm);
  }

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    air,
                                                    "world");

  G4PVPlacement *physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicWaterBox = new G4LogicalVolume(solidWaterBox,
                                                       waterMaterial,
                                                       "waterBox");

  G4PVPlacement *physiWaterBox = new G4PVPlacement(0,
                                                   G4ThreeVector(),
                                                   logicWaterBox,
                                                   "waterBox",
                                                   logicWorld,
                                                   0,
                                                   false,
                                                   0);
  logicWorld->SetVisAttributes(&visInvWhite);
  logicWaterBox->SetVisAttributes(&visInvWhite);

  G4Tubs *solidTrackingVol = new G4Tubs("TrackingVol", 0, 39 * nm, 104 * nm, 0., 2. * CLHEP::pi); // volume to track radicals in 9nm larger than chromatin
  G4LogicalVolume *logicTrackingVol = new G4LogicalVolume(solidTrackingVol,
                                                          waterMaterial,
                                                          "TrackingVol");
  G4PVPlacement *physiTrackingVol = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      logicTrackingVol,
                                                      "TrackingVol",
                                                      logicWaterBox,
                                                      0,
                                                      false,
                                                      0);
  logicTrackingVol->SetVisAttributes(&visInvWhite);

  G4Tubs *solidChromatinSegment = new G4Tubs("chromatinSegment", 0, 30 * nm, 95 * nm, 0., 2. * CLHEP::pi);
  G4LogicalVolume *logicChromatinSegment = new G4LogicalVolume(solidChromatinSegment,
                                                               waterMaterial,
                                                               "chromatinSegment");
  G4PVPlacement *physiChromatinSegment = new G4PVPlacement(0,
                                                           G4ThreeVector(),
                                                           logicChromatinSegment,
                                                           "chromatinSegment",
                                                           logicTrackingVol,
                                                           0,
                                                           false,
                                                           0);
  logicChromatinSegment->SetVisAttributes(&visPink);

  // G4cout << G4BestUnit(solidChromatinSegment->GetCubicVolume(), "Volume") << G4endl;

  G4double sugarRadius{2.9 * angstrom};

  G4Orb *solidSugar = new G4Orb("sugar", sugarRadius);
  G4LogicalVolume *logicSugar0 = new G4LogicalVolume(solidSugar,
                                                     waterMaterial,
                                                     "sugar0");

  G4LogicalVolume *logicSugar1 = new G4LogicalVolume(solidSugar,
                                                     waterMaterial,
                                                     "sugar1");

  G4Tubs *solidHistone = new G4Tubs("histone",
                                    0,
                                    5.5 * nanometer,
                                    2.75 * nanometer,
                                    0 * degree,
                                    360 * degree);
  G4LogicalVolume *logicHistone = new G4LogicalVolume(solidHistone,
                                                      waterMaterial,
                                                      "histone");

  ifstream f("sugarPos.csv", ios::in);

  std::vector<G4ThreeVector> fPositions0;
  std::vector<G4ThreeVector> fPositions1;
  std::vector<G4ThreeVector> fPositionsBase0;
  std::vector<G4ThreeVector> fPositionsBase1;

  while (!f.eof())

  {
    double x1, y1, z1, x2, y2, z2, base_x1, base_y1, base_z1, base_x2, base_y2, base_z2;
    f >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> base_x1 >> base_y1 >> base_z1 >> base_x2 >> base_y2 >> base_z2;
    fPositions0.push_back(G4ThreeVector(x1,
                                        y1,
                                        z1));
    fPositions1.push_back(G4ThreeVector(x2,
                                        y2,
                                        z2));
    fPositionsBase0.push_back(G4ThreeVector(base_x1,
                                            base_y1,
                                            base_z1));
    fPositionsBase1.push_back(G4ThreeVector(base_x2,
                                            base_y2,
                                            base_z2));
    // G4cout << x1 << G4endl;
  }

  for (size_t i = 0; i < fPositions1.size(); ++i) // 0 and 1 are the same size
  {
    new G4PVPlacement(0,
                      fPositions0[i],
                      logicSugar0,
                      "sugar0",
                      logicChromatinSegment,
                      false,
                      i,
                      false);
    fpGun->AddMolecule("Deoxyribose",
                       fPositions0[i],
                       1 * picosecond);

    G4double R = G4UniformRand();
    if (R < 0.5) // Add bases randomly
    {
      fpGun->AddMolecule("Cytosine",
                         fPositionsBase0[i],
                         1 * picosecond);
    }
    else
    {
      fpGun->AddMolecule("Adenine",
                         fPositionsBase0[i],
                         1 * picosecond);
    }

    new G4PVPlacement(0,
                      fPositions1[i],
                      logicSugar1,
                      "sugar1",
                      logicChromatinSegment,
                      false,
                      i,
                      false);
    fpGun->AddMolecule("Deoxyribose",
                       fPositions1[i],
                       1 * picosecond);

    if (R < 0.5)
    {
      fpGun->AddMolecule("Guanine", // Add bases randomly
                         fPositionsBase1[i],
                         1 * picosecond);
    }
    else
    {
      fpGun->AddMolecule("Thymine ",
                         fPositionsBase1[i],
                         1 * picosecond);
    }
  }

  // Histones
  ifstream f2("histonePositions.csv", ios::in);

  std::vector<G4ThreeVector> fPositions;
  std::vector<G4RotationMatrix *> fRotations;

  while (!f2.eof())

  {
    double x1, y1, z1, rotx, roty, rotz;
    f2 >> x1 >> y1 >> z1 >> rotx >> roty >> rotz;
    fPositions.push_back(G4ThreeVector(x1,
                                       y1,
                                       z1));
    fRotations.push_back(new G4RotationMatrix(roty, // geant convention is y,z,x
                                              rotz,
                                              rotx));
  }
  for (size_t i = 0; i < fPositions.size(); ++i)
  {
    new G4PVPlacement(fRotations[i],
                      fPositions[i],
                      logicHistone,
                      "histone",
                      logicChromatinSegment,
                      false,
                      false);
  }

  logicSugar0->SetVisAttributes(&visInvWhite);
  logicSugar1->SetVisAttributes(&visInvRed);
  logicHistone->SetVisAttributes(&visBlue);

  G4UserLimits *l = new G4UserLimits();

  G4double percent = 0.1;
  // Sets a max Step length in the chromatin fibre, sugar and histone volumes:
  G4double maxStep = sugarRadius * 2 * percent; // set max step length to a percentage of sugar diameter
  l->SetMaxAllowedStep(maxStep);

  logicSugar0->SetUserLimits(l);
  logicSugar1->SetUserLimits(l);
  logicChromatinSegment->SetUserLimits(l);
  logicHistone->SetUserLimits(l);

  return physiWorld;
}
