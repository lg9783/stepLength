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
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"
#include <G4Scheduler.hh>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "G4AnalysisManager.hh"
#include "G4DNAMolecule.hh"
#include "G4MoleculeTable.hh"
#include "G4OH.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TransportationManager.hh"
#include "CommandLineParser.hh"


using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction() 
    : G4UserTimeStepAction()
{
    // AddTimeStep(1*picosecond, 0.35*picosecond);
    // AddTimeStep(10*picosecond, 1*picosecond);
    // AddTimeStep(100*picosecond, 3*picosecond);
    // AddTimeStep(1000*picosecond, 10*picosecond);
    // AddTimeStep(10000*picosecond, 100*picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(const TimeStepAction& other) 
    : G4UserTimeStepAction(other)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction&
TimeStepAction::operator=(const TimeStepAction& rhs)
{
    if (this == &rhs) return *this;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void TimeStepAction::UserPreTimeStepAction()
{
// G4cout << "time step" << G4Scheduler::Instance()->GetGlobalTime()<< G4endl;
G4TrackManyList* allTrackList = G4ITTrackHolder::Instance()->GetMainList();
G4TrackManyList::iterator it = allTrackList->begin();
G4TrackManyList::iterator end = allTrackList->end();


for(;it!=end;++it)
{
G4Track* track = *it; // track can be an OH, e_aq, H2, ...
G4String name = GetMolecule(track)->GetName();

auto localPos = track->GetPosition();

G4VTouchable* newTouchable = CreateTouchableAtPoint(localPos);

// G4cout << name << G4endl;


G4String volumeName = newTouchable->GetSolid()->GetName();

// G4cout << localPos << G4endl;



if ((volumeName =="histone")||(volumeName=="world")||(volumeName=="waterBox"))

{
// G4cout<<"Removed " <<GetMolecule(track)->GetName() <<track->GetTrackID() << " " << track->GetPosition()<<" " << volumeName << G4endl;//<< " "<<track->GetCreatorProcess()->GetProcessName()<<G4endl;

  // Histones act as perfect scavengers for all radiolysis products
  // Do not want to track radiolysis products in the other volumes to save computational time, this leaves volumes only immediately around the bases. more than 9nm from literature
  track->SetTrackStatus(fKillTrackAndSecondaries);

}
delete newTouchable;
}
}


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::UserReactionAction(const G4Track& trackA,
    const G4Track& trackB,
    const std::vector<G4Track*>* pProducts)
{
      //check for the case "no product"
      // G4cout << GetMolecule(trackA)->GetName() <<"+" << GetMolecule(trackB)->GetName() << G4endl;
    // if(!pProducts)
    // {
    //     return;
    // }

    //   if(pProducts && (GetMolecule((*pProducts)[0])->
    // GetDefinition() != G4DamagedDeoxyribose::Definition()))
    // {
    //     return;
    // }
    
  if ((GetMolecule((&trackA))->GetDefinition() == G4Deoxyribose::Definition())||(GetMolecule((&trackB))->GetDefinition() == G4Deoxyribose::Definition()))
  {
    
    const G4Track* DNAElement = nullptr;
    const G4Track* radical    = nullptr;
    if(GetMolecule(&trackA)->
    GetDefinition() == G4Deoxyribose::Definition())
    {
        DNAElement = &trackA;
        radical    = &trackB;
    }
    else 
    {
        DNAElement = &trackB;
        radical    = &trackA;
    }
    
    if(GetMolecule(radical)->GetDefinition() != G4OH::Definition())
    {
        return;
    }

    CommandLineParser *parser = CommandLineParser::GetParser();
    Command *command(0);
    if ((command = parser->GetCommandIfActive("-out")) == 0)
      return;


    G4ThreeVector localPos = DNAElement->GetPosition();

    G4VTouchable* newTouchable = CreateTouchableAtPoint(localPos);

    G4String volumeName = newTouchable->GetVolume()->GetName();
    G4int copyNum =newTouchable->GetReplicaNumber();

   delete newTouchable;


    G4int flagVolume = 0;

    if (volumeName=="sugar0")
    {
     flagVolume = 1; 
    }
    else if (volumeName=="sugar1")
    {
     flagVolume = 2; 
    }

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    G4int eventID =  G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    analysisManager->FillNtupleIColumn(2, 0, eventID);
    analysisManager->FillNtupleDColumn(2, 1, DNAElement->GetStep()->GetPreStepPoint()->GetPosition().x() / nanometer);
    analysisManager->FillNtupleDColumn(2, 2, DNAElement->GetStep()->GetPreStepPoint()->GetPosition().y() / nanometer);
    analysisManager->FillNtupleDColumn(2, 3, DNAElement->GetStep()->GetPreStepPoint()->GetPosition().z() / nanometer);
    analysisManager->FillNtupleIColumn(2, 4, flagVolume);
    analysisManager->FillNtupleIColumn(2, 5, copyNum);


    analysisManager->AddNtupleRow(2);
  }

}


G4Navigator* TimeStepAction::GetNavigator() {
  static G4ThreadLocal G4Navigator* theNavigator = 0;
  if (!theNavigator) theNavigator = new G4Navigator;

  // Make sure current world volume is the one in use  
  G4VPhysicalVolume* theWorld =
    G4TransportationManager::GetTransportationManager()->
      GetNavigatorForTracking()->GetWorldVolume();

  if (theNavigator->GetWorldVolume() != theWorld)
    theNavigator->SetWorldVolume(theWorld);

  return theNavigator;
}

G4VTouchable* TimeStepAction::CreateTouchableAtPoint(const G4ThreeVector& pos) {
  G4VTouchable* touchable = new G4TouchableHistory;
  GetNavigator()->LocateGlobalPointAndUpdateTouchable(pos, touchable, false);
  return touchable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
