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
/// \file StackingAction.cc
/// \brief Implementation of the StackingAction class

#include "StackingAction.hh"
#include "G4DNAChemistryManager.hh"
#include "G4StackManager.hh"
#include "G4ITTransportationManager.hh"
#include "G4ITTrackHolder.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
    : G4UserStackingAction()
{
    fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::NewStage()
{

    if (stackManager->GetNTotalTrack() == 0)
    {
        if (fpEventAction->GetEdepTracking() > 0)
        {
            G4ITTransportationManager::GetTransportationManager()->SetWorldForTracking(
                G4ITTransportationManager::GetTransportationManager()->GetParallelWorld("ChemistryWorld"));
            G4DNAChemistryManager::Instance()->Run(); // starts chemistry
        }
        else
        {
            auto &fDelayedList = G4ITTrackHolder::Instance()->GetDelayedLists();

            if (!fDelayedList.empty())
            {
                fDelayedList.clear();
            }

        }
    }
}
