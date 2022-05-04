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


#include "RunAction.hh"
#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "globals.hh"
#include <map>
#include "CommandLineParser.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
using namespace G4DNAPARSER;

RunAction::RunAction()
    : G4UserRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run*)
{
    CreateNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run*)
{
    WriteNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::CreateNtuple()
{
      CommandLineParser *parser = CommandLineParser::GetParser();
        Command *command(0);
        if ((command = parser->GetCommandIfActive("-out")) == 0)
            return;

      // Open an output file
  G4String fileName {"output.root"};
  if (command->GetOption().empty() == false)
  {
    fileName = command->GetOption();
  }

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetVerboseLevel(0);
    analysisManager->SetNtupleDirectoryName("ntuple");

    // open output file
    //
    G4bool fileOpen = analysisManager->OpenFile(fileName);
    if (!fileOpen)
    {
        G4cout << "\n---> HistoManager::book(): cannot open " << fileName<< G4endl;
        return;
    }

    analysisManager->SetFirstNtupleId(0);
    analysisManager->CreateNtuple("EventEdep", "EventEdep");
  analysisManager->CreateNtupleDColumn("Edep_J");
  analysisManager->CreateNtupleIColumn("EventNo");
  analysisManager->CreateNtupleDColumn("TrackLengthChomatin");

  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Direct", "Direct");
  analysisManager->CreateNtupleIColumn(1, "EventNo");
  analysisManager->CreateNtupleIColumn(1, "copyNum");
  analysisManager->CreateNtupleDColumn(1, "eDep_eV");
  analysisManager->CreateNtupleDColumn(1, "x");
  analysisManager->CreateNtupleDColumn(1, "y");
  analysisManager->CreateNtupleDColumn(1, "z");
  analysisManager->CreateNtupleIColumn(1, "Strand");
    analysisManager->FinishNtuple(1);

    // For chemistry

  analysisManager->CreateNtuple("Indirect", "Indirect");

  analysisManager->CreateNtupleIColumn(2, "EventNo");
  analysisManager->CreateNtupleDColumn(2, "x");
  analysisManager->CreateNtupleDColumn(2, "y");
  analysisManager->CreateNtupleDColumn(2, "z");
  analysisManager->CreateNtupleIColumn(2, "Strand");
  analysisManager->CreateNtupleIColumn(2, "copyNum");


    analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("Events", "Events");

  analysisManager->CreateNtupleIColumn(3, "EventNo");
  analysisManager->CreateNtupleDColumn(3, "Energy");
  analysisManager->CreateNtupleDColumn(3, "Origin_x");
  analysisManager->CreateNtupleDColumn(3, "Origin_y");
  analysisManager->CreateNtupleDColumn(3, "Origin_z");
  analysisManager->CreateNtupleDColumn(3, "Direction_x");
  analysisManager->CreateNtupleDColumn(3, "Direction_y");
  analysisManager->CreateNtupleDColumn(3, "Direction_z");

  analysisManager->FinishNtuple(3);

    G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void RunAction::WriteNtuple()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    analysisManager->Clear();
    G4cout << "\n----> Histograms are saved" << G4endl;
}
