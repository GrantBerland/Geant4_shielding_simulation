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
// $Id: EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"


#include <fstream>


EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fEdep_det1(0.),
  fEdep_det2(0.),
  fEdep_det3(0.),
  fEdep_det4(0.)
{}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
  fEdep_det1 = 0.;
  fEdep_det2 = 0.;
  fEdep_det3 = 0.;
  fEdep_det4 = 0.;
}

void EventAction::AddEdep(G4double edep)
{
  edep += fEdep;
}

void EventAction::AddEdep_multiple(G4String solid, G4double edep)
{
/*
  if (solid == "Electronics") {fEdep += edep;}
  if (solid == "Shielding") {fEdep_det1 += edep;}
  if (solid == "Electronics") {fEdep_det2 += edep;}
  if (solid == "Electronics") {fEdep_det3 += edep;}
  if (solid == "Electronics") {fEdep_det4 += edep;}
*/
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  /*
  fRunAction->AddEdep(fEdep);

  G4AnalysisManager* man = G4AnalysisManager::Instance();

  G4double init_energy = event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();

  man->FillH1(1, init_energy);
  man->FillH1(2, fEdep);
  man->FillH1(3, fEdep_det1);
  man->FillH1(4, fEdep_det2);
  man->FillH1(5, fEdep_det3);
  man->FillH1(6, fEdep_det4);
*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
