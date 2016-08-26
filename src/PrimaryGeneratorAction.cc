//once again, fairly simple - copied off the UOW tute
#include "PrimaryGeneratorAction.hh"

//include class definition for the particle gun
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

//include header file for particle 
#include "G4Proton.hh"
//#include "G4Electron.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include <math.h>

PrimaryGeneratorAction::PrimaryGeneratorAction(){
//we need to use dynamic memory allocation for the G4ParticleGun object
  gun = new G4ParticleGun(1);
  gun->SetParticleDefinition(G4Proton::ProtonDefinition());
  gun->SetParticleEnergy(201.*MeV);
  
  //X will be the lateral dimension, and Z the vertical
  //G4double X = -10*cm+20*cm*G4UniformRand();
  //G4double Z = -10*cm+20*cm*G4UniformRand();

  //Assume beam window is 30 cm from the first vial's centre
  //gun->SetParticlePosition(G4ThreeVector(X,-30*cm,Z));
  gun->SetParticlePosition(G4ThreeVector(0, -30*cm, -10*mm));

  //hopefully G4UniformRand gives a new random number each call.

  gun->SetParticleMomentumDirection(G4ThreeVector(0,1,0));
  //momentum directional cosines


}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
//free the dynamically allocated memory
  delete gun;

}

//this method will be called by the RunManager at the beginning of each particle history
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  
  //X will be the lateral dimension, and Z the vertical
  //G4double X = -10*cm+20*cm*G4UniformRand();
  //G4double Z = -10*cm+20*cm*G4UniformRand();

  //G4cout << "Vertex X = " << X << ", Vertex Z = " << Z << G4endl;

  //Assume beam window is 30 cm from the first vial's centre
  //gun->SetParticlePosition(G4ThreeVector(X,-30*cm,Z));
  gun->SetParticlePosition(G4ThreeVector(0, -30*cm, -10*mm));

  //hopefully G4UniformRand gives a new random number each call.

  gun->SetParticleMomentumDirection(G4ThreeVector(0,1,0));
  //momentum directional cosines

  //this gets the particle firing
  gun->GeneratePrimaryVertex(anEvent);

}
