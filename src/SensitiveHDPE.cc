#include "SensitiveHDPE.hh"
#include "G4Step.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4SteppingManager.hh"
#include <iterator>
#include "G4TrackVector.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

SensitiveHDPE::SensitiveHDPE(G4String name) : G4VSensitiveDetector(name){
  //runnum = 0;
  name = "HDPE_nistabs_log";
}

SensitiveHDPE::~SensitiveHDPE(){;}

/*This method is invoked at the beginning of each event. The argument of this method is an object of the G4HCofThisEvent class. Hits collections, where hits produced in this particular event are stored, can be associated to the G4HCofThisEvent object in this method. The hits collections associated with the G4HCofThisEvent  object during this method can be used for ``during the event processing'' digitization.*/
void SensitiveHDPE::Initialize(G4HCofThisEvent* HCE){
  runnum = 0;
  TrackID = 0;
  ParentID = 0;
  Counter = 0;
  EdepThisEventUnquenched = 0;


  MultipleIn = false;
  KinEnIn = 0.;
  EvtType = 0;
  TrackLenInVol = 0.;
  TrackLenPrior = 0.;
  VertexX = 0.;
  VertexZ = 0.;
  VertexR = 0.;

  pVertX = 0;
  pVertZ = 0;
  pStopX = 0;
  pStopY = 0;
  pStopZ = 0;

}

G4int SensitiveHDPE::GetEvtType(){return EvtType;}
G4double SensitiveHDPE::GetKinEnIn(){return KinEnIn;}
G4double SensitiveHDPE::GetTrackLenInVol(){return TrackLenInVol;}
G4double SensitiveHDPE::GetTrackLenPrior(){return TrackLenPrior;}
G4double SensitiveHDPE::GetVertexX(){return VertexX;}
G4double SensitiveHDPE::GetVertexZ(){return VertexZ;}
G4double SensitiveHDPE::GetVertexR(){return VertexR;}
G4double SensitiveHDPE::GetEdep(){return EdepThisEventUnquenched;}

/*This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. The first argument of this method is a G4Step  object of the current step. The second argument is a G4TouchableHistory object for the ``Readout geometry'' described in the next section. The second argument is NULL for the case ``Readout geometry'' is not assigned to this sensitive detector. In this method, one or more G4VHit objects should be constructed if the current step is meaningful for your detector.*/
G4bool SensitiveHDPE::ProcessHits(G4Step* theStep, G4TouchableHistory*){

  //Here's the plan: I'll classify events into different categories:
  //1. Multiple particles entered the vial and deposited energy.
  //2. Only the primary proton entered the vial and deposited energy.
  //3. Only a gamma ray entered the vial and deposited energy.
  //4. Only a secondary proton entered the vial and deposited energy.
  //5. Only an electron entered the vial and deposited energy.
  //6. Only a neutron entered the vial and deposited energy.
  //7. Only an some other particle entered the vial and deposited energy.

  //For each of these, I'll save the total energy deposit in the vial for that
  //event and the kinetic energy of the incident particle (or the sum if there
  //were multiple).

  //Check whether the particle is stepping into the volume from outside...
  if((theStep->GetTrack()->GetTrackID()!=TrackID)&&
     (theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()!="HDPE_nistabs_log")
     ){
    //stepping into material for the first time
    if(TrackID!=0){
      MultipleIn = true;
      EvtType = 1;
    }      
    
    TrackID = theStep->GetTrack()->GetTrackID();
    KinEnIn += theStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
    //Collect the vertex locations of the entering particle.
    VertexX = theStep->GetTrack()->GetVertexPosition().getX();
    VertexZ = theStep->GetTrack()->GetVertexPosition().getZ();
    VertexR = sqrt(VertexX*VertexX + VertexZ*VertexZ);


    TrackLenPrior += theStep->GetTrack()->GetTrackLength();

    string thename = (theStep->GetPreStepPoint()->GetMaterial()->GetName());
    //G4cout << "Entering Material " << thename << endl;
    //G4cout << "Track ID = " << theStep->GetTrack()->GetTrackID() << endl;
    //thename = (theStep->GetTrack()->GetDefinition()->GetParticleName());
    //G4cout << "Particle name = " << thename << endl;
    //G4cout << "Particle Energy = "
    //	   << theStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()
    //	   << " Mev" << endl;
    //G4cout << "Logical vol at vertex = "
    //	   << theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()
    //	   << endl;
    if(theStep->GetTrack()->GetTrackID() != 1){
    //thename = theStep->GetTrack()->GetCreatorProcess()->GetProcessName();
    //G4cout << "Creator process was " << thename << endl;
      thename = (theStep->GetTrack()->GetDefinition()->GetParticleName());
      if(thename == "gamma"){
	EvtType = 3;
      }
      else if(thename == "proton"){
	EvtType = 4;
      }
      else if(thename == "e-"){
	EvtType = 5;
      }
      else if(thename == "neutron"){
	EvtType = 6;
      }
      else{
	EvtType = 7;
      }
    }
    else{
      //G4cout << "PRIMARY PARTICLE" << endl;
      EvtType = 2;
    }

  }

  //Now accumulate the path length in the scintillator.
  if((theStep->GetTrack()->GetTrackID()==TrackID)&&
     (theStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()!="HDPE_nistabs_log")
     ){
    //The particle that we're tracking, get the step length.
    TrackLenInVol+=theStep->GetStepLength();
  }

 
  //Need to add in alternate Edep collection, simple Edep per step.
  EdepThisEventUnquenched += theStep->GetTotalEnergyDeposit()*MeV;
  //G4cout << "Energy deposit so far this event = " << EdepThisEventUnquenched
  //	 << G4endl;

  //Check if last step.
  if(theStep->GetPostStepPoint()->GetKineticEnergy()==0){
    if(theStep->GetTrack()->GetParticleDefinition()->GetParticleName() 
       == "proton"){
      //cout << "A Proton has Stopped" << endl;
      //get run action pointer
      RunAction* myRunAction = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
      pVertX = theStep->GetTrack()->GetVertexPosition().getX();
      pVertZ = theStep->GetTrack()->GetVertexPosition().getZ();
      pStopX = theStep->GetPostStepPoint()->GetPosition().getX();
      pStopZ = theStep->GetPostStepPoint()->GetPosition().getZ();
      pStopY = theStep->GetPostStepPoint()->GetPosition().getY();

      myRunAction->TallyHDPE(pVertX, pVertZ, pStopX, pStopZ, pStopY);

    }
  }

  return true;  
}


/*This method is invoked at the end of each event. The argument of this method is the same object as the previous method. Hits collections occasionally created in your sensitive detector can be associated to the G4HCofThisEvent object.*/
void SensitiveHDPE::EndOfEvent(G4HCofThisEvent*)
{
  //only need to do anything if there were some interactions.
  if(EdepThisEventUnquenched!=0){
    //get run action pointer
    RunAction* myRunAction = (RunAction*)(G4RunManager::GetRunManager()->GetUserRunAction());
    
    if(myRunAction){
      //Tally unquenched Edeps.
      //G4cout << "Calling RunAction->TallyEdepNoQuench" << G4endl;
      myRunAction->TallyEdepNoQuench(EdepThisEventUnquenched);

    }
  }
}
