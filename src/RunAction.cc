//
#include "RunAction.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
//
#include "DetectorConstruction.hh"
#include "G4VSensitiveDetector.hh"
#include <math.h>
#include "TH1D.h"
#include "TFile.h"
#include <sstream>

RunAction::RunAction(DetectorConstruction* DC){
  //
  
  //take the DetectorConstruction pointer given when this object is created (in main) and copy to local member
  myDC = DC;
  //
  
}

RunAction::~RunAction(){

}

void RunAction::BeginOfRunAction(const G4Run* aRun){
  
  //Arguments:Name of file,type of rewrite option(UPDATE means append),comments
  RootOP = new TFile("EdepNoQuench.root","RECREATE","LS Sim Output");

  //arguments: Name, Title, number of bins, lower level range, upper
  //level range.
  //
  //All Energies are passed in MeV.
  //
  EdepNoQuench = new TH1D("Edep with no Quench","Energy Deposition in Scintillant per Event, no quench correction",210000,0.,210000.);
  //
  //EdepQuenchPerEvt = new TH1D("Edep with Quench per Evt","Energy Deposition in Scintillant per Event, quench correction applied per event",100000,0.,1000.);
  //
  //EdepQuenchPerInt = new TH1D("Edep with Quench per Int","Energy Deposition in Scintillant per Event, quench correction applied per Interaction",100000,0.,1000.);
  EdepRear = new TH1D("Edep with no Quench, rear vial","Energy Deposition in rear Scintillant per Event, no quench correction",210000,0.,210000.);

  //EdepTreeRear = new TTree("ResultsRear",
  //			   "Particles causing Edep in rear vial");
  //EdepTreeFront = new TTree("ResultsFront",
  //			    "Particles causing Edep in front vial");

  EdepTree = new TTree("Results", "Particles causing Edep in either vial");

  EdepTree->Branch("TypeFront", &IntTypeFront);
  EdepTree->Branch("KinEnFront", &KinEnInFront);
  EdepTree->Branch("EdepFront", &totalEdepFront);
  EdepTree->Branch("TrLenPriorFront", &TrackLenPriorFront);
  EdepTree->Branch("TrLenInVolFront", &TrackLenInVolFront);
  EdepTree->Branch("VertXFront", &VertXFront);
  EdepTree->Branch("VertZFront", &VertZFront);
  EdepTree->Branch("VertRFront", &VertRFront);
  EdepTree->Branch("TypeRear", &IntTypeRear);
  EdepTree->Branch("KinEnRear", &KinEnInRear);
  EdepTree->Branch("EdepRear", &totalEdepRear);
  EdepTree->Branch("TrLenPriorRear", &TrackLenPriorRear);
  EdepTree->Branch("TrLenInVolRear", &TrackLenInVolRear);
  EdepTree->Branch("VertXRear", &VertXRear);
  EdepTree->Branch("VertZRear", &VertZRear);
  EdepTree->Branch("VertRRear", &VertRRear);

  HDPETree = new TTree("HDPEResults", "Protons stopping in the HDPE");

  //HDPETree->Branch("Type", &ParticleType);
  HDPETree->Branch("VertX", &VertX);
  HDPETree->Branch("VertZ", &VertZ);
  HDPETree->Branch("StopX", &StopX);
  HDPETree->Branch("StopY", &StopY);
  HDPETree->Branch("StopZ", &StopZ);
  HDPETree->Branch("StopR", &StopR);

}



void RunAction::TallyEdepNoQuench(const G4double thisEdep){

  //G4cout << "Inside TallyEdepNoQuench, Edep = " << thisEdep << G4endl;

  G4double  thisEdepkeV = thisEdep*1000;

  //Fill ROOT histogram
  EdepNoQuench->Fill(thisEdepkeV);

}

void RunAction::TallyHDPE(G4double& theVertX, G4double& theVertZ,
			  G4double& theStopX, G4double& theStopZ,
			  G4double& theStopY){
  
  VertX = theVertX;
  VertZ = theVertZ;
  StopX = theStopX;
  StopY = theStopY;
  StopZ = theStopZ;
  StopR = sqrt((StopX - VertX)*(StopX - VertX) + 
	       (StopZ - VertZ)*(StopZ - VertZ));


  //Fill ROOT histogram
  HDPETree->Fill();

}


void RunAction::TallyEvtData(G4int typeFr, G4double KinEnFr, G4double EdepFr,
			     G4double TrLenPriorFr, G4double TrLenInFr,
			     G4double VxFr, G4double VzFr, G4double VrFr,
			     G4int typeR, G4double KinEnR, G4double EdepR,
			     G4double TrLenPriorR, G4double TrLenInR,
			     G4double VxR, G4double VzR, G4double VrR){
  IntTypeFront = typeFr;
  KinEnInFront = KinEnFr;
  totalEdepFront = EdepFr;
  TrackLenPriorFront = TrLenPriorFr;
  TrackLenInVolFront = TrLenInFr;
  VertXFront = VxFr;
  VertZFront = VzFr;
  VertRFront = VrFr;
  IntTypeRear = typeR;
  KinEnInRear = KinEnR;
  totalEdepRear = EdepR;
  TrackLenPriorRear = TrLenPriorR;
  TrackLenInVolRear = TrLenInR;
  VertXRear = VxR;
  VertZRear = VzR;
  VertRRear = VrR;

  EdepTree->Fill();

}

//task to be carried out at the end of the run
void RunAction::EndOfRunAction(const G4Run* aRun){
  //get the number of primary particles being simulated for this run
  G4double NumberOfEvents = aRun->GetNumberOfEventToBeProcessed();
 
  // Name the histograms
  G4String SpectName;
  G4String Beginning = "EdepNoQuenchFront";
  //stringstream s;
  //G4double RunNumber = aRun->GetRunID();
  //s << RunNumber;
  //G4String RunNum;
  //s >> RunNum;
  //SpectName = Beginning+RunNum;
  
  EdepNoQuench->SetName(Beginning);
  EdepNoQuench->Scale(1/NumberOfEvents);
  EdepNoQuench->Write();

  Beginning = "EdepNoQuenchRear";
  //SpectName = Beginning+RunNum;
  EdepRear->SetName(Beginning);
  EdepRear->Scale(1/NumberOfEvents);
  EdepRear->Write();
 
  //Beginning = "EdepQuenchPerEvtFromRunNum";
  //SpectName = Beginning+RunNum;
  //EdepQuenchPerEvt->SetName(SpectName);
  //EdepQuenchPerEvt->Scale(1/NumberOfEvents);
  //EdepQuenchPerEvt->Write();

  // Beginning = "EdepQuenchPerIntFromRunNum";
  //SpectName = Beginning+RunNum;
  //EdepQuenchPerInt->SetName(SpectName);
  //EdepQuenchPerInt->Scale(1/NumberOfEvents);
  //EdepQuenchPerInt->Write();
  RootOP->Write();
  RootOP->Close();
  

}
