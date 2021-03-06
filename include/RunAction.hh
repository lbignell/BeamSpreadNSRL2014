#ifndef RunAction_hh
#define RunAction_hh 1
//
#include "G4UserRunAction.hh"
#include "G4UnitsTable.hh"
#include "DetectorConstruction.hh"

#include "TH1D.h"
#include <fstream>
#include "TTree.h"

//declare the DetectorConstruction class as we will define a pointer later
class DetectorConstruction;
class TFile;

//needed for using standard libraries
using namespace std;

//run action class, carries out tasks at the begin and end of each run
//the concept of a run incorporates a fixed geometry, fixed beam conditions, simulation of number of primaries
//begins with /run/beamOn command and finishes with tracking of last secondary to zero energy 

class RunAction : public G4UserRunAction {
//
public:
//run action class needs pointer ot the detector construction class in order to get details of the readout geometry
//accepts pointer to detector construction class

  RunAction(DetectorConstruction*);
  ~RunAction();

    //

private:
//an ofstream to access the output file
  ofstream outfile;

  TH1D *EdepNoQuench; TH1D* EdepRear;

  TTree* EdepTree;
  TTree* HDPETree;
  //TBranch* TypeBranch;
  //TBranch* KinEnInBranch;
  //TBranch* EdepBranch;
  G4int IntTypeFront;
  G4double KinEnInFront;
  G4double totalEdepFront;
  G4double TrackLenPriorFront;
  G4double TrackLenInVolFront;
  G4double VertXFront;
  G4double VertZFront;
  G4double VertRFront;
  G4int IntTypeRear;
  G4double KinEnInRear;
  G4double totalEdepRear;
  G4double TrackLenPriorRear;
  G4double TrackLenInVolRear;
  G4double VertXRear;
  G4double VertZRear;
  G4double VertRRear;
  TFile* RootOP;

  G4double VertX;
  G4double VertZ;
  G4double StopX;
  G4double StopY;
  G4double StopZ;
  G4double StopR;



//local pointer for detector construction class
    DetectorConstruction* myDC;

public:

//note argument of these methods is a pointer to a G4Run object
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void TallyEdepNoQuench(const G4double);
  //void TallyEdepQuenchPerEvt(const G4double);
  //void TallyEdepQuenchPerInt(const G4double);
  //Arguments; vertex X, vertex Z, StopX, StopZ, StopY, Edep in HDPE.
  void TallyHDPE(G4double&, G4double&, G4double&, G4double&,
		 G4double&);
  //Argument list: Interaction type, Kinetic Energy in, Total Edep, Track length
  //prior to entering liquid, track length in liquid.
  //void TallyEvtDataRear(G4int, G4double, G4double, G4double, G4double);
  //void TallyEvtDataFront(G4int, G4double, G4double, G4double, G4double);
  //Arg list: Interaction type, Kinetic Energy in, Total Edep, Track length
  //prior to entering liquid, track length in liquid, vertex X, vertex Z, vertex
  //R. The first lot is for the front detector, the second is for the rear
  //detector.
  //Note to self: I should really modify this so that I'm passing a pointer to a
  //vector or something similar...
  void TallyEvtData(G4int, G4double, G4double, G4double, G4double,
		    G4double, G4double, G4double,
		    G4int, G4double, G4double, G4double, G4double,
		    G4double, G4double, G4double);

};

#endif
