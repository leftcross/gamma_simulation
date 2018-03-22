#include <fstream>
#include <iostream>
#include "TMap.h"
#include "TObjString.h"
#include "TTree.h"
#include "TList.h"
#include "TParameter.h"
#include "TNamed.h"
#include "TString.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>


void Macros(TString fName="muon",Int_t num = 20,double numberRun=3000){

  std::cout<<"Run started \n";
  double halfx=100,halfy=100,halfz=1;
  double centrex=0,centrey=0,centrez=198;
  TString rootfile="muon_cosmic";

  double labdepth=6.01; // Depth of the lab below ground
  
  srand (time(NULL));

      for (int k=0+10; k < num+10; k++){
        
        std::cout<<"Iteration "<<k<<std::endl;
	stringstream ss;
	ss<<(k+1);
	TString strn=ss.str();
        TString FileNameIn = fName+"_"+strn+".mac";
        std::fstream macroFile(FileNameIn,std::fstream::out);

    

        macroFile << "# verbosities\n";
        macroFile << "/control/verbose 0 \n";
        macroFile << "/run/verbose 0\n";
        macroFile << "/tracking/verbose 0 \n";

        macroFile << "/random/setSeeds" << " " << rand()%200 << " " << rand()%300 << endl;

        macroFile <<"#\n";
        macroFile <<"# draw flags\n";
        macroFile <<"#/damic/draw/drawColours custom\n";
        macroFile <<"#/damic/draw/gammaColour yellow\n";
        macroFile <<"#/damic/draw/drawTracks all\n";
        macroFile <<"#/damic/draw/drawHits 1\n";
        macroFile <<"#\n";

        macroFile  <<"# file output\n";
        macroFile  <<"/damic/saveHits 0\n";
        macroFile  <<"#\n\n";

        macroFile  <<"# gun: particle\n";
        macroFile  <<"/damic/gun/verbose 0\n";
        macroFile  <<"/damic/gun/energytype Muon\n";
        macroFile  <<"/damic/gun/particle  mu-\n";
        macroFile  <<"/damic/gun/depth "<<labdepth<<" km\n";
        macroFile  <<"# gun: distribution\n";
        macroFile  <<"/damic/gun/angtype Muon\n\n";

        macroFile  <<"# gun: source\n";
        macroFile  <<"/damic/gun/type Volume\n";
        macroFile  <<"/damic/gun/shape Cube\n";
        macroFile  <<"/damic/gun/halfx "<<  halfx <<" cm\n";
        macroFile  <<"/damic/gun/halfy "<<  halfy <<" cm\n";
        macroFile  <<"/damic/gun/halfz "<<  halfz <<" cm\n";
        macroFile  <<"/damic/gun/centre "<<centrex<<" "<<centrey<<" "<<centrez<<" "<<" cm\n\n";

        macroFile  <<"/damic/histogramfile "<<rootfile<<"_"<<(k+1)<<std::endl;

        macroFile  <<"/run/beamOn "<<numberRun<<std::endl;
   
        macroFile.close();

      }
      std::cout<<"Done\n";
    }
