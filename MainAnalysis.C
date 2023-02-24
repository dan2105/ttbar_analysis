#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TProof.h"
#include <string>
#include "TString.h"

void MainAnalysis( int option = 0){

//Vou adicionar o proof aqui e vou colocar a leitura do arquivo daqui TString path = .....
    TString path = "/Users/danielernani/Desktop/workspace/MyTopAnalysis/";

    if(option==11){
        TChain* ttbar_inclusive = new TChain("nominal");
        ttbar_inclusive->AddFile(path+"data/out_data17_lowmu.root");
        ttbar_inclusive->Process(path+"source/top_analysis_run2.C+","top_data17_13TeV_lowmu");
     }

}
