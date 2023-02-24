/*
Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "TopCommonDilepton/TopCommonDileptonEventSaver.h"
#include "TopCommonDilepton/NeutrinoWeighter.h"
#include "TopCommonDilepton/TopDileptonReconstruction.h"
#include "TopCommonDilepton/TTbarNNLORecursiveRew.h"
#include "TopCommonDilepton/TTbarNNLOReweighter.h"

#ifdef __CINT__

#pragma extra_include "TopCommonDilepton/TopCommonDileptonEventSaver.h";
#pragma extra_include "TopCommonDilepton/NeutrinoWeighter.h";
#pragma extra_include "TopCommonDilepton/TopDileptonReconstruction.h";

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#pragma link C++ class NeutrinoWeighter+;
#pragma link C++ class TopCommonDileptonEventSaver+;
#pragma link C++ class TopDileptonReconstruction+;
#pragma link C++ class TTbarNNLORecursiveRew+;
#pragma link C++ class TTbarNNLOReweighter+;

//#pragma link C++ class top::CustomObjectLoader+;
//#pragma link C++ class top::HowtoExtendAnalysisTopLoader+;
//#pragma link C++ class top::CustomEventSaver+;

#endif
