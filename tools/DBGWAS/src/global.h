//
// Created by Leandro Ishi Soares de Lima on 19/07/16.
//

#ifndef KSGATB_GLOBAL_H
#define KSGATB_GLOBAL_H
#include <gatb/gatb_core.hpp>
#include "Utils.h"

//global vars
extern Graph* graph;
extern vector< UnitigIdStrandPos >* nodeIdToUnitigId;
extern vector< Strain >* strains;
extern const char* STR_STRAINS_FILE;
extern const char* STR_KSKMER_SIZE;

//TODO: we should put this back. I put it out and forced the output folder to be always ./output because gemma forcibly uses this directory. If there are 2 executions of the tool, this could bugs because of gemma
//extern const char* STR_OUTPUT;
//TODO: we should put this back. I put it out and forced the output folder to be always ./output because gemma forcibly uses this directory. If there are 2 executions of the tool, this could bugs because of gemma

extern const char* STR_NBCORES;
extern const char* STR_MAX_NEIGHBOURHOOD;
extern const char* STR_SKIP1;
extern const char* STR_SKIP2;
extern const char* STR_NEWICK_PATH;
extern const char* STR_SFF;

//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
//extern const char* STR_COUNT_MODE;
//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option

extern string pathToExecParent;
extern bool skip1;
extern bool skip2;
extern bool presenceAbsenceCountMode;
extern boost::variant< int, double > SFF;

void populateParser (Tool *tool);

#endif //KSGATB_GLOBAL_H
