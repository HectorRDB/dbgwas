//
// Created by Leandro Ishi Soares de Lima on 04/02/17.
//

#include "global.h"
#include <string>

using namespace std;

const char* STR_STRAINS_FILE = "-strains";
const char* STR_MAX_NEIGHBOURHOOD = "-nh";
const char* STR_KSKMER_SIZE = "-k";

//TODO: we should put this back. I put it out and forced the output folder to be always ./output because gemma forcibly uses this directory. If there are 2 executions of the tool, this could bugs because of gemma
//const char* STR_OUTPUT = "-output";
//TODO: we should put this back. I put it out and forced the output folder to be always ./output because gemma forcibly uses this directory. If there are 2 executions of the tool, this could bugs because of gemma

const char* STR_NBCORES = "-nb-cores";
const char* STR_SKIP1 = "-skip1";
const char* STR_SKIP2 = "-skip2";
const char* STR_NEWICK_PATH = "-newick";
const char* STR_SFF = "-SFF";

//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
//const char* STR_COUNT_MODE = "-count";
//TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option

string pathToExecParent = "";
bool skip1 = false;
bool skip2 = false;
bool presenceAbsenceCountMode = false;
boost::variant< int, double > SFF;

//global vars used by both programs
Graph* graph;
vector< UnitigIdStrandPos >* nodeIdToUnitigId;
vector< Strain >* strains = NULL;

void populateParser (Tool *tool) {
  // We add some custom arguments for command line interface
  //TODO: we should put this back. I put it out and forced the output folder to be always ./output because gemma forcibly uses this directory. If there are 2 executions of the tool, this could bugs because of gemma
  //tool->getParser()->push_front (new OptionOneParam (STR_OUTPUT, "Path to the folder where the final and temporary files will be stored",  false, "output"));
  //TODO: we should put this back. I put it out and forced the output folder to be always ./output because gemma forcibly uses this directory. If there are 2 executions of the tool, this could bugs because of gemma

  tool->getParser()->push_front (new OptionOneParam (STR_KSKMER_SIZE, "K-mer size",  false, "31"));
  tool->getParser()->push_front (new OptionOneParam (STR_STRAINS_FILE, "A text file describing the strains containing 3 collumns: 1) ID of the strain; 2) Phenotype (0/1/NA); 3) Path to a multi-fasta file containing the sequences of the strain. This file needs a header. Check the sample_example folder for an example.",  true));
  tool->getParser()->push_front (new OptionOneParam (STR_NEWICK_PATH, "Path to a newick tree file",  true));
  tool->getParser()->push_front (new OptionOneParam (STR_MAX_NEIGHBOURHOOD, "Denotes the maximum neighbourhood that can be viewed in the final visualization",  false, "5"));

  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
  //tool->getParser()->push_front (new OptionOneParam (STR_COUNT_MODE, "The count mode. If \"01\", then the count mode is seen as presence/absence. If \"Freq\", then the count mode is seen as the frequency",  false, "01"));
  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option

  tool->getParser()->push_front (new OptionOneParam (STR_SFF, "Denotes the Significant Features Filter - the features selected to create a visualisation around them. If it is a float number n, then only the features with q-value<=n are selected. If it is an integer, then only the n first features are selected. Take a look at the output/significant_patterns.txt file to get a list of features ordered by q-value to better choose this parameter (re-run the tool with -skip2 in order to directly produce the visualisation of the features selected by your parameter).",  false, "100"));

  tool->getParser()->push_front (new OptionNoParam (STR_SKIP1, "Skips Step 1, running only Steps 2 and 3. Assumes that in the output folder, we have the files bugwas_input*, frequency_unitig_to_total_pheno0_pheno1_NA_count, gemma_input*, graph.edges.dbg, graph.h5, and graph.nodes",  false));
  tool->getParser()->push_front (new OptionNoParam (STR_SKIP2, "Skips Steps 1 and 2, running only Step 3. Assumes that in the output folder, we have the files graph.nodes, graph.edges.dbg, frequency_unitig_to_total_pheno0_pheno1_NA_count, significant_unitigs",  false));

}
