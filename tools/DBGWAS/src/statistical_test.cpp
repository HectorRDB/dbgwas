//
// Created by Leandro Ishi Soares de Lima on 10/07/17.
//

#include "statistical_test.h"
#include "global.h"

// The Tool constructor allows to give a name to our tool.
// This name appears when one gets help in the command line or in the final output
statistical_test::statistical_test ()  : Tool ("statistical_test") //give a name to our tool
{
  populateParser(this);
}

void statistical_test::execute () {
  //parameters
  checkParametersStatisticalTest(this);
  if (skip2) return;
  //string outputFolder = getInput()->getStr(STR_OUTPUT);
  string outputFolder("output");
  string newickTreeFilePath = getInput()->getStr(STR_NEWICK_PATH);

  //execute the statistical test
  //To do so, we need to create an output folder inside the outputFolder
  createFolder(string("output/output"));
  //We also need to create this other folder... bizzare...
  createFolder(string("output/bugwas_out__PCloadings/output"));

  //execute the statistical test itself...
  executeCommand(string("Rscript --vanilla DBGWAS.R ") + outputFolder +
                 " " + outputFolder + "/bugwas_input.id_phenotype " +
                 newickTreeFilePath + " " + outputFolder +"/bugwas_out_ " +
                     pathToExecParent + "gemma.0.93b 0.01");





  //sort the file by q-value and output it to output/patterns.txt
  //read
  auto patterns = PatternFromStats::readFile(outputFolder + "/bugwas_out__DBGWAS_patterns.txt");
  //sort
  sort(patterns.begin(), patterns.end());
  //write
  PatternFromStats::writeFile(outputFolder + "/patterns.txt", patterns);
}