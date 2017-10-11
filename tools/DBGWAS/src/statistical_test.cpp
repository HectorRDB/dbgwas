/*
## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
## Centre National de la Recherche Scientifique>

## 1. This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero General Public License as published
## by the Free Software Foundation version 3 of the  License and under the
## terms of article 2 below.
## 2. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General
## Public License for more details.
## You should have received a copy of the GNU Affero General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 3. Communication to the public by any means, in particular in the form of
## a scientific paper, a poster, a slideshow, an internet page, or a patent,
## of a result obtained directly or indirectly by running this program must
## cite the following paper :
##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix,
##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic
##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
## -------------------------------------------------------------------------

## Authors (alphabetically): Jacob L., Jaillard M., Lima L.
*/

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
  double mafFilter = getInput()->getDouble(STR_MAF_FILTER);

  //execute the statistical test
  //To do so, we need to create an output folder inside the outputFolder
  createFolder(string("output/output"));
  //We also need to create this other folder... bizzare...
  createFolder(string("output/bugwas_out__PCloadings/output"));

  //create the command line
  stringstream ssCommand;
  ssCommand << "Rscript --vanilla --no-save --no-restore DBGWAS.R "
            << outputFolder << " "
            << outputFolder << "/bugwas_input.id_phenotype "
            << outputFolder << "/bugwas_out_ "
            << pathToExecParent << "gemma.0.93b "
            << mafFilter << " ";
  if (hasNewickFile)
    ssCommand << newickTreeFilePath << " ";
  ssCommand << "2>&1";

  //execute the command line
  executeCommand(ssCommand.str());

  //sort the file by q-value and output it to output/patterns.txt
  //read
  auto patterns = PatternFromStats::readFile(outputFolder + "/bugwas_out__DBGWAS_patterns.txt");
  //sort
  sort(patterns.begin(), patterns.end());
  //write
  PatternFromStats::writeFile(outputFolder + "/patterns.txt", patterns);

  //print some stats
  cout << "################################################################################" << endl;
  cout << "Stats: " << endl;
  cout << "Total number of patterns: " << patterns.size() << endl;
  cout << "################################################################################" << endl;
}