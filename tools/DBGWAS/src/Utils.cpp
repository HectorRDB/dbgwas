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

#include "Utils.h"
#include "global.h"
#include "Blast.h"

using namespace std;

char complement(char b)
{
  switch(b)
  {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';

    case 'a': return 't';
    case 't': return 'a';
    case 'g': return 'c';
    case 'c': return 'g';

    case 'N': return 'N';
    case '*': return '*';
  }
  return '?';
}

string reverse_complement(const string &seq)
{
  string s(seq.begin(),seq.end());
  string::iterator pos;

  reverse(s.begin(), s.end());

  for(pos=s.begin();pos!=s.end();++pos)
    *pos=complement(*pos);

  return s;
}


//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile) {
  vector<string> allReadFilesNames;
  string tempStr;

  ifstream readsFileStream;
  openFileForReading(readsFile, readsFileStream);
  while (getline(readsFileStream, tempStr)) {
    if (tempStr.size() > 0)
      allReadFilesNames.push_back(tempStr);
  }
  readsFileStream.close();

  return allReadFilesNames;
}

//this function also populates strains if needed
void checkStrainsFile(const string &strainsFile) {
  vector<Strain> localStrains;
  bool header=true;
  ifstream input;
  openFileForReading(strainsFile, input);
  set<string> allIds;

  for(string line; getline( input, line ); )
  {
    //parse header
    if (header) {
      header=false;
      continue;
    }

    //ignore empty lines
    if (line.size()==0)
      continue;

    //create the strain
    stringstream ss;
    ss << line;
    string id, pheno, path;
    ss >> id >> pheno >> path;

    //check for duplicated IDs
    if (allIds.find(id)!=allIds.end()) {
      stringstream ss;
      ss << "Duplicated IDs in " << strainsFile << ": " << id << endl;
      fatalError(ss.str());
    }
    allIds.insert(id);

    //check for disallowed phenotypes
    bool phenoIsNumber;
    double phenoAsNumber;
    std::tie(phenoIsNumber, phenoAsNumber) = is_number(pheno);

    //allowed phenotypes are only "NA" or numbers between 0 and 1
    if (!phenoIsNumber && pheno!="NA") {
      stringstream ss;
      ss << "Phenotype not allowed: " << pheno << " . The only allowed values for phenotypes are real numbers or NA." << endl;
      fatalError(ss.str());
    }

    //check if the path is ok
    ifstream file;
    openFileForReading(path, file);
    if (!file.is_open()) {
      stringstream ss;
      ss << "Error opening file " << path << " in " << strainsFile << endl;
      fatalError(ss.str());
    }
    file.close();


    //add the strain if it is different from NA
    if (pheno=="NA" && keepNA==false) {
      cerr << "[WARNING] Skipping strain " << id << " because its phenotype is NA and " << STR_KEEP_NA << " is not set." << endl;
    }else {
      Strain strain(id, pheno, path);
      localStrains.push_back(strain);
    }
  }
  input.close();

  //in the end, check if strain is null. If it is, populate it
  if (strains==NULL)
    strains = new vector<Strain>(localStrains);
}



string readFileAsString(const char* fileName) {
  std::ifstream t;
  openFileForReading(fileName, t);
  std::string str;

  t.seekg(0, std::ios::end);
  str.reserve(t.tellg());
  t.seekg(0, std::ios::beg);

  str.assign((std::istreambuf_iterator<char>(t)),
             std::istreambuf_iterator<char>());
  t.close();
  return str;
}


void copyDirectoryRecursively(const fs::path& sourceDir, const fs::path& destinationDir)
{
  if (!fs::exists(sourceDir) || !fs::is_directory(sourceDir))
  {
    throw std::runtime_error("Source directory " + sourceDir.string() + " does not exist or is not a directory");
  }
  if (fs::exists(destinationDir))
  {
    throw std::runtime_error("Destination directory " + destinationDir.string() + " already exists");
  }
  if (!fs::create_directory(destinationDir))
  {
    throw std::runtime_error("Cannot create destination directory " + destinationDir.string());
  }

  for (const auto& dirEnt : fs::recursive_directory_iterator{sourceDir})
  {
    const auto& path = dirEnt.path();
    auto relativePathStr = path.string();
    boost::replace_first(relativePathStr, sourceDir.string(), "");
    fs::copy(path, destinationDir / relativePathStr);
  }
}


int getNbLinesInFile(const string &filename) {
  std::ifstream file;
  openFileForReading(filename.c_str(), file);

  // Number of lines in the file
  int n = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');

  file.close();

  return n;
}

void checkExecutables(Tool *tool) {
  //configure the global vars of the executable paths
  gemmaPath = tool->getInput()->getStr(STR_GEMMA_PATH);
  blastPath = tool->getInput()->getStr(STR_BLAST_PATH);
  phantomjsPath = tool->getInput()->getStr(STR_PHANTOMJS_PATH);
  RscriptPath = tool->getInput()->getStr(STR_RSCRIPT_PATH);

  //fix the variables
  boost::replace_all(gemmaPath, "<DBGWAS_lib>", DBGWAS_lib);
  boost::replace_all(blastPath, "<DBGWAS_lib>", DBGWAS_lib);
  boost::replace_all(phantomjsPath, "<DBGWAS_lib>", DBGWAS_lib);
  boost::replace_all(RscriptPath, "<DBGWAS_lib>", DBGWAS_lib);

  //check if the executables work
  //check gemma
  executeCommand(gemmaPath, false, gemmaPath + " does not work, but it is required. You can install a version of GEMMA that works on your system and tell DBGWAS to use it through the parameter " + STR_GEMMA_PATH + " ."); //if it returns an exit status != 0, then it does not work and we issue a fatal error

  //check Rscript
  executeCommand(RscriptPath+" --version", false, RscriptPath + " does not work, but it is required. You can install a version of R that works on your system and tell DBGWAS to use it through the parameter " + STR_RSCRIPT_PATH + " .");

  //check blast, if the user wants to annotate the subgraphs
  if (tool->getInput()->get(STR_NUCLEOTIDE_DB)) {
    executeCommand(blastPath+"/makeblastdb -version", false, string("You want to annotate the output subgraphs with ") + STR_NUCLEOTIDE_DB + " , but " + blastPath+"/makeblastdb does not work. You can install a version of the Blast suite that works on your system" +
        " and tell DBGWAS to use it through the parameter " + STR_BLAST_PATH + " .");
    executeCommand(blastPath+"/blastn -version", false, string("You want to annotate the output subgraphs with ") + STR_NUCLEOTIDE_DB + " , but " + blastPath+"/blastn does not work. You can install a version of the Blast suite that works on your system" +
        " and tell DBGWAS to use it through the parameter " + STR_BLAST_PATH + " .");
  }
  if (tool->getInput()->get(STR_PROTEIN_DB)) {
    executeCommand(blastPath+"/makeblastdb -version", false, string("You want to annotate the output subgraphs with ") + STR_PROTEIN_DB + " , but " + blastPath+"/makeblastdb does not work. You can install a version of the Blast suite that works on your system" +
        " and tell DBGWAS to use it through the parameter " + STR_BLAST_PATH + " .");
    executeCommand(blastPath+"/blastx -version", false, string("You want to annotate the output subgraphs with ") + STR_PROTEIN_DB + " , but " + blastPath+"/blastx does not work. You can install a version of the Blast suite that works on your system" +
        " and tell DBGWAS to use it through the parameter " + STR_BLAST_PATH + " .");
  }

  //check phantomjs
  if (tool->getInput()->get(STR_NO_PREVIEW) == 0) {
    executeCommand(phantomjsPath+" --version", false, phantomjsPath + " does not work, DBGWAS cant produce the components preview on the summary output page." +
                                             " Either you get a version of phantomjs that works on your system (DBGWAS was "
                                             "tested with version 2.1.1) and tell DBGWAS to use it through the parameter " + STR_PHANTOMJS_PATH +
                                             " , or you can choose to not produce the components preview "
                                             "through the parameter " + STR_NO_PREVIEW + " .");
  }
}


void checkParametersBuildDBG(Tool *tool) {
  //check if we skip or not
  skip1 = tool->getInput()->get(STR_SKIP1) != 0;
  skip2 = tool->getInput()->get(STR_SKIP2) != 0;
  keepNA = tool->getInput()->get(STR_KEEP_NA) != 0;
  hasNewickFile = tool->getInput()->get(STR_NEWICK_PATH) != 0;

  if (skip2) skip1=true;

  if (skip1) {
    cerr << "Skipping Step 1!" << endl;
    return;
  }

  //check executables and parameters
  checkExecutables(tool);

  //check the count mode
  //TODO: seeveral questions are still unclear if we use the Freq count mode (how to run bugwas, the coloring, etc...). For now I am disabling this option
  /*
  string countMode = tool->getInput()->getStr(STR_COUNT_MODE);
  if (countMode!="01" && countMode!="Freq") {
    stringstream ss;
    ss << "Wrong value for parameter " << STR_COUNT_MODE << ". Value found: " << countMode << " . Values accepted: 01 or Freq.";
    fatalError(ss.str());
  }
   */

  //check the strains file
  string strainsFile = tool->getInput()->getStr(STR_STRAINS_FILE);
  checkStrainsFile(strainsFile);

  //check output
  string outputFolderPath = stripLastSlashIfExists(tool->getInput()->getStr(STR_OUTPUT));
  boost::filesystem::path p(outputFolderPath.c_str());
  if (boost::filesystem::exists(p)) {
    stringstream ss;
    ss << "Could not create dir " << outputFolderPath << " - path already exists. Remove it and re-run the tool, or use -skip1 or -skip2 parameters.";
    fatalError(ss.str());
  }
  createFolder(p.string());
}

//parse SFF
void parseSFF(const string &SFFString) {
  qOrPValue = SFFString[0];
  string SFFNumber = SFFString.substr(1);
  if (qOrPValue!='q' && qOrPValue!='p')
    fatalError(string("Error on ") + string(STR_SFF) + " parameter. First argument must be p or q.");

  if (SFFNumber.find(".")==string::npos) {
    //. not found in SFFNumber : integer
    //get the first n significant patterns
    int n;
    {
      stringstream ss;
      ss << SFFNumber;
      ss >> n;
      if (ss.fail())
        fatalError(string("Error on ") + string(STR_SFF) + " parameter. Second argument must be an integer or a double.");
    }

    SFF=n;
  }else {
    //double
    //get all sequence in which the p/q-value is <= n
    double n;
    {
      stringstream ss;
      ss << SFFNumber;
      ss >> n;
      if (ss.fail())
        fatalError(string("Error on ") + string(STR_SFF) + " parameter. Second argument must be an integer or a double.");
    }
    SFF = n;
  }
}


void checkParametersStatisticalTest(Tool *tool) {
  if (skip2) {
    cerr << "Skipping Step 2!" << endl;
    return;
  }

  //check executables and parameters
  checkExecutables(tool);

  //check if newickTreeFilePath exists
  if (hasNewickFile) {
    string newickTreeFilePath = tool->getInput()->getStr(STR_NEWICK_PATH);
    boost::filesystem::path p(newickTreeFilePath.c_str());
    if (!boost::filesystem::exists(p)) {
      stringstream ss;
      ss << "Error locating newick tree file path: " << newickTreeFilePath;
      fatalError(ss.str());
    }
  }

  //parse SFF
  string SFFString = tool->getInput()->getStr(STR_SFF);
  parseSFF(SFFString);
}


void checkParametersGenerateOutput(Tool *tool) {
  //check executables and parameters
  checkExecutables(tool);

  //create the output folder for step 3
  string outputFolder = stripLastSlashIfExists(tool->getInput()->getStr(STR_OUTPUT))+string("/step3");
  createFolder(outputFolder);

  //create the tmp folder of step 3
  string tmpFolder = outputFolder+string("/tmp");
  createFolder(tmpFolder);


  //remove old visualisation folder and create new
  string visualisationFolder = stripLastSlashIfExists(tool->getInput()->getStr(STR_OUTPUT))+string("/visualisations");
  removeOldAndCreateFolder(visualisationFolder, "maybe previous visualisations?");
  string componentsFolder = visualisationFolder + string("/components");
  createFolder(componentsFolder);

  //do the same for the textual output
  string textualOutputFolder = stripLastSlashIfExists(tool->getInput()->getStr(STR_OUTPUT))+string("/textualOutput");
  removeOldAndCreateFolder(textualOutputFolder, "maybe previous textual output?");
  string textualComponentsFolder = textualOutputFolder + string("/components");
  createFolder(textualComponentsFolder);


  //check the nucleotide DB
  if (tool->getInput()->get(STR_NUCLEOTIDE_DB)) {
    //build the nucleotide DB
    nucleotideDBPath = Blast::makeblastdb("nucl", tool->getInput()->getStr(STR_NUCLEOTIDE_DB), outputFolder);
    thereIsNucleotideDB=true;
  }

  //check the protein DB
  if (tool->getInput()->get(STR_PROTEIN_DB)) {
    //build the protein DB
    proteinDBPath = Blast::makeblastdb("prot", tool->getInput()->getStr(STR_PROTEIN_DB), outputFolder);
    thereIsProteinDB=true;
  }

  //get the -no-preview parameter
  noPreview = tool->getInput()->get(STR_NO_PREVIEW) != 0;

  //get the phenotype threshold
  phenotypeThreshold = tool->getInput()->getDouble(STR_PHENOTYPE_THRESHOLD);

  //parse SFF
  string SFFString = tool->getInput()->getStr(STR_SFF);
  parseSFF(SFFString);
}


void fatalError (const string &message) {
  cerr << endl << endl << "[FATAL ERROR] " << message << endl << endl;
  cerr.flush();
  exit(1);
}


void executeCommand(const string &command, bool verbose, const string &messageIfItFails) {
  // run a process and create a streambuf that reads its stdout and stderr
  if (verbose)
    cerr << "Executing " << command << "..." << endl;

  //create the process
  redi::ipstream proc(command, redi::pstreams::pstdout | redi::pstreams::pstderr);
  string line;

  // read child's stdout
  while (getline(proc.out(), line)) {
    if (verbose)
      cout << line << endl;
  }
  // read child's stderr
  while (getline(proc.err(), line)) {
    if (verbose)
      cerr << line << endl;
  }

  //check exit status
  proc.close();
  if (proc.rdbuf()->exited()) {
    if (proc.rdbuf()->status() != 0) {
      stringstream ss;
      ss << "Error executing " << command << ". Exit status: " << proc.rdbuf()->status() << endl;
      if (messageIfItFails != "")
        ss << "Message: " << messageIfItFails << endl;
      fatalError(ss.str());
    }
    if (verbose)
      cerr << "Executing " << command << " - Done!" << endl;
  }
  else
    fatalError("On executeCommand()");
}

//strips all last "/" if exists in the parameter
string stripLastSlashIfExists (string path) {
  while(path.size()>0 && path.back()=='/')
    path.pop_back();
  return path;
}


void openFileForReading(const string &filePath, ifstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

void openFileForWriting(const string &filePath, ofstream &stream) {
  stream.open(filePath);
  if (!stream.is_open()) {
    stringstream ss;
    ss << "Error opening file " << filePath;
    fatalError(ss.str());
  }
}

void createFolder(const string &path) {
  boost::filesystem::path folder(path.c_str());

  if (boost::filesystem::exists(folder))
    return;

  if (!boost::filesystem::create_directories(folder)) {
    stringstream ss;
    ss << "Could not create dir " << path << " - unknown reasons...";
    fatalError(ss.str());
  }
}

void removeOldAndCreateFolder(const string &path, const string &reason){
  boost::filesystem::path folder(path.c_str());
  if (boost::filesystem::exists(folder)) {
    cerr << "[WARNING] Removing " << path << " because path already exists (" << reason << "). " << endl;
    boost::filesystem::remove_all(folder);
  }
  createFolder(folder.string());
}


void GetSignificantPatterns::operator()(int &n) const
{
  n = min(n, (int)patterns.size());
  int i=0;
  for (const auto &pattern : patterns) {
    if (i < n)
      significantPatterns.push_back(pattern);
    i++;
  }
}

void GetSignificantPatterns::operator()(double &qOrPValueThreshold) const
{
  for (const auto &pattern : patterns) {
    if ((qOrPValue=='p' && pattern.pValue <= qOrPValueThreshold) || (qOrPValue=='q' && pattern.qValue <= qOrPValueThreshold))
      significantPatterns.push_back(pattern);
  }
}


string getDirWhereDBGWASIsInstalled() {
  char* path = NULL;
  int length, dirname_length;
  string toReturn;

  length = wai_getExecutablePath(NULL, 0, &dirname_length);
  if (length > 0)
  {
    path = (char*)malloc(length + 1);
    if (!path)
      //error, malloc did not work
      fatalError("Error on Utils.cpp::getDirWhereDBGWASIsInstalled()");

    //get the executable path
    wai_getExecutablePath(path, length, &dirname_length);

    //get the executable dir
    path[length] = '\0';
    path[dirname_length] = '\0';

    //save it in toReturn
    toReturn = string(path);

    //free the memory
    free(path);
  }
  else {
    fatalError("Error on Utils.cpp::getDirWhereDBGWASIsInstalled()");
  }

  return toReturn;
}


//tries to parse s, and returns a pair<bool, double>
//the first value indicates if s was successfully parsed into a double
//the second value indicates the double (it is only valid if the first is true)
tuple<bool, double> is_number(const std::string& s) {
  double number;
  try {
    number = std::stod(s);
  }
  catch(...) {
    return make_tuple(false, number);
  }
  return make_tuple(true, number);
}


void Strain::createPhenotypeCounter(const string &filePath, vector< Strain >* strains) {
  PhenoCounter phenoCounter;
  for (const auto &strain : (*strains))
    phenoCounter.add(strain.phenotype, 1);
  //serialize phenoCounter
  ofstream phenoCounterFile;
  {
    openFileForWriting(filePath, phenoCounterFile);
    boost::archive::text_oarchive boostOutputArchive(phenoCounterFile);
    //serialization itself
    boostOutputArchive & phenoCounter;
  } //boostOutputArchive is closed on destruction
  phenoCounterFile.close();
}