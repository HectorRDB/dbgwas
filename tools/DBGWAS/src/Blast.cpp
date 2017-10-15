//
// Created by Leandro Ishi Soares de Lima on 22/09/17.
//

#include "Blast.h"
#include "Utils.h"
#include "global.h"
#include <algorithm>

//parse a string and build a BlastRecord from it
BlastRecord BlastRecord::parseString (const string &str) {
  BlastRecord record;
  stringstream stream;
  stream << str;

  //read
  string header;
  stream >> record.nodeId >> header >> record.qcovs >> record.bitscore >> record.pident >> record.evalue;

  //escape '
  boost::replace_all(header, "'", "\\'");

  //parse DBGWAS_index_tag
  try {
    record.DBGWAS_index_tag = extractValue(header, "DBGWAS_index_tag");
  } catch (const ValueNotFound &e) {
    record.DBGWAS_index_tag = header;
  }
  if (record.DBGWAS_index_tag.size()==0) {
    cerr << "[WARNING] DBGWAS_index_tag of " << header << " is empty! Setting to EMPTY" << endl;
    record.DBGWAS_index_tag = "EMPTY";
  }


  //parse DBGWAS_graph_tag
  try {
    record.DBGWAS_graph_tag = extractValue(header, "DBGWAS_graph_tag");
  } catch (const ValueNotFound &e) {
    record.DBGWAS_graph_tag = header;
  }
  if (record.DBGWAS_graph_tag.size()==0) {
    cerr << "[WARNING] DBGWAS_graph_tag of " << header << " is empty! Setting to EMPTY" << endl;
    record.DBGWAS_graph_tag = "EMPTY";
  }

  return record;
}

//parse the header and extract the value corresponding to the given tag
string BlastRecord::extractValue (const string &header, const string &tag) {
  if (header.find(tag) != string::npos) {
    //found the tag, extract the value
    auto posDBGWASIndexTag = header.find(tag);
    auto posSemicolon = header.find_first_of("; \t\n\r", posDBGWASIndexTag+1);
    if (posSemicolon==string::npos) posSemicolon = header.length();
    auto posStartValue = posDBGWASIndexTag+string(tag).length()+1;
    string value = header.substr(posStartValue, posSemicolon-posStartValue);
    boost::trim(value);
    return value;
  }

  //tag not found
  throw ValueNotFound();
}


vector<BlastRecord> Blast::blast (const string &command, const string &queryPath, const string &dbPath, int nbCores) {
  string outFilePath = queryPath+"."+command+"Out";

  //build the command line
  stringstream ss;
  ss << "./" << command << " -query " << queryPath << " -db " << dbPath << " -out " << outFilePath << " -num_threads " << nbCores << " -outfmt '6 qseqid sseqid qcovs bitscore pident evalue'";
  string commandLine=ss.str();
  executeCommand(commandLine);

  //read output
  vector<BlastRecord> records;
  {
    auto recordsAsStringVector = getVectorStringFromFile(outFilePath);
    for (const auto& recordString : recordsAsStringVector) {
      records.push_back(BlastRecord::parseString(recordString));
    }
  }

  return records;
}


string Blast::makeblastdb (const string &dbtype, const string &originalDBPath, const string &outputFolderPath) {
  string concatenatedDBPath;
  {
    stringstream ss;
    ss << outputFolderPath << "/" << dbtype << "_db";
    concatenatedDBPath = ss.str();
  }

  //concatenate all DBs into one
  {
    //TODO: I don't know why <(echo) is not allowed in the cat command, so we have to use this ugly thing...
    //create an empty file with only a newline
    string newlineFilepath = outputFolderPath + "/newline";
    executeCommand(string("echo > ") + newlineFilepath);
    //TODO: I don't know why <(echo) is not allowed in the cat command, so we have to use this ugly thing...

    string catCommand;
    stringstream ss;
    ss << "cat " << originalDBPath << " > " << concatenatedDBPath;
    catCommand = ss.str();
    boost::replace_all(catCommand, ",", string(" ") + newlineFilepath + " ");
    executeCommand(catCommand);
  }

  //replace spaces for underscores in the concatenated files, as this could create some problems...
  string fixedDBPath(concatenatedDBPath);
  fixedDBPath += "_fixed";
  {
    string commandLineFixSpaces = string("tr ' ' '_' <") + concatenatedDBPath + " >" + fixedDBPath;
    executeCommand(commandLineFixSpaces);
  }

  //create the DB using the fixed FASTA
  {
    string commandLineMakeblastdb = string("./makeblastdb -dbtype ") + dbtype + " -in " + fixedDBPath;
    executeCommand(commandLineMakeblastdb);
  }

  //return the fixed db path
  return fixedDBPath;
}


void AnnotationRecord::SetOfNodesAndEvalue::addNode(int node, long double evalue) {
  nodes.insert(node);
  evalue=min(minEvalue, evalue);
}


string AnnotationRecord::getSQLRepresentation() const {
  stringstream ss;
  if (annotations.size()==0)
    ss << UNIQUE_SYMBOL_MARKER << "No annotations found" << UNIQUE_SYMBOL_MARKER << " ";
  else {
    for (const auto & tagAndSetOfNodesAndEvalue : annotations)
      ss << UNIQUE_SYMBOL_MARKER << tagAndSetOfNodesAndEvalue.first << UNIQUE_SYMBOL_MARKER << " ";
  }
  return ss.str();
}

string AnnotationRecord::getHTMLRepresentationForIndexPage(int componentId) const {
  stringstream ss;
  if (annotations.size()==0)
    ss << "<b>No annotations found.</b>";
  else {
    ss << "<script> buildHandsonTableForAnnotationIndexPage(" << componentId << ", [";
    for (const auto & tagAndSetOfNodesAndEvalue : annotations)
      ss << "[\\\'" << tagAndSetOfNodesAndEvalue.first << "\\\', " << tagAndSetOfNodesAndEvalue.second.getHTMLRepresentationForIndexPage() << "], ";
    ss << "]) </script>";
  }
  return ss.str();
}


//transform to a javascript array
string AnnotationRecord::SetOfNodesAndEvalue::getHTMLRepresentationForIndexPage () const {
  stringstream ss;
  ss << nodes.size() << ", " << minEvalue << ", [";
  for (const auto &node : nodes)
    ss << "\\\'n" << node << "\\\', ";
  ss << "]";
  return ss.str();
}

set<string> AnnotationRecord::getAllAnnotationsNames() const {
  set<string> allAnnotationsNames;
  for (const auto & tagAndSetOfNodesAndEvalue : annotations)
    allAnnotationsNames.insert(tagAndSetOfNodesAndEvalue.first);
  return allAnnotationsNames;
}


//get an HTML representation of the annotation component for the graph page
string AnnotationRecord::getHTMLRepresentationForGraphPage() const {
  stringstream ss;
  ss << "[";
  for (const auto & tagAndSetOfNodesAndEvalue : annotations)
    ss << "['" << tagAndSetOfNodesAndEvalue.first << "', " << tagAndSetOfNodesAndEvalue.second.getHTMLRepresentationForIndexPage() << "], ";
  ss << "]";
  return ss.str();
}


//add an annotation to this set
void AnnotationRecord::addAnnotation(const string &tag, int node, long double evalue) {
  annotations[tag].addNode(node, evalue);
  nodeId2Annotation[node].insert(tag);
}