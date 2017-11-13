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
  ss << dirWhereDBGWASIsInstalled << DBGWAS_lib << "/" << command << " -query " << queryPath << " -db " << dbPath << " -out " << outFilePath << " -num_threads " << nbCores << " -outfmt '6 qseqid sseqid qcovs bitscore pident evalue'";
  string commandLine=ss.str();
  executeCommand(commandLine, false);

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
    string commandLineMakeblastdb = dirWhereDBGWASIsInstalled + DBGWAS_lib + string("/makeblastdb -dbtype ") + dbtype + " -in " + fixedDBPath;
    executeCommand(commandLineMakeblastdb);
  }

  //return the fixed db path
  return fixedDBPath;
}


void AnnotationRecord::SetOfNodesAndEvalue::addNode(int node, long double evalue) {
  nodes.insert(node);
  minEvalue = min(minEvalue, evalue);
}

//get a representation of this annotation to be added to the SQL string in the index page
string AnnotationRecord::getSQLRepresentationForIndexPage() const {
  stringstream ss;
  if (annotations.size()==0)
    ss << UNIQUE_SYMBOL_MARKER << "No annotations found" << UNIQUE_SYMBOL_MARKER << " ";
  else {
    for (const auto & indexAndSetOfNodesAndEvalue : annotations)
      ss << UNIQUE_SYMBOL_MARKER << annotationIndex[indexAndSetOfNodesAndEvalue.first] << UNIQUE_SYMBOL_MARKER << " ";
  }
  return ss.str();
}

//get JS array representation of the annotation component for the index page
string AnnotationRecord::getAnnotationsForHOTForIndexPage(int componentId) const {
  stringstream ss;
  ss << "[";
  for (const auto & indexAndSetOfNodesAndEvalue : annotations)
    ss << "['" << annotationIndex[indexAndSetOfNodesAndEvalue.first] << "', " << indexAndSetOfNodesAndEvalue.second.getHTMLRepresentationForIndexPage() << "], ";
  ss << "]";

  return ss.str();
}


//transform to a javascript array
string AnnotationRecord::SetOfNodesAndEvalue::getHTMLRepresentationForGraphPage () const {
  stringstream ss;
  ss << scientific;
  ss << nodes.size() << ", " << minEvalue << ", [";
  for (const auto &node : nodes)
    ss << "'n" << node << "', ";
  ss << "]";
  return ss.str();
}

//transform to a javascript array
string AnnotationRecord::SetOfNodesAndEvalue::getHTMLRepresentationForIndexPage () const {
  stringstream ss;
  ss << scientific;
  ss << nodes.size() << ", " << minEvalue;
  return ss.str();
}

set<string> AnnotationRecord::getAnnotationIndexAsSet() const {
  set<string> allAnnotationsNames;
  for (const auto & annotation : annotationIndex)
    allAnnotationsNames.insert(annotation);
  return allAnnotationsNames;
}


//get a JS representation of the annotation component for the graph page, with the annotation index, nb of nodes and evalue
string AnnotationRecord::getJSRepresentationAnnotIdNbNodesEvalueForGraphPage() const {
  stringstream ss;
  ss << "[";
  for (const auto & indexAndSetOfNodesAndEvalue : annotations)
    ss << "[" << indexAndSetOfNodesAndEvalue.first << ", " << indexAndSetOfNodesAndEvalue.second.getHTMLRepresentationForGraphPage() << "], ";
  ss << "]";
  return ss.str();
}

//get all annotations names as a JS vector
string AnnotationRecord::getAnnotationIndexAsJSVector() const {
  stringstream ss;
  ss << "[";
  for (const auto & annotation : annotationIndex)
    ss << "'" << annotation << "', ";
  ss << "]";
  return ss.str();
}


//add an annotation to this set
void AnnotationRecord::addAnnotation(const string &tag, int node, long double evalue) {
  if (find(annotationIndex.begin(), annotationIndex.end(), tag) == annotationIndex.end())
    //did not find the tag, push it
    annotationIndex.push_back(tag);

  //find the index
  auto index = find(annotationIndex.begin(), annotationIndex.end(), tag)-annotationIndex.begin();

  //add to annotations
  annotations[index].addNode(node, evalue);
  nodeId2Annotations[node].addAnnotation(index, evalue);
}

void AnnotationRecord::AnnotationsAndEvalue::addAnnotation(int annotation, long double evalue) {
  if (annotation2Evalue.count(annotation)>0) //it is already present in the map
    annotation2Evalue[annotation] = min(annotation2Evalue[annotation], evalue);
  else
    annotation2Evalue[annotation] = evalue;
}


//get all the annotations IDs from a node as JS vector
string AnnotationRecord::getAllAnnotationsIDsFromANodeAsJSVector(int node) {
  stringstream ss;
  ss << "[";
  for (const auto &annotationAndEvalue : nodeId2Annotations[node].annotation2Evalue)
    ss << annotationAndEvalue.first << ", ";
  ss << "]";
  return ss.str();
}

//get a dictionary in JS where the key is the node id and the value is a pair annotation and evalue
string AnnotationRecord::getJSRepresentationNodeId2AnnotationsEvalueForGraphPage() const {
  stringstream ss;
  ss << scientific;

  ss << "{";

  for (const auto &nodeId2AnnotationsIt : nodeId2Annotations) {
    ss << "'n" << nodeId2AnnotationsIt.first << "': [";

    for (const auto &annotationEvalue : nodeId2AnnotationsIt.second.annotation2Evalue)
      ss << "['" << annotationEvalue.first << "', " << annotationEvalue.second << "]";

    ss << "]";
  }

  ss << "}";

  return ss.str();
};