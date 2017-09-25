//
// Created by Leandro Ishi Soares de Lima on 22/09/17.
//

#include "Blast.h"
#include "Utils.h"

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


vector<BlastRecord> Blast::blast (const string &command, const string &queryPath, const string &dbPath) {
  string outFilePath = queryPath+"."+command+"Out";

  //build the command line
  stringstream ss;
  ss << command << " -query " << queryPath << " -db " << dbPath << " -out " << outFilePath << " -outfmt '6 qseqid sseqid qcovs bitscore pident evalue'";
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


string Blast::makeblastdb (const string &dbtype, const string &originalDBPath) {
  string fixedDBPath = originalDBPath + ".DBGWAS.fasta";

  //replace spaces for underscores in the FASTA file, as this could create some problems...
  cout << "[WARNING] Copying and replacing spaces for _ in your DB " << originalDBPath << endl;
  string commandLineFixSpaces = string("tr ' ' '_' <") + originalDBPath + " >" + fixedDBPath;
  executeCommand(commandLineFixSpaces);

  //create the DB using the fixed FASTA
  string commandLineMakeblastdb = string("makeblastdb -dbtype ") + dbtype +  " -in " + fixedDBPath;
  executeCommand(commandLineMakeblastdb);

  //return the fixed db path
  return fixedDBPath;

}