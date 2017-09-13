//
// Created by Leandro Ishi Soares de Lima on 22/09/17.
//

#ifndef DBGWAS_BLAST_H
#define DBGWAS_BLAST_H

#include <sstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;


struct ValueNotFound{};

class BlastRecord {
public:
    int nodeId; //the query id
    string DBGWAS_index_tag, DBGWAS_graph_tag; //the tags, if available
    double qcovs, bitscore, pident, evalue; //the blast fields we care about

    //constructors
    BlastRecord(){}
    BlastRecord(int nodeId, const string &DBGWAS_index_tag, const string &DBGWAS_graph_tag,
                double qcovs, double bitscore, double pident, double evalue):
        nodeId(nodeId), DBGWAS_index_tag(DBGWAS_index_tag), DBGWAS_graph_tag(DBGWAS_graph_tag), qcovs(qcovs),
        bitscore(bitscore), pident(pident), evalue(evalue) {}

    //parse a string and build a BlastRecord from it
    static BlastRecord parseString (const string &str);
private:
    //parse the header and extract the value corresponding to the given tag
    static string extractValue (const string &header, const string &tag);
};

class Blast {
public:
    //Blast (using command) the file in queryPath agains the db on dbPath and return the results as a vector<BlastRecord>
    static vector<BlastRecord> blast (const string &command, const string &queryPath, const string &dbPath);
};



#endif //DBGWAS_BLAST_H
