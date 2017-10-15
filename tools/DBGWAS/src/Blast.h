//
// Created by Leandro Ishi Soares de Lima on 22/09/17.
//

#ifndef DBGWAS_BLAST_H
#define DBGWAS_BLAST_H

#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>

using namespace std;


struct ValueNotFound{};

class BlastRecord {
public:
    int nodeId; //the query id
    string DBGWAS_index_tag, DBGWAS_graph_tag; //the tags, if available
    double qcovs, bitscore, pident; //the blast fields we care about
    long double evalue; //the blast fields we care about

    //constructors
    BlastRecord(){}
    BlastRecord(int nodeId, const string &DBGWAS_index_tag, const string &DBGWAS_graph_tag,
                double qcovs, double bitscore, double pident, long double evalue):
        nodeId(nodeId), DBGWAS_index_tag(DBGWAS_index_tag), DBGWAS_graph_tag(DBGWAS_graph_tag), qcovs(qcovs),
        bitscore(bitscore), pident(pident), evalue(evalue) {}

    //parse a string and build a BlastRecord from it
    static BlastRecord parseString (const string &str);
private:
    //parse the header and extract the value corresponding to the given tag
    static string extractValue (const string &header, const string &tag);
};


//For each set of annotations of a component, record the annotations' name, the nodes mapping to each and the lowest e-value of a hit to it
class AnnotationRecord {
private:
    //stores all the nodes that map to an annotation and their minimum e-value
    class SetOfNodesAndEvalue {
    private:
        set<int> nodes;
        long double minEvalue;
    public:
        SetOfNodesAndEvalue():nodes(), minEvalue(std::numeric_limits<long double>::max()){}

        //add a new node to this set
        void addNode(int node, long double evalue);

        //get the nb of nodes mapping here
        int getNbOfNodes () const { return nodes.size(); }

        //get the evalue
        long double getMinEvalue () const { return minEvalue; }

        //transform to a javascript array
        string getHTMLRepresentationForGraphPage () const;

        //transform to a javascript array
        string getHTMLRepresentationForIndexPage () const;
    };
    map<string, SetOfNodesAndEvalue> annotations;
    map<int, set<string> > nodeId2Annotation;
public:
    AnnotationRecord():annotations(){}

    //add an annotation to this set
    void addAnnotation(const string &tag, int node, long double evalue);

    //get a representation of this annotation to be added to the SQL string in the index page
    string getSQLRepresentation() const;

    //get an HTML representation of the annotation component for the index page
    string getHTMLRepresentationForIndexPage(int componentId) const;

    //get an HTML representation of the annotation component for the graph page
    string getHTMLRepresentationForGraphPage() const;

    //get all the annotations names
    set<string> getAllAnnotationsNames() const;

    //get all the annotations names from a node
    set<string> getAllAnnotationsNamesFromANode(int node) {
        return nodeId2Annotation[node];
    }

};



class Blast {
public:
    //Blast (using command) the file in queryPath agains the db on dbPath and return the results as a vector<BlastRecord>
    static vector<BlastRecord> blast (const string &command, const string &queryPath, const string &dbPath, int nbCores);

    //make the blast DB
    //returns the path to the blast DB
    static string makeblastdb (const string &dbtype, const string &originalDBPath, const string &outputFolderPath);
};



#endif //DBGWAS_BLAST_H
