//
// Created by Leandro Ishi Soares de Lima on 05/06/16.
//

#ifndef KISSPLICE_UTILS_H
#define KISSPLICE_UTILS_H

#include <gatb/gatb_core.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "pstream.h"
#include <ctime>
#include <sstream>
#include <iostream>
#include <limits>
#include <iterator>
#include <set>
#include <boost/algorithm/string/replace.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

using namespace std;
namespace fs = boost::filesystem;

char complement(char b);
string reverse_complement(const string &seq);
//Read all strings in the readsFile file and return them as a vector of strings
vector<string> getVectorStringFromFile(const string &readsFile);

//this function also populates strains if needed
void checkStrainsFile(const string &strainsFile);

string readFileAsString(const char* fileName);

void copyDirectoryRecursively(const fs::path& sourceDir, const fs::path& destinationDir);

int getNbLinesInFile(const string &filename);

void checkParametersBuildDBG(Tool *tool);
void checkParametersGenerateOutput(Tool *tool);
void checkParametersStatisticalTest(Tool *tool);
void fatalError (const string &message);
void executeCommand(const string &command, bool verbose=true);
void openFileForReading(const string &filePath, ifstream &stream);
void openFileForWriting(const string &filePath, ofstream &stream);
void createFolder(const string &path);



//global vars used by all programs
class UnitigIdStrandPos {
public:
    int unitigId;
    char strand; //TODO: change for bool
    int pos;
    int unitigSize;
    int kmerSize; //the size of the kmer (just for the reverseStrand() function) //TODO: static global var, sth like that

    UnitigIdStrandPos(int unitigId=0, char strand='?',int pos=0, int unitigSize=0, int kmerSize=0):
        unitigId(unitigId), strand(strand), pos(pos), unitigSize(unitigSize), kmerSize(kmerSize){}
    int getUnitigId () const {checkValidity(); return unitigId; }
    char getStrand () const {checkValidity(); return strand; }
    int getPos() const {checkValidity(); return pos; }
    string toString() const {
      stringstream ss;
      try {
        ss << getUnitigId() << " " << getStrand() << " " << getPos();
      }catch (const runtime_error &e) {
        ss << "-1 ? -1";
      }
      return ss.str();
    }
    void reverseStrand () {
      strand = (strand=='F' ? 'R' : 'F');
      pos = unitigSize-pos-kmerSize;
    }

    //check if the unitig is valid or not
    void checkValidity() const {
      if (strand=='?')
        throw runtime_error("Invalid unitig.");
    }
};

struct Strain {
    string id, phenotype, path;
    Strain(const string &id, const string &phenotype, const string &path) : id(id), phenotype(phenotype) {
      //transfor to canonical path
      boost::filesystem::path boostPath(boost::filesystem::canonical(path));
      this->path = boostPath.string();
    }

    static void createReadsFile(const string &readsFile, vector< Strain >* strains) {
      ofstream fout;
      openFileForWriting(readsFile, fout);

      for (const auto &strain : (*strains))
        fout << strain.path << endl;

      fout.close();
    }

    static void createIdPhenoFile(const string &filePath, vector< Strain >* strains) {
      ofstream fout;
      openFileForWriting(filePath, fout);
      fout << "ID\tpheno" << endl;

      for (const auto &strain : (*strains))
        fout << strain.id << "\t" << strain.phenotype << endl;

      fout.close();
    }

    //create a file with 2 ints. The first is the nb of pheno0 strains and the second is the nb of pheno1 strains
    static void createFileWithAmountOfStrainsInEachPheno(const string &filePath, vector< Strain >* strains) {
      ofstream fout;
      openFileForWriting(filePath, fout);

      int pheno0Count=0;
      int pheno1Count=0;
      for (const auto &strain : (*strains)) {
        if (strain.phenotype == "0")
          pheno0Count++;
        if (strain.phenotype == "1")
          pheno1Count++;
      }

      fout << pheno0Count << " " << pheno1Count;

      fout.close();
    }
};




struct PatternFromStats {
    int pattern;
    long double qValue;
    long double weight;
    long double normalizedWeight;
    long double waldStatistic;

    bool operator < (const PatternFromStats& other) const {
      return this->qValue < other.qValue;
    }

    static vector<PatternFromStats> readFile(const string &filename, bool header=false) {
      vector<PatternFromStats> patterns;

      //read
      {
        ifstream patternStream;
        openFileForReading(filename, patternStream);
        patternStream >> setprecision(std::numeric_limits<long double>::digits10 + 1);
        PatternFromStats pattern;

        //remove header
        if (header) {
          string tmp;
          getline(patternStream, tmp);
        }
        while (patternStream >> pattern.pattern >> pattern.qValue >> pattern.weight >> pattern.waldStatistic)
          patterns.push_back(pattern);
        patternStream.close();
      }

      //normalize
      vector<long double>  normalizedWeights;
      for (const auto &pattern : patterns)
        normalizedWeights.push_back(pattern.weight);

      //remove the minimum weight of everyone
      long double minWeight = *min_element(normalizedWeights.begin(), normalizedWeights.end());
      transform(normalizedWeights.begin(), normalizedWeights.end(), normalizedWeights.begin(), [&](long double weight) {
          return weight-minWeight;
      });

      //transform to [0,1]
      long double maxWeight = *max_element(normalizedWeights.begin(), normalizedWeights.end());
      transform(normalizedWeights.begin(), normalizedWeights.end(), normalizedWeights.begin(), [&](long double weight) {
          return weight/maxWeight;
      });

      //write back
      auto normalizedWeightsIt = normalizedWeights.begin();
      for (auto &pattern : patterns) {
        pattern.normalizedWeight = *normalizedWeightsIt;
        ++normalizedWeightsIt;
      }

      return patterns;
    }

    static void writeFile(const string &filename, const vector<PatternFromStats> &patterns) {
      ofstream patternSortedStream;
      openFileForWriting(filename, patternSortedStream);
      patternSortedStream << setprecision(std::numeric_limits<long double>::digits10 + 1);
      patternSortedStream << "pattern q-value weight wald_statistic" << endl;
      for (const auto &pattern : patterns)
        patternSortedStream << pattern.pattern << " " << pattern.qValue << " " << pattern.weight << " " << pattern.waldStatistic << endl;
      patternSortedStream.close();
    }
};





//get the significant patterns according to SFF
class GetSignificantPatterns
    : public boost::static_visitor<>
{
private:
    const vector<PatternFromStats> &patterns;
    vector<PatternFromStats> &significantPatterns;
public:
    GetSignificantPatterns(const vector<PatternFromStats> &patterns, vector<PatternFromStats> &significantPatterns) :
        patterns(patterns), significantPatterns(significantPatterns){}

    void operator()(int &n) const;

    void operator()(double &qValue) const;

};
#endif //KISSPLICE_UTILS_H
