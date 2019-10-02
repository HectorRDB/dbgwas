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

#ifndef KSGATB_PHENOCOUNTER_H
#define KSGATB_PHENOCOUNTER_H

#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <vector>
#include "Utils.h"

using namespace std;

//class used to count how many times a unitig is seen overall, in phenotype 0, 1 or NA
class PhenoCounter {
private:
    vector<double> validPhenotypes;
    int NACount;
public:
    PhenoCounter():validPhenotypes(), NACount(0){}

    //return # phenotypes <= threshold
    int getPheno0(double threshold) const {
      return std::count_if(validPhenotypes.begin(), validPhenotypes.end(), [&](double phenotype) {
          return phenotype <= threshold;
      });
    }

    //return # phenotypes > threshold
    int getPheno1(double threshold) const { return validPhenotypes.size() - getPheno0(threshold); }

    //return # phenotypes == NA
    int getNA() const { return NACount; }

    //get the total # of appearances of this unitig (total pheno counter = unitig counter)
    int getTotal() const { return validPhenotypes.size() + NACount; }

    //get the nb of valid phenotypes
    int getNbOfValidPhenos() const { return validPhenotypes.size(); }

    //add the given phenotype count times to this object
    void add(const string &phenotype, int count);


    //method to allow boost serialization
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & validPhenotypes;
      ar & NACount;
    }
};


#endif //KSGATB_PHENOCOUNTER_H
