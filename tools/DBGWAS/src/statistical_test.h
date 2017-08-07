//
// Created by Leandro Ishi Soares de Lima on 10/07/17.
//

#ifndef CDBGGWAS_STATISTICAL_TEST_H
#define CDBGGWAS_STATISTICAL_TEST_H

#include <gatb/gatb_core.hpp>
using namespace std;

class statistical_test : public Tool {
public:

    // Constructor
    statistical_test();

    // Actual job done by the tool is here
    void execute();
};


#endif //CDBGGWAS_STATISTICAL_TEST_H
