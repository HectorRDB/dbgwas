//
// Created by Leandro Ishi Soares de Lima on 30/06/17.
//

#ifndef KSGATB_PHENOCOUNTER_H
#define KSGATB_PHENOCOUNTER_H


//class used to count how many times a unitig is seen overall, in phenotype 0, 1 or NA
class PhenoCounter {
private:
    int pheno0;
    int pheno1;
    int NA;
public:
    PhenoCounter(int pheno0=0, int pheno1=0, int NA=0):pheno0(pheno0), pheno1(pheno1), NA(NA){}
    int getPheno0() const { return pheno0; }
    void increasePheno0(int n) { pheno0+=n; }

    int getPheno1() const { return pheno1; }
    void increasePheno1(int n) { pheno1+=n; }

    int getNA() const { return NA; }
    void increaseNA(int n) { NA+=n; }

    int getTotal() const { return getPheno0() + getPheno1() + getNA(); }
    int getSumOfBothPhenos() const { return getPheno0() + getPheno1(); }
};


#endif //KSGATB_PHENOCOUNTER_H
