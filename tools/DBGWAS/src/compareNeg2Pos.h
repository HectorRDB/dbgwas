//
// Created by Leandro Ishi Soares de Lima on 27/06/18.
//

#ifndef DBGWAS_COMPARENEG2POS_H_H
#define DBGWAS_COMPARENEG2POS_H_H

#include <iostream>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/misc/svg.h>

typedef seqan::FragmentStore<> TStore;
typedef seqan::String<char> TSequence;                             // sequence type
typedef seqan::StringSet<TSequence> TStringSet;                    // container for strings
typedef seqan::StringSet<TSequence, seqan::Dependent<> > TDepStringSet;   // dependent string set
typedef seqan::Graph<seqan::Alignment<TDepStringSet> > TAlignGraph;       // alignment graph
typedef seqan::Value<TStore::TContigStore>::Type  TContig;

// function declarations
void makeAssembly(TStore &store);
bool findExact(const TSequence &text, const TSequence &pattern);
TSequence selectContig(const TStore &store, const TStringSet &signif);
int homologyScore(const TSequence &seq1, const TSequence &seq2);

#endif //DBGWAS_COMPARENEG2POS_H_H
