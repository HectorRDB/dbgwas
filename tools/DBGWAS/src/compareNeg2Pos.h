//
// Created by Leandro Ishi Soares de Lima on 27/06/18.
//

#ifndef DBGWAS_COMPARENEG2POS_H_H
#define DBGWAS_COMPARENEG2POS_H_H

#include <iostream>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/misc/svg.h>

using namespace seqan;

typedef FragmentStore<> TStore;
typedef String<char> TSequence;                             // sequence type
typedef StringSet<TSequence> TStringSet;                    // container for strings
typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
typedef seqan::Graph<Alignment<TDepStringSet> > TAlignGraph;       // alignment graph
typedef Value<TStore::TContigStore>::Type  TContig;

// function declarations
TStore makeAssembly(TStore store);
bool findExact(const TSequence &text, const TSequence &pattern);
TSequence selectContig(const TStore &store, const TStringSet &signif);
int homologyScore(TSequence seq1, TSequence seq2);

#endif //DBGWAS_COMPARENEG2POS_H_H
