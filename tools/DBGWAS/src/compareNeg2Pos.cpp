#include <compareNeg2Pos.h>

using namespace seqan;

void makeAssembly(TStore &store) {
    ConsensusAlignmentOptions options;
    options.useContigID = false;
    consensusAlignment(store, options);
}

bool findExact(const TSequence &text, const TSequence &pattern){ 
// Brute force pattern matching for every position
// from a seqAn example... (I tried t use the text.find(pattern) function, but it is not applicabel to TSequence objects...)

    bool found = false;

    // Number of consecutive matching characters per position
    String<int> score;
    resize(score, length(text) - length(pattern) + 1);

    for (unsigned i = 0; i < length(text) - length(pattern) + 1; ++i)
    {
        int localScore = 0;
        for (unsigned j = 0; j < length(pattern); ++j)
            if (text[i + j] == pattern[j])
                ++localScore;
        score[i] = localScore;
    }

    for (unsigned i = 0; i < length(score); ++i)
        if (score[i] == (int)length(pattern))
            found = true;

    return found;
}

TSequence selectContig(const TStore &store, const TStringSet &signif) { 
    int nb_contigs = length(store.contigStore);
    // std::cout << "NB of contigs: " << nb_contigs << std::endl;

    if(nb_contigs == 1){
        // single contig ==> no need to make a choice
        TSequence contig = store.contigStore[0].seq;
        // std::cout << "Length of returned (unique) contigs: " << length(contig) << std::endl;
        return contig;
    }else if(nb_contigs < 1){
        // Here we have a problem ==> I should return an error...
        return 0;
    }else{
        // choose the best contig to represent the Neg/Pos population
        // Rule : filter ou contigs without significant unitigs, and take the longer...
        TSequence contig = "A"; // initialize with the shortest !

        for (int i = 0; i < nb_contigs ; i++)
        {
            TSequence seq = store.contigStore[i].seq;
            bool valid = false;

            if(length(signif) == 0){
                valid = true; // there is no significant unitig : consider all contigs...
                // std::cout << "No signf. : all contigs are candidated" << std::endl;
            }else{
                // check if there is a match betwwen a signif. sequence and the contig:
                for (int j =0; j < length(signif) ; j++)
                {
                    bool found = findExact(seq,signif[j]);
                    if(found)
                    {
                        valid = true;
                        // std::cout << "Signif\t" << j << "\twas found in conig\t"<< i << std::endl;
                    }else{
                        // std::cout << "No match found between contig\t" << i << "\tand signif\t" << j << std::endl;
                    }
                } 
            }

            if(valid)
            {
                // std::cout << "Current contig\t" << i << "\t of length\t" << length(seq) << "\twhile current max =\t" << length(contig) << std::endl;
                if(length(seq) > length(contig)){
                    contig = seq; 
                }
            }
        }
        // std::cout << "Length of returned contigs (among several):\t" << length(contig) << std::endl;
        return contig;
    }
}

int homologyScore(const TSequence &seq1, const TSequence &seq2){
    TStringSet homologySet;
    appendValue(homologySet, seq1);
    appendValue(homologySet, seq2);

    TAlignGraph alignG(homologySet);     
    int score = globalAlignment(alignG, Score<int, Simple>(1, -1, -1), AlignConfig<true, true, true, true>());         
    std::cout << alignG << std::endl;

    return score;  
}