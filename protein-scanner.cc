#include <stddef.h> // needed for re2.h
#include <stdio.h>
#include <utility>
#include <vector>
#include "re2/re2.h"
#include "protein-scanner.h"
#include "protein-markov-model.h"

using std::string;

void ProteinScanner::Scan(const char *s,
                          std::vector< std::pair<string, double> > *results) {
    re2::StringPiece text(s);
    string match;
    while(RE2::FindAndConsume(&text, protein_regex, &match)) {
        if(!RE2::FullMatch(match, dna_regex)) { // will recognize dna matches
            double prob = detector.GetProbability(match);
            if(prob > threshold) {
                results->push_back(make_pair(match, prob));
            }
        }
    }
}
