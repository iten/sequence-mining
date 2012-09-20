#include <stddef.h> // needed for re2.h
#include <stdio.h>
#include "re2/re2.h"
#include "protein-scanner.h"
#include "protein-markov-model.h"

using std::string;

void ProteinScanner::Scan(const char *s) {
    re2::StringPiece text(s);
    string match;
    while(RE2::FindAndConsume(&text, protein_regex, &match)) {
        if(!RE2::FullMatch(match, dna_regex)) { // will recognize dna matches
            double prob = detector.GetProbability(match.data());
            if(prob > 0.5) {
                printf("%s %lf\n", match.data(),
                       detector.GetProbability(match.data()));
            }
        }
    }
}
