#ifndef PROTEIN_SCANNER_H_
#define PROTEIN_SCANNER_H_
#include <stddef.h> // needed for re2.h
#include "re2/re2.h"
#include "protein-markov-model.h"

class ProteinScanner {
public:
    ProteinScanner(const char *regex, const char *pos_matrix_file,
                   const char *neg_matrix_file, 
                   double prior) : protein_regex(regex), dna_regex("[TGCAtgca]+") {
        detector.SetPrior(prior);
        detector.LoadPositiveMatrix(pos_matrix_file);
        detector.LoadNegativeMatrix(neg_matrix_file);
    }
    void Scan(const char *s);
private:
    re2::RE2 protein_regex;
    re2::RE2 dna_regex;
    ProteinMarkovModel detector;
};
#endif // PROTEIN_SCANNER_H_
