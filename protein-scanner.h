#ifndef PROTEIN_SCANNER_H_
#define PROTEIN_SCANNER_H_
#include <stddef.h> // needed for re2.h
#include <string>
#include <vector>
#include <utility>
#include "re2/re2.h"
#include "protein-markov-model.h"

class ProteinScanner {
public:
    ProteinScanner(const char *regex, const char *pos_matrix_file,
                   const char *neg_matrix_file, 
                   double prior, double detection_threshold) :
    protein_regex(regex),
    dna_regex("[TGCAtgca]+"),
    threshold(detection_threshold) {
        detector.SetPrior(prior);
        detector.LoadPositiveMatrix(pos_matrix_file);
        detector.LoadNegativeMatrix(neg_matrix_file);
    }
void Scan(const char *s,
          std::vector< std::pair<std::string, double> > *results);
private:
    re2::RE2 protein_regex;
    re2::RE2 dna_regex;
    double threshold;
    ProteinMarkovModel detector;
};
#endif // PROTEIN_SCANNER_H_
