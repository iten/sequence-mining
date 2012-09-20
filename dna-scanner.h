#ifndef DNA_SCANNER_H_
#define DNA_SCANNER_H_
#include <stddef.h> // for re2.h
#include "re2/re2.h"
class DnaScanner {
public:
    DnaScanner(const char *regex, size_t min_pure_len, size_t min_impure_len, 
           size_t min_output_len, double min_dna_fraction) 
    : possible_regex(regex), min_pure_len(min_pure_len), 
        min_impure_len(min_impure_len), min_output_len(min_output_len),
        min_dna_fraction(min_dna_fraction) {}
    void Scan(const char *s);
    static void GetNucleicChars(const std::string &s, std::string *output);
private:
    re2::RE2 possible_regex;
    size_t min_pure_len;
    size_t min_impure_len;
    size_t min_output_len;
    double min_dna_fraction;
};
#endif // DNA_SCANNER_H_
