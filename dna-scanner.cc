#include <stddef.h>
#include <stdio.h>
#include <vector>
#include <string>
#include "re2/re2.h"
#include "dna-scanner.h"

using std::string;

// static
// 
void DnaScanner::GetNucleicChars(const string &text, string *output) {
    for(size_t i = 0; i < text.size(); ++i) {
        switch(text[i]) {
        case 'A':
        case 'G':
        case 'T':
        case 'C':
        case 'U':
        case 'a':
        case 'g':
        case 't':
        case 'c':
        case 'u':
            *output += text[i];
            break;
        default:
            break;
        }
    }
}

void DnaScanner::Scan(const char *s, std::vector<string> *results) {
    re2::StringPiece text(s);
    string match;
    string current_seq;
    while(RE2::FindAndConsume(&text, possible_regex, &match)) {
        string dna_only;
        GetNucleicChars(match, &dna_only);
        double dna_fraction = ((float) dna_only.size())/match.size();
        if(dna_fraction >= 0.99999f && dna_only.size() >= min_pure_len) {
            current_seq += dna_only;
        } else if(dna_fraction >= min_dna_fraction &&
                  dna_only.size() >= min_impure_len) {
            current_seq += dna_only;
        } else { 
            if(current_seq.size() > min_output_len) {
                results->push_back(current_seq);
            }
            current_seq = "";
        }
    }
    if(current_seq.size() > min_output_len) {
        results->push_back(current_seq);
    }
}
