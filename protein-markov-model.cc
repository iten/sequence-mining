#include <string>
#include <stdio.h>
#include "protein-markov-model.h"

using std::string;

namespace {
int AminoAcidToIndex(char amino_acid) {
    switch(amino_acid) {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'D':
    case 'd':
        return 2;
    case 'E':
    case 'e':
        return 3;
    case 'F':
    case 'f':
        return 4;
    case 'G':
    case 'g':
        return 5;
    case 'H':
    case 'h':
        return 6;
    case 'I':
    case 'i':
        return 7;
    case 'K':
    case 'k':
        return 8;
    case 'L':
    case 'l':
        return 9;
    case 'M':
    case 'm':
        return 10;
    case 'N':
    case 'n':
        return 11;
    case 'P':
    case 'p':
        return 12;
    case 'Q':
    case 'q':
        return 13;
    case 'R':
    case 'r':
        return 14;
    case 'S':
    case 's':
        return 15;
    case 'T':
    case 't':
        return 16;
    case 'V':
    case 'v':
        return 17;
    case 'W':
    case 'w':
        return 18;
    case 'Y':
    case 'y':
        return 19;
    default:
        return -1;
    }
}
}

bool ProteinMarkovModel::LoadPositiveMatrix(string filename) {
    FILE *f = fopen(filename.c_str(), "r");
    if(f == NULL) {
        perror("Opening positive matrix failed");
        return false;
    }
    for(int i = 0; i < 20; ++i) {
        for(int j = 0; j < 20; ++j) {
            double tmp;
            int ret = fscanf(f, "%lf\n", &tmp);
            if(ret == 0) {
                perror("Failed to parse matrix from file");
                return false;
            }
            if(ret == EOF) {
                fprintf(stderr, "Matrix file not long enough\n");
                return false;
            }
            transitions_pos[i][j] = tmp;
        }
    }
    fclose(f);
    return true;
}

bool ProteinMarkovModel::LoadNegativeMatrix(string filename) {
    FILE *f = fopen(filename.c_str(), "r");
    if(f == NULL) {
        perror("Opening negative matrix failed");
        return false;
    }
    for(int i = 0; i < 20; ++i) {
        for(int j = 0; j < 20; ++j) {
            double tmp;
            int ret = fscanf(f, "%lf\n", &tmp);
            if(ret == 0) {
                perror("Failed to parse matrix from file");
                return false;
            }
            if(ret == EOF) {
                fprintf(stderr, "Matrix file not long enough\n");
                return false;
            }
            transitions_neg[i][j] = tmp;
        }
    }
    fclose(f);
    return true;
}

void ProteinMarkovModel::SetPrior(double in) {
    prior = in;
}

// Note: we can safely assume the string is ASCII and not UTF-8 at this point,
// since we've already eliminated any other characters.
double ProteinMarkovModel::GetProbability(string sequence) {
    double prob_seq_pos = 1.0;
    double prob_seq_neg = 1.0;
    // Start iteration from the second character since we look at the previous
    // character to form the transition probability.
    for(size_t i = 1; i < sequence.size(); ++i) {
        int prev_idx = AminoAcidToIndex(sequence[i - 1]);
        int curr_idx = AminoAcidToIndex(sequence[i]);
        if((prev_idx == -1) | (curr_idx == -1)) {
            fprintf(stderr, "WARNING: unable to handle transition %c|%c\n",
                    sequence[i - 1], sequence[i]);
            return 0.0;
        }
        double prob_transition_pos = transitions_pos[prev_idx][curr_idx];
        double prob_transition_neg = transitions_neg[prev_idx][curr_idx];
        prob_seq_pos *= prob_transition_pos;
        prob_seq_neg *= prob_transition_neg;
    }
    double prob_seq = prob_seq_pos*prior + prob_seq_neg*(1.0-prior);
    return prob_seq_pos*prior/prob_seq;
}
