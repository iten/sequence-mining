#ifndef PROTEIN_MARKOV_MODEL_H_
#define PROTEIN_MARKOV_MODEL_H_
#include <string>
class ProteinMarkovModel {
public:
    ProteinMarkovModel() { prior = 0.01; }
    ProteinMarkovModel(double prior) : prior(prior) {}
    double GetProbability(std::string sequence);
    void SetPrior(double prior);
    bool LoadPositiveMatrix(std::string filename);
    bool LoadNegativeMatrix(std::string filename);
private:
    double transitions_pos[20][20]; // Positive transition matrix
    double transitions_neg[20][20]; // Negative transition matrix
    double prior;
};
#endif // PROTEIN_MARKOV_MODEL_H_
