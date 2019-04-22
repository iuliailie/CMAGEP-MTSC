#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "Gene.h"
#include "Global.h"

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
class Chromosome
{


//{private:
//   struct PointerCompare {
//      bool operator()(const Chromosome* l, const Chromosome* r) {
//        return *l < *r;
//      }
//    };

public:


    Chromosome(int geneNb,char linkFunction);
    Chromosome();
    ~Chromosome()
    {
        delete [] Values;

        for(int i=0; i<numberOfGenes; i++)
            delete this->Genes->at(i);

        delete Genes;
        //delete optimizedParameters;
        this->numberOfGenes=0;
        this->fitness=0;
        this->optimizedFitness=0;
        this->maxFitnessPossible=0;
        this->linkFunction='\0';
        this->OrfSum=0;
    }

    int numberOfGenes=0;
    double fitness=1e50;
    double optimizedFitness=1e50;

    double maxFitnessPossible=0;

    char linkFunction;

    double * Values;
    vector<double> *optimizedParameters;
    vector<Gene *> *Genes;
    int OrfSum=0;
    void computeChromosomeValue(double *p,map<char,FunctionStructure> *localMap,int j);
    void getChromosomeMathExpression();
    string mathExpression;
    string simplifiedMathExpression;
    string optimizedMathExpression;
    Chromosome *copyChromosome();


    //mutation methods
    void Mutate();
    void doInversion();
    void doISTransposition();
    void doRootTransposition();
    void doOnePointRecombination(Chromosome * otherChrome);
    void doTwoPointRecombination(Chromosome * other);
    void GeneSwap(Chromosome * Chrom2);
    void GeneSwap(Chromosome * Chrom2,int geneNb);

    void doTwoPointRecombinationOld(Chromosome *other);
    //utils
//
//  bool operator< (const Chromosome left) const
//            {
//
//
//                return left.optimizedFitness >= optimizedFitness;
//            }
//
// bool operator<(const Chromosome& chrom) const
//    {
//        return (this->optimizedFitness.compare(chrom.optimizedFitness) == -1);
//    }



//struct by_OptimizedFitness {
//    bool operator()(Chromosome const &a, Chromosome const &b) const {
//        return a.optimizedFitness < b.optimizedFitness;
//    }


};
struct Cmp
{
    bool operator()(Chromosome* p1, Chromosome *p2) const
    {
        if(p1->optimizedFitness != p2->optimizedFitness)
            return p1->optimizedFitness < p2->optimizedFitness;

        // if (p1->y != p2->y) return p1->y < p2->y;
        return false;
    }
};


#endif // CHROMOSOME_H
