#ifndef FITNESSFUNCTION_H
#define FITNESSFUNCTION_H

#include "Global.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

enum FitnessType {NOH=0,rMSE=1,RSRE=2,rRAE=3,ME=4,AIC=5,CM=6,AICx=7,sOS=8,CEM=9,MAE=10};

class FitnessFunction
{
public:
    FitnessType type;
    double *predicted;
    double * actual;
    int sampleSize;
    double precision;
    int selectionRange;
    int programSize;
    int maxProgramSize;
    int minProgramSize;
    double maxFitness;
    double ComputeFitness();
    double NumberOfHitsFitness();
    double  relativeMeanSquaredError();
    double  rootRelativeSquaredError();
    double  RelativeAbsoluteError();
    double ModellingEfficiency();
    double sumOfSquares();
    double AkaikeIndex();
    double ShannonComplexity(double *ODP);
    double * ordinal_pattern_loop();
    double xAkaikeIndex();
    double ComplexityCorrectedEM();
    double MeanAbsoluteError();
    FitnessFunction();

private:

};

#endif // FITNESSFUNCTION_H
