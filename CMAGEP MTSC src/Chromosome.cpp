#include "Chromosome.h"

#include <string>
#include <iostream>
#include <vector>

using namespace std;
//bool operator< (Chromosome &cC1, Chromosome &cC2)
//{
//    return cC1.optimizedFitness <=cC2.optimizedFitness;
//}

Chromosome::Chromosome(int geneNb, char _linkFunction) {
    numberOfGenes = geneNb;
    linkFunction = _linkFunction;
    Genes = new std::vector<Gene *>();
    Values = new double[number_of_fitness_cases];
    optimizedParameters = new std::vector<double>();
}

Chromosome::Chromosome() {
    Genes = new std::vector<Gene *>();
    Values = new double[number_of_fitness_cases];
    optimizedParameters = new std::vector<double>();
}

Chromosome *Chromosome::copyChromosome() {
    Chromosome *crom = new Chromosome();
    crom->numberOfGenes = this->numberOfGenes;
    crom->linkFunction = this->linkFunction;
    crom->fitness = this->fitness;
    crom->OrfSum = this->OrfSum;
    crom->optimizedFitness = this->optimizedFitness;
    crom->optimizedMathExpression = this->optimizedMathExpression;
    crom->mathExpression = this->mathExpression;
    crom->simplifiedMathExpression = this->simplifiedMathExpression;
    crom->maxFitnessPossible = this->maxFitnessPossible;
    crom->optimizedParameters = this->optimizedParameters;

    for (int i = 0; i < crom->numberOfGenes; i++) {
        Gene *copiedGene = 0;
        copiedGene = this->Genes->at(i)->copyGene();
        crom->Genes->push_back(copiedGene);
    }

    return crom;
}

void Chromosome::computeChromosomeValue(double *p, map<char, FunctionStructure> *localMap, int j) {
    // cout<<"for value="<<p[0];
    // cout<<endl;
    vector<double> vector1, vector2, vector3;//initialize working vectors
    this->OrfSum = 0;

//
//    for(int i=0;i<3;i++)
//    {
////        for (int tip=0;tip<3;tip++)
////        { i+=tip;
//            cout<<*(p+i)<<" ";
//
////        }
//
//    cout<<endl;
//    }
//    cout<<endl<<"done"<<endl;

    for (int CromIterator = 0; CromIterator < this->numberOfGenes; CromIterator++) {
        Gene *gene;
        gene = (this->Genes->at(CromIterator));
        Node *node = &(gene->NodeMap->at(0));
        //cout<<gene->geneString<<endl;
        double val = computeTreeValue(node, p, localMap);

        if (node->param != 1) {
            val = node->param * val;
            cout << node->param << "----" << endl;
        }

        this->OrfSum = this->OrfSum + gene->ORF;

        if (vector1.size() == 0)//if  we are at the first gene, the chromosome value is equal to the gene value
            vector1.insert(vector1.begin(), val);
        else//else we apply the linking function to the values from the former genes and the current one
        {
            vector2.insert(vector2.begin(), val);
            vector3.insert(vector3.begin(), 0);
            transform(vector1.begin(), vector1.end(), vector2.begin(), vector3.begin(),
                      localMap->at(this->linkFunction).pFunction);
            vector1.at(0) = vector3.at(0);
        }
    }

//    if(this->Values->size()==nr)
//        this->Values->clear();
//cout<<endl<<vector1.at(0);
    this->Values[j] = vector1.at(0);
    vector1.clear();
    vector2.clear();
    vector3.clear();
}


void Chromosome::Mutate() {
    string functii = functionSet;
    string terminale = terminals;

//    split(functii,functionsVector,' ');
//    split(terminale,terminalsVector,' ');
//    split(functii+' '+terminale,allElementsVector,' ');

    // select a random gene
    int geneNoToMutate = rand() % (this->numberOfGenes);
    // select a random position in the gene
    vector<string> tempHeadVector = this->Genes->at(geneNoToMutate)->headVector;
    vector<string> tempTailVector = this->Genes->at(geneNoToMutate)->tailVector;


    // in the head we can change to a function or a terminal
    int indexInGeneToMutate = rand() % (tempHeadVector.size() + tempTailVector.size());

    string temp;

    if (indexInGeneToMutate < this->Genes->at(geneNoToMutate)->headVector.size()) {
        // note that the first position in the head must be a function


        if (indexInGeneToMutate == 0) {
            // mutate to a (different) function because the head can only remain a function
            temp = functionsVector.at(rand() % (functionsVector.size()));
            tempHeadVector.at(indexInGeneToMutate) = temp;
            //(index,1,1,temp);//replace for length 1 for 1 time starting from index
        } else {
            // mutate to a (different) function or a terminal
            temp = allElementsVector.at(rand() % ((allElementsVector).size()));

            tempHeadVector.at(indexInGeneToMutate) = temp;
        }
    } else {
        // in the tail we can only mutate to a terminal

        temp = terminalsVector.at(rand() % terminalsVector.size());

        tempTailVector.at(indexInGeneToMutate - tempHeadVector.size()) = temp;
    }

    string tempHead, tempTail;

    tempHead = getStringFromVector(tempHeadVector);
    tempTail = getStringFromVector(tempTailVector);

    Gene *old = this->Genes->at(geneNoToMutate);
    this->Genes->at(geneNoToMutate) = new Gene(tempHead, tempTail);

    delete old;
    // return true;.
}


void Chromosome::doInversion() {
    //randomly choose for which gene will the inversion be performed

    int geneNoToModify = rand() % (this->numberOfGenes);

    // select a random position in the gene
    vector<string> tempHeadVector = this->Genes->at(geneNoToModify)->headVector;
    vector<string> tempTailVector = this->Genes->at(geneNoToModify)->tailVector;
    vector<string> tempStringVector = this->Genes->at(geneNoToModify)->geneComponentsVector;

    //select subsection to reverse

    int i1 = rand() % this->Genes->at(geneNoToModify)->geneComponentsVector.size();
    int i2 = rand() % this->Genes->at(geneNoToModify)->geneComponentsVector.size();

    std::vector<string> subStrVector;
    for (int i = i1; i < i2; i++)
        subStrVector.push_back(tempStringVector.at(i));

    std::reverse(subStrVector.begin(), subStrVector.end());

    int subInd = 0;
    for (int i = i1; i < i2; i++)
        tempStringVector.at(i) = subStrVector.at(subInd++);

    for (int j = 0; j < tempHeadVector.size(); j++)
        tempHeadVector.at(j) = tempStringVector.at(j);

    for (int j = tempHeadVector.size(); j < tempTailVector.size(); j++)
        tempTailVector.at(j - tempHeadVector.size()) = tempStringVector.at(j);

    string tempHead, tempTail;

    tempHead = getStringFromVector(tempHeadVector);
    tempTail = getStringFromVector(tempTailVector);

    Gene *old = this->Genes->at(geneNoToModify);
    this->Genes->at(geneNoToModify) = new Gene(tempHead, tempTail);

    delete old;
}


void Chromosome::doISTransposition() {

    //inserting genetic sequence in head

    // select a random gene


    int geneNoToModify = rand() % (this->numberOfGenes);

    // select a random position in the gene


    int k = rand() % this->numberOfGenes;

    Gene *randomGene = this->Genes->at(geneNoToModify);

    vector<string> tempHeadVector = randomGene->headVector;
    vector<string> tempTailVector = randomGene->tailVector;
    vector<string> tempStringVector = randomGene->geneComponentsVector;


    unsigned long i1 = rand() % tempHeadVector.size();//seq. to insert start
    unsigned long i2 = rand() % (tempStringVector.size() - tempHeadVector.size() - 1) + 1;// seq. to insert length
    unsigned long insP;
    if (tempHeadVector.size() - (i2) > 0)
        insP = rand() % (tempHeadVector.size() - (i2));
    else insP = 0;;//insertion point

    std::vector<string> subStrVector;
    //generate new vector with sequence to insert
    for (int i = i1; i < i1 + i2; i++)
        subStrVector.push_back(tempStringVector.at(i));


    int subInd = 0;
    for (int i = insP; i < insP + i2; i++)//introduce new seq. in gene
        tempStringVector.at(i) = subStrVector.at(subInd++);


    for (int j = 0; j < tempHeadVector.size(); j++)
        tempHeadVector.at(j) = tempStringVector.at(j);

    for (int j = tempHeadVector.size(); j < tempTailVector.size(); j++)
        tempTailVector.at(j - tempHeadVector.size()) = tempStringVector.at(j);


    string tempHead, tempTail;

    tempHead = getStringFromVector(tempHeadVector);
    tempTail = getStringFromVector(tempTailVector);

    Gene *old = this->Genes->at(geneNoToModify);
    this->Genes->at(geneNoToModify) = new Gene(tempHead, tempTail);

    delete old;


}



void Chromosome::doRootTransposition() {

    int geneNoToModify = rand() % (this->numberOfGenes);

    // select a random position in the gene


    int k = rand() % this->numberOfGenes;

    Gene *randomGene = this->Genes->at(geneNoToModify);

    vector<string> tempHeadVector = randomGene->headVector;
    vector<string> tempTailVector = randomGene->tailVector;
    vector<string> tempStringVector = randomGene->geneComponentsVector;


    unsigned long i1 = rand() % tempHeadVector.size();//seq. to insert start
    unsigned long i2 = rand() % (tempStringVector.size() - tempHeadVector.size() - 1) + 1;// seq. to insert length
    unsigned long insP;

    if (tempHeadVector.size() - (i2) > 0)
        insP = rand() % (tempHeadVector.size() - (i2));
    else insP = 0;;//insertion point

    std::vector<string> subStrVector;
    //generate new vector with sequence to insert
    for (int i = i1; i < i1 + i2; i++)
        subStrVector.push_back(tempStringVector.at(i));


    int subInd = 0;
    for (int i = insP; i < insP + i2; i++)//introduce new seq. in gene
        tempStringVector.at(i) = subStrVector.at(subInd++);


    for (int j = 0; j < tempHeadVector.size(); j++)
        tempHeadVector.at(j) = tempStringVector.at(j);

    for (int j = tempHeadVector.size(); j < tempTailVector.size(); j++)
        tempTailVector.at(j - tempHeadVector.size()) = tempStringVector.at(j);


    string tempHead, tempTail;

    tempHead = getStringFromVector(tempHeadVector);
    tempTail = getStringFromVector(tempTailVector);

    Gene *old = this->Genes->at(geneNoToModify);
    this->Genes->at(geneNoToModify) = new Gene(tempHead, tempTail);

    delete old;


}


void Chromosome::GeneSwap(Chromosome *Chrom2) {
    // Choose a random gene from both chromosomes, and swap them
    int g1 = rand() % this->numberOfGenes;
    int g2 = rand() % Chrom2->numberOfGenes;
    Gene *gene1 = new Gene(this->Genes->at(g1)->_head, this->Genes->at(g1)->_tail);
    Gene *gene2 = new Gene(Chrom2->Genes->at(g2)->_head, Chrom2->Genes->at(g2)->_tail);
    Gene *oldGene1, *oldGene2;
    oldGene1 = this->Genes->at(g1);
    oldGene2 = Chrom2->Genes->at(g2);
    Chrom2->Genes->at(g2) = gene1;
    this->Genes->at(g1) = gene2;
    gene1 = 0;
    gene2 = 0;
    delete gene1, gene2;
    delete oldGene1, oldGene2;
}

void Chromosome::GeneSwap(Chromosome *Chrom2, int geneNb) {
    // Choose a random gene from both chromosomes, and swap them
    int g1 = geneNb;
    int g2 = geneNb;
    Gene *gene1 = new Gene(this->Genes->at(g1)->_head, this->Genes->at(g1)->_tail);
    Gene *gene2 = new Gene(Chrom2->Genes->at(g2)->_head, Chrom2->Genes->at(g2)->_tail);
    delete this->Genes->at(g1), Chrom2->Genes->at(g2);
    Chrom2->Genes->at(g2) = gene1;
    this->Genes->at(g1) = gene2;
    gene1 = 0;
    gene2 = 0;//clean the pointers
    delete gene1, gene2;
}


void Chromosome::doOnePointRecombination(Chromosome *otherChrome) {

    int g = getARandom(0, this->numberOfGenes - 1);//choosing the gene
    int rP = getARandom(0, this->Genes->at(0)->geneComponentsVector.size() - 1);//choose crossover point

    vector<string> geneS1V = this->Genes->at(g)->geneComponentsVector;
    vector<string> geneS2V = otherChrome->Genes->at(g)->geneComponentsVector;

    vector<string> newStringVec1;
    vector<string> newStringVec2;



    for (int i=0;i<rP;i++)
    {
        newStringVec1.push_back(geneS1V.at(i));
        newStringVec2.push_back(geneS2V.at(i));
    }

    for (int i=rP;i<geneS2V.size();i++)
    {
        newStringVec1.push_back(geneS2V.at(i));
        newStringVec2.push_back(geneS1V.at(i));
    }


    vector<string> geneS1VHead = this->Genes->at(g)->headVector;
    vector<string> geneS1VTail = this->Genes->at(g)->tailVector;

    vector<string> geneS2VHead = otherChrome->Genes->at(g)->headVector;
    vector<string> geneS2VTail = otherChrome->Genes->at(g)->tailVector;



    for (int i=0;i<geneS1VHead.size();i++)
    {
        geneS1VHead.at(i)=newStringVec1.at(i);
        geneS2VHead.at(i)=newStringVec2.at(i);
    }

    for (int i=geneS1VHead.size();i<geneS1VHead.size();i++)
    {
        geneS1VTail.at(i)=newStringVec1.at(i);
        geneS2VTail.at(i)=newStringVec2.at(i);
    }


    string tempHead, tempTail;

    tempHead = getStringFromVector(geneS1VHead);
    tempTail = getStringFromVector(geneS1VTail);

    Gene *oldG1 = this->Genes->at(g);
    Gene *oldG2 = otherChrome->Genes->at(g);

    this->Genes->at(g) = new Gene(getStringFromVector(geneS1VHead), getStringFromVector(geneS1VTail));
    otherChrome->Genes->at(g)=new Gene(getStringFromVector(geneS2VHead), getStringFromVector(geneS2VTail));

    delete oldG1,oldG2;


}
void Chromosome::doTwoPointRecombination(Chromosome *other) {
    vector<string> newStringVec1;
    vector<string> newStringVec2;

    //we generate a string which will contain all the gene
    // strings for ease of double crossover

    for (int i = 0; i <this->numberOfGenes; i++)
    {
        for (int indVec=0; indVec<this->Genes->at(i)->geneComponentsVector.size();indVec ++) {
            newStringVec1.push_back(this->Genes->at(i)->geneComponentsVector.at(indVec));
            newStringVec2.push_back(other->Genes->at(i)->geneComponentsVector.at(indVec));
        }

        //fullString2 += other->Genes->at(i)->geneString;
    }


}

void Chromosome::doTwoPointRecombinationOld(Chromosome *other) {
    string fullString1 = "";
    string fullString2 = "";

    for (int i = 0; i <this->numberOfGenes; i++) //we generate a string which will contain all the genes' string for ease of double crossover
    {
        fullString1 += this->Genes->at(i)->geneString;
        fullString2 += other->Genes->at(i)->geneString;
    }

    int p1 = getARandom(0, fullString1.length()-2);//we choose the two crossover points, the first can be at most the char before the last
    int p2 = getARandom(p1 + 1, fullString1.length() - 1);
    int hL = this->Genes->at(0)->headLength;//get headLength for future use
    string substr1 = fullString1.substr(p1,
                                        p2 - p1);//get the substrings which will be exchanged betweent the two chromes
    string substr2 = fullString2.substr(p1, p2 - p1);
    fullString1.replace(p1, p2 - p1, substr2);
    fullString2.replace(p1, p2 - p1, substr1);

    for (int i = 0; i < this->numberOfGenes; i++) {
        delete this->Genes->at(i), other->Genes->at(i);
        string sH1, sH2, sT1, sT2;
        sH1 = fullString1.substr((2 * hL + 1) * i, hL + (2 * hL + 1) * i - (2 * hL + 1) * i);
        sT1 = fullString1.substr(hL + (2 * hL + 1) * i, 2 * hL + 1 + (2 * hL + 1) * i - (hL + (2 * hL + 1) * i));
        this->Genes->at(i) = new Gene(sH1, sT1);
        sH2 = fullString2.substr((2 * hL + 1) * i, hL + (2 * hL + 1) * i - (2 * hL + 1) * i);
        sT2 = fullString2.substr(hL + (2 * hL + 1) * i, 2 * hL + 1 + (2 * hL + 1) * i - (hL + (2 * hL + 1) * i));
        other->Genes->at(i) = new Gene(sH2, sT2);
    }



}

void Chromosome::getChromosomeMathExpression() {
    string ChromExpression = "";

    for (int k = 0; k < this->numberOfGenes; k++) {
        string geneExpression = printExpresie(&this->Genes->at(k)->NodeMap->at(0));

        if (k < this->numberOfGenes - 1)
            ChromExpression += geneExpression + this->linkFunction;
        else
            ChromExpression += geneExpression;
    }

    this->mathExpression = ChromExpression;
}





