#ifndef GENE_H
#define GENE_H

#include "Global.h"
#include <string>
#include <iostream>
#include <vector>
#include<map>

using namespace std;
class Gene
{
public:
    Gene();

    Gene(string head, string tail);
    ~Gene()
    {
        this->headLength=0;
        this->tailLength=0;
        this->ORF=0;
        this->_head.clear();
        this-> _tail.clear();
        this-> geneString.clear();
        this->value=0;
        NodeMap->clear();
        delete this->NodeMap;
    }


    int headLength;
    int tailLength;
    int ORF;
    double value;
    string _head;
    string _tail;
    string geneString;

    vector<string> headVector;
    vector<string> tailVector;
    vector<string> geneComponentsVector;

    map<int,Node> *NodeMap;

    void getOrfSize(map<char, FunctionStructure> *localMap);
    void GenerateNodeMap(map<char,FunctionStructure> *localMap);
    void getStringVectors();

    string toString();

    Gene  *copyGene();

private:

};

#endif // GENE_H
