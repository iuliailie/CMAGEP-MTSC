#include "Gene.h"
#include "Global.h"
#include <string>
#include <iostream>
#include <vector>

using namespace std;
Gene::Gene(string head, string tail)
{
    _head=head;
    _tail=tail;
    this->getStringVectors();
    headLength=headVector.size();
    tailLength=tailVector.size();
    geneString=_head+_tail;
    this->NodeMap=new  map<int,Node>();
}


Gene::Gene()
{
}


Gene  *Gene::copyGene()//we just need a new address because the rest of properties will be computed in the new cycle
{
    Gene *copiedGene=new Gene();
    copiedGene->_head=this->_head;
    copiedGene->_tail=this->_tail;
    copiedGene->headLength=this->headVector.size();
    copiedGene->tailLength=this->tailVector.size();
    copiedGene->geneString=copiedGene->_head+copiedGene->_tail;
    copiedGene->getStringVectors();
    copiedGene->NodeMap=new  map<int,Node>();
    return copiedGene;
}




void Gene::getOrfSize(map<char, FunctionStructure> *localMap)
{
    string geneString=this->geneString;
    ORF=0;
    int orf=1;//we start with an orf size of one, assuming we only have a terminal to begin with

    for(int i=0; i<geneComponentsVector.size(); i++)
    {
        if(isFunction(this->geneComponentsVector.at(i))) //if the current character is a function we increase the orf size with the arity of the function
            //deducting 1 for the function itself whcih can be considered as a parameter for anther function
            orf+=localMap->at(this->geneComponentsVector.at(i).at(0)).Arity-1;
        else
            orf--;//if it's a terminal we deduct the number of needed parameters with one

        if(orf==0)//when we don't need anymore characters to fill the need of parameters for all the function we return the current position in the gene adding one because the count starts at 0
        {
            ORF =i+1;
            break;
        }
    }

    if(ORF<=0)
        ORF=-1;//in case of a problem we return -1
}
string Gene::toString()
{
    return _head+' '+_tail;
}


void Gene::getStringVectors()

{
    split(this->_head,this->headVector,' ');
    split(this->_tail,this->tailVector,' ');
    split(this->toString(),this->geneComponentsVector,' ');
}




void Gene::GenerateNodeMap(map<char,FunctionStructure> *localMap)//we create the binary tree here for the current gene
{
    // this->NodeMap=new  map<int,Node>();
    this->getOrfSize(localMap);
    int sizeOrf=this->ORF;
    string geneString=this->toString();


    std::vector<std::string> geneStringVector= this->geneComponentsVector;
          

    int *Mark=new int[sizeOrf];

    for(int i=0; i<sizeOrf; i++)
        Mark[i]=0;

    //  Node * nodejshd=&this->NodeMap[0];
    for(int i=sizeOrf-1; i>=0; i--)
    {
        if(this->NodeMap->size()<sizeOrf)
        {
            int j=i-1;

            if(!(isFunction(geneStringVector.at(i)))&& Mark[i]==0 && sizeOrf!=1)
            {
                while(j>0)
                {
                    if(isFunction(geneStringVector.at(j))&&(Mark[j]==0))
                        break;

                    j--;
                }

                Node *leftNode, *rightNode, *node;

                if(localMap->at(geneStringVector.at(j).at(0)).Arity==1)
                {
                    leftNode=new Node(geneStringVector.at(i));
                    rightNode=new Node("0");
                    this->NodeMap->insert(std::pair<int,Node>(i,*leftNode));
                    Mark[i]++;
                }
                else
                {
                    leftNode=new Node(geneStringVector.at(i));
                    this->NodeMap->insert(std::pair<int,Node>(i,*leftNode));
                    Mark[i]++;
                    map<int,Node>::iterator it;

                    if((it=this->NodeMap->find((i-1)))==this->NodeMap->end())
                    {
                        rightNode=new Node(geneStringVector.at(i-1));
                        this->NodeMap->insert(std::pair<int,Node>(i-1,*rightNode));
                    }
                    else
                        rightNode=&this->NodeMap->at(i-1);

                    Mark[i-1]++;
                }

                node=new Node(geneStringVector.at(j),rightNode,leftNode);

                if(localMap->at(geneStringVector.at(j).at(0)).Arity==1)
                {
                    node->type=UnOperator;
                }

                Mark[j]++;
                this->NodeMap->insert(std::pair<int,Node>(j,*node));
                delete node,leftNode,rightNode;
                // i--;
            }
            else
                if(isFunction(geneStringVector.at(i)) && Mark[i]==1)
                {
                    while(j>0)
                    {
                        if(isFunction(geneStringVector.at(j))&&(Mark[j]==0))
                            break;

                        j--;
                    }

                    Node *leftNode, *rightNode, *node;

                    if(localMap->at(geneStringVector.at(j).at(0)).Arity==1)
                    {
                        leftNode=&this->NodeMap->at(i);
                        rightNode=new Node("0");
                        Mark[i]++;
                    }
                    else
                        if(sizeOrf>1)
                        {
                            leftNode=&this->NodeMap->at(i);
                            Mark[i]++;
                            map<int,Node>::iterator it;

                            if((it=this->NodeMap->find((i-1)))==this->NodeMap->end())
                            {
                                rightNode=new Node(geneStringVector.at(i-1));
                                this->NodeMap->insert(std::pair<int,Node>(i-1,*rightNode));
                            }
                            else
                                rightNode=&this->NodeMap->at(i-1);

                            Mark[i-1]++;
                        }

                    node=new Node(geneStringVector.at(j),rightNode,leftNode);

                    if(localMap->at(geneStringVector.at(j).at(0)).Arity==1)
                    {
                        node->type=UnOperator;
                    }

                    this->NodeMap->insert(std::pair<int,Node>(j,*node));
                    Mark[j]++;
                    delete node,leftNode,rightNode;
                    // i--;
                }
                else
                    if(sizeOrf==1)
                    {
                        Node *leftNode;
                        leftNode=new Node(geneStringVector.at(0));
                        this->NodeMap->insert(std::pair<int,Node>(0,*leftNode));
                        Mark[0]++;
                        delete leftNode;
                    }
        }
    }

    delete [] Mark;
}



