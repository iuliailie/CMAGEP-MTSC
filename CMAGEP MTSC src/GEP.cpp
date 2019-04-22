#include "GEP.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <fstream>

using namespace std;

GEP::GEP()
{
    this->statistics=new vector<std::array<double, 7>>();
}

//vector<Chromosome *>  *GEP::generatePopulation2(int HeadLength,
//        int geneNumber, char linkingFunction ,  int populationSize)
//{
//    int  headLength=HeadLength;
//    int popSize=populationSize;
//    string functii=functionSet;
//    string terminale=terminalsToUse+constants;
//    vector<Chromosome *> *ChromList=new vector<Chromosome *>();
//
//    while(popSize>0)
//    {
//        int nrGene=geneNumber;
//        Chromosome *crom=new Chromosome();
//        crom->numberOfGenes=geneNumber;
//        crom->linkFunction=linkingFunction;
//
//        while(nrGene>0)//while we still have genes to generate
//        {
//            string head="",tail="";//reinitialise working tail and head
//            int index=0;//index for head and tail string
//
//            while(index<headLength) //while the head isnt't filled with values
//            {
//                if(head.size()<=0)//if it's the first position in head we must fill it with a function to make sure we get diversity
//                    head.push_back((functii).at(rand()%functii.size()));
//                else
//                    head.push_back((functii+terminale).at(rand()%((terminale+functii).size())));//else we can choose all values from functions and terminals
//
//                index++;
//            }
//
//            index=0;
//
//            while(index<headLength+1)//the biggest arity is 2 so t=h+1 from t=h*(n-1)+1
//            {
//                tail.push_back((terminale).at(rand()%terminale.size()));//we fill the tail only with values from the terminals string
//                index++;
//            }
//
//            Gene *gene=new Gene(head,tail);//we create a new gene with the head and tail recently generated
//            crom->Genes->push_back(gene);//we add the new gene to the vector of genes of the Chromosome
//            nrGene--;
//        }
//
//        ChromList->push_back(crom);
//        popSize--;
//    }
//
//    return ChromList;
//}


vector<Chromosome *>  *GEP::generatePopulation(int HeadLength,
                                               int geneNumber, char linkingFunction ,  int populationSize)
{
    int  headLength=HeadLength;
    int popSize=populationSize;
    string functii=functionSet;
    string terminale=terminalsToUse+constants;
    vector<Chromosome *> *ChromList=new vector<Chromosome *>();

    while(popSize>0)
    {
        int nrGene=geneNumber;
        Chromosome *crom=new Chromosome();
        crom->numberOfGenes=geneNumber;
        crom->linkFunction=linkingFunction;

        while (nrGene > 0) {  //while we still have genes to generate
            string head = "", tail = "";//reinitialise working tail and head
            int indexGeneConstruct = 0;//indexGeneConstruct for head and tail string

            while (indexGeneConstruct < headLength) //while the head isnt't filled with values
            {


                if (head.size() <=
                    0)//if it's the first position in head we must fill it with a function to make sure we get diversity
// head.push_back((functii).at(rand()%functii.size()));
                {
                    head = head + functionsVector.at(rand() % functionsVector.size());
                    indexGeneConstruct++;
                } else {
                    head = head + ' ' + (allElementsVector).at(rand() % allElementsVector.size());
                    indexGeneConstruct++;

                }
            }


            indexGeneConstruct = 0;

            while (indexGeneConstruct <
                   (headLength + 1))//the largest arity is 2 so t=h+1 from t=h*(n-1)+1 with n largest arity
            {
                if (tail.size() == 0)
                    tail = tail + (terminalsToUseVector).at(rand() % terminalsToUseVector.size());
                else tail = tail + ' ' + (terminalsToUseVector).at(rand() % terminalsToUseVector.size());
                indexGeneConstruct++;

            }


// cout << head << "---" << tail << endl;
            Gene *gene = new Gene(head, tail);
//we create a new gene with the head and tail recently generated
            crom->Genes->push_back(gene);//we add the new gene to the vector of genes of the Chromosome
            nrGene--;
        }

        ChromList->push_back(crom);
        popSize--;
    }

    return ChromList;
}




int GEP::randSelector(double totalFitness,vector<Chromosome *> *ChromList)//random wheighted selection
{
    //we avoid selecting the chromosomes which have fitness=0
    double randFitness=getARandom(1.0,totalFitness);
    int currentFitness=0;

    for(unsigned k=0; k<ChromList->size(); k++)
    {
        currentFitness+=ChromList->at(k)->fitness;

        if(currentFitness<=randFitness)
        {
            return k;
        }
    }

    return 0;
}

int GEP::randSelectorTournament(vector<Chromosome *> *ChromList)//tournament selection
{
    //we avoid selecting the chromosomes which have fitness=0
    double c1=getARandom(0,ChromList->size()-1);
    double c2=getARandom(0,ChromList->size()-1);

    if(ChromList->at(c1)->fitness<=ChromList->at(c2)->fitness)
    {
        return c1;
    }
    else
        return c2;

    return 0;
}


int GEP::randSelectorOpt(double totalFitness,vector<Chromosome *> *ChromList)//random wheighted selection
{
    //we avoid selecting the chromosomes which have fitness=0
    double randFitness=getARandom(1.0,totalFitness);
    int currentFitness=0;

    for(unsigned k=0; k<ChromList->size(); k++)
    {
        currentFitness+=ChromList->at(k)->optimizedFitness;

        if(currentFitness<=randFitness)
        {
            return k;
        }
    }

    return 0;
}

int GEP::randSelectorTournamentOpt(vector<Chromosome *> *ChromList)//tournament selection
{
    //we avoid selecting the chromosomes which have fitness=0
    double c1=getARandom(0,ChromList->size()-1);
    double c2=getARandom(0,ChromList->size()-1);

    if(ChromList->at(c1)->optimizedFitness<ChromList->at(c2)->optimizedFitness)
        return c1;


    return c2;
    return 0;
}




vector<Chromosome *> *GEP::SelectChromosomes(vector<Chromosome *> *ChromList)//selecting the new generation of chromosomes before mutations
{
    vector<Chromosome *>  * SelectedChromes=new vector<Chromosome *>;
    double bestValue=1000;//for saving the best fitted chrome;
    double worstValue=0;//for saving the worst fitted chorme;
    double average=0,minComplexity=1e50,maxComplexity=0,averageComplexity=0;
    int orf=0;
    int k=0;
    double totalFitness=0;//total fitness of all chromosomes

    for(unsigned int i=0; i<ChromList->size(); i++)
    {
        if(bestValue>ChromList->at(i)->fitness)
        {
            bestValue=ChromList->at(i)->fitness;
            orf=ChromList->at(i)->OrfSum;
            k=i;
        }
        else
            if(bestValue==ChromList->at(i)->fitness && ChromList->at(i)->OrfSum<=orf)
            {
                bestValue=ChromList->at(i)->fitness;
                k=i;
                orf=ChromList->at(i)->OrfSum;
            }

        if(worstValue<ChromList->at(i)->fitness)
        {
            worstValue=ChromList->at(i)->fitness;
        }

        totalFitness+=ChromList->at(i)->fitness;

        //complexityStats
        if(minComplexity>ChromList->at(i)->OrfSum)
        {
            minComplexity=ChromList->at(i)->OrfSum;
        }

        if(maxComplexity<ChromList->at(i)->OrfSum)
        {
            maxComplexity=ChromList->at(i)->OrfSum;
        }

        averageComplexity+=ChromList->at(i)->OrfSum;
        //  cout<<i<<" has fitness "<<ChromList->at(i)->fitness<<endl;
    }

    average=totalFitness/ChromList->size();
    averageComplexity/=ChromList->size();
//saving one run stats
    array<double, 7> oneStepStats;
    oneStepStats[0]=bestValue;
    oneStepStats[1]=worstValue;
    oneStepStats[2]=average;
    oneStepStats[3]=0;
    oneStepStats[4]=minComplexity;
    oneStepStats[5]=maxComplexity;
    oneStepStats[6]=averageComplexity;
    this->statistics->push_back(oneStepStats);

    cout<<ChromList->at(k)->fitness<<endl;

    Chromosome *copiedChrom=0;
    copiedChrom=ChromList->at(k)->copyChromosome();
    SelectedChromes->push_back(copiedChrom);//we save the best fitted chromosome in the first position of the newly selected generation

//here we begin the random weighted selection
    for(unsigned int i=1; i<ChromList->size(); i++)
    {
        //  int sel=randSelector(totalFitness,ChromList);
        int sel=randSelectorTournament(ChromList);
        Chromosome *copiedChromNew=0;
        copiedChromNew=ChromList->at(sel)->copyChromosome();
        SelectedChromes->push_back(copiedChromNew);

    }

    return SelectedChromes;
}



vector<Chromosome *> *GEP::SelectChromosomesOnOptimizedFitness(vector<Chromosome *> *ChromList, bool useOptimization)//selecting the new generation of chromosomes before mutations
{
    vector<double> *vectorGol=new vector<double>();
    vector<Chromosome *>  * SelectedChromes=new vector<Chromosome *>;
    double bestValue=1e50;//for saving the best fitted chrome;
    double worstValue=0;//for saving the worst fitted chorme;
    double averageFitness=0, minComplexity=1e50,maxComplexity=0,averageComplexity=0;
    int orf=0;
    int k=0;
    double totalFitness=0;//total optimizedFitness of all chromosomes
    int n=ChromList->size();
    //!! sorting list of chromosomes according to their fitnesses
    sort(ChromList->begin(),ChromList->end(),Cmp());
    //saving one run stats
    bestValue=ChromList->at(0)->optimizedFitness;
    worstValue=ChromList->at(ChromList->size()-1)->optimizedFitness;

    for(int i=0; i<ChromList->size(); i++)
    {
        averageFitness+= ChromList->at(i)->optimizedFitness;

        //complexityStats
        if(minComplexity>ChromList->at(i)->OrfSum)
        {
            minComplexity=ChromList->at(i)->OrfSum;
        }

        if(maxComplexity<ChromList->at(i)->OrfSum)
        {
            maxComplexity=ChromList->at(i)->OrfSum;
        }

        averageComplexity+=ChromList->at(i)->OrfSum;
    }

    averageFitness=averageFitness/ChromList->size();
    averageComplexity=averageComplexity/ChromList->size();
    array<double, 7> oneStepStats;
    oneStepStats[0]=bestValue;
    oneStepStats[1]=worstValue;
    oneStepStats[2]=averageFitness;
    oneStepStats[3]=0;
    oneStepStats[4]=minComplexity;
    oneStepStats[5]=maxComplexity;
    oneStepStats[6]=averageComplexity;
    this->statistics->push_back(oneStepStats);

    int indice=0;

    if(useOptimization)
    {
               for(indice=0; indice<0.05*ChromList->size();indice++)
            {
                 //  ChromList->at(indice)->mathExpression= "((x_8)/(7))+sqrt((x_5))+sqrt((x_16)) ";
               ChromList->at(indice)->simplifiedMathExpression=simplifySolution(ChromList->at(indice)->mathExpression);
                cout<<ChromList->at(indice)->simplifiedMathExpression<<" "<<ChromList->at(indice)->optimizedFitness<<endl;

                if( ChromList->at(indice)->simplifiedMathExpression.at(0)!='!')
                {
                    tuple<string,double, vector<double >*> tup;
                    cout<<ChromList->at(indice)->simplifiedMathExpression<<" "<<ChromList->at(indice)->optimizedFitness;

                    tuple<string,int,vector<double>*> functAndParam=addParametersInPosition(ChromList->at(indice)->simplifiedMathExpression);

                    if(ChromList->at(indice)->optimizedParameters->size()>0)
                    {
                        tup=CallCma2(functAndParam,ChromList->at(indice)->simplifiedMathExpression,
                                     (int)this->fitnessType,ChromList->at(indice)->optimizedParameters);
                    }
                    else

                        tup= CallCma2(functAndParam,ChromList->at(indice)->simplifiedMathExpression,
                                      (int)this->fitnessType,vectorGol);//altfel se poate ca structura sa fie diferita prin evolutie si parametrii sa nu se mai potriveasca

                    if(get<1>(tup)<=ChromList->at(indice)->optimizedFitness)
                    {
                        ChromList->at(indice)->optimizedMathExpression=get<0>(tup);
                        ChromList->at(indice)->optimizedFitness=get<1>(tup);
                        ChromList->at(indice)->optimizedParameters=get<2>(tup);
                    }
                    else
                    {
                        ChromList->at(indice)->optimizedMathExpression= ChromList->at(indice)->simplifiedMathExpression;
                    }

                    cout<<ChromList->at(indice)->optimizedMathExpression<<"  "<<ChromList->at(indice)->optimizedFitness<<endl;
                }
                
                else
                    {
                        ChromList->at(indice)->optimizedMathExpression= ChromList->at(indice)->mathExpression;
                        ChromList->at(indice)->simplifiedMathExpression=ChromList->at(indice)->mathExpression;
                    }

                    ChromList->at(indice)->maxFitnessPossible=-1e06;
                    cout<<"\n ";
            }

           sort(ChromList->begin(),ChromList->begin()+indice,Cmp());
    }
    else
    {
        ChromList->at(0)->optimizedMathExpression= ChromList->at(indice)->mathExpression;
        ChromList->at(0)->simplifiedMathExpression=ChromList->at(indice)->mathExpression;
    }

    Chromosome *copiedChrom=0;
    copiedChrom=ChromList->at(0)->copyChromosome();//!! saving the fist chromosome accroding to its fitness beacuse Chrom list is sorted
    SelectedChromes->push_back(copiedChrom);//we save the best fitted chromosome in the first position of the newly selected generation
    cout<<"------------    "<<SelectedChromes->at(0)->optimizedFitness<<"     -----best fitness-------------"<<endl;
    cout<<"------------    "<<SelectedChromes->at(0)->optimizedMathExpression<<"     ----------------------"<<endl;

//here we begin the random weighted selection

    for(unsigned int i=1; i<ChromList->size(); i++)
    {
        int sel=randSelectorTournamentOpt(ChromList);
        Chromosome *copiedChromNew=0;
        copiedChromNew=ChromList->at(sel)->copyChromosome();
        SelectedChromes->push_back(copiedChromNew);
    }

    delete vectorGol;
    return SelectedChromes;
}
vector<Chromosome *> *GEP::doOptimizationOnly(vector<Chromosome *> *ChromList, bool useOptimization)//selecting the new generation of chromosomes before mutations
{
    vector<double> *vectorGol=new vector<double>();
    vector<Chromosome *>  * SelectedChromes=new vector<Chromosome *>;
    double bestValue=1e50;//for saving the best fitted chrome;
    double worstValue=0;//for saving the worst fitted chorme;
    double averageFitness=0, minComplexity=1e50,maxComplexity=0,averageComplexity=0;
    int orf=0;
    int k=0;
    double totalFitness=0;//total optimizedFitness of all chromosomes
    int n=ChromList->size();
    //!! sorting list of chromosomes according to their fitnesses
    sort(ChromList->begin(),ChromList->end(),Cmp());
    //saving one run stats
    bestValue=ChromList->at(0)->optimizedFitness;
    worstValue=ChromList->at(ChromList->size()-1)->optimizedFitness;

    for(int i=0; i<ChromList->size(); i++)
    {
        averageFitness+= ChromList->at(i)->optimizedFitness;

        //complexityStats
        if(minComplexity>ChromList->at(i)->OrfSum)
        {
            minComplexity=ChromList->at(i)->OrfSum;
        }

        if(maxComplexity<ChromList->at(i)->OrfSum)
        {
            maxComplexity=ChromList->at(i)->OrfSum;
        }

        averageComplexity+=ChromList->at(i)->OrfSum;
    }

    averageFitness=averageFitness/ChromList->size();
    averageComplexity=averageComplexity/ChromList->size();
    array<double, 7> oneStepStats;
    oneStepStats[0]=bestValue;
    oneStepStats[1]=worstValue;
    oneStepStats[2]=averageFitness;
    oneStepStats[3]=0;
    oneStepStats[4]=minComplexity;
    oneStepStats[5]=maxComplexity;
    oneStepStats[6]=averageComplexity;
    this->statistics->push_back(oneStepStats);
//!! working with the 10 best chromosomes
    int indice=0;
//        vector<string> *geneComponentsVector=new vector<string>();
//                string STRING;
//                ifstream infiles;
//                infiles.open("/Users/iuliailie/CMAGEP-v2/Gecco/toOptimize.txt");
//
//                while(infiles.good()) // To get you all the lines.
//                {
//                    getline(infiles,STRING);
//                    // Saves the line in STRING.
//                    geneComponentsVector->push_back(STRING);
//
//                }
    /* vector<string> *geneComponentsVector=new vector<string>();
        string STRING;
        ifstream infiles;
        infiles.open("/Users/iuliailie/CMAGEP-v2/FisierRezultateResp.txt");

        while(infiles.good()) // To get you all the lines.
        {
            getline(infiles,STRING);
            // Saves the line in STRING.
            geneComponentsVector->push_back(STRING);

        } */

    if(useOptimization)
    {
        /* string OptResultFile="currentNode";
        OptResultFile=OptResultFile.substr (0,2);
        stringstream ss;
        ss<<"/Users/iuliailie/CMAGEP-v2/Results/OptimResults"<<OptResultFile<<".txt";
        const char *c_name = ss.str().c_str();
        string  writingFile=c_name;
        fstream outfile; */
        for(indice=0; indice<=5; indice++)
        {
            //ChromList->at(indice)->mathExpression="3.48*A^2+2.5*B";  //for standard results optimizationsd
            //for standard results optimizations
            // ChromList->at(indice)->mathExpression=geneComponentsVector->at(indice);
            ChromList->at(indice)->simplifiedMathExpression=simplifySolution(ChromList->at(indice)->mathExpression);

            if(ChromList->at(indice)->simplifiedMathExpression.at(0)!='!')
            {
                tuple<string,double, vector<double >*> tup;
                cout<<ChromList->at(indice)->simplifiedMathExpression<<" "<<ChromList->at(indice)->optimizedFitness;
                tuple<string,int,vector<double>*> functAndParam=addParametersInPosition(ChromList->at(indice)->simplifiedMathExpression);

                if(ChromList->at(indice)->optimizedParameters->size()>0)
                {
//                      for( int indexParam=0;indexParam<ChromList->at(indice)->optimizedParameters->size();indexParam++)
//                    {
//                        cout<<ChromList->at(indice)->optimizedParameters->at(indexParam)<<" ";
//                    }
                    tup=CallCma2(functAndParam,ChromList->at(indice)->simplifiedMathExpression,
                                 (int)this->fitnessType,ChromList->at(indice)->optimizedParameters);//pastrez parametrii doar pentru primul cromozom
                }
                else
                    tup= CallCma2(functAndParam,ChromList->at(indice)->simplifiedMathExpression,
                                  (int)this->fitnessType,vectorGol);//altfel se poate ca structura sa fie diferita prin evolutie si parametrii sa nu se mai potriveasca

                if(get<1>(tup)<=ChromList->at(indice)->optimizedFitness)
                {
                    ChromList->at(indice)->optimizedMathExpression=get<0>(tup);
                    ChromList->at(indice)->optimizedFitness=get<1>(tup);
                    ChromList->at(indice)->optimizedParameters=get<2>(tup);
                }
                else
                {
                    ChromList->at(indice)->optimizedMathExpression= ChromList->at(indice)->simplifiedMathExpression;
                }

                cout<<get<0>(tup)<<"  "<<get<1>(tup)<<endl;
            }
            else
            {
                ChromList->at(indice)->optimizedMathExpression= ChromList->at(indice)->mathExpression;
                ChromList->at(indice)->simplifiedMathExpression=ChromList->at(indice)->mathExpression;
            }

            ChromList->at(indice)->maxFitnessPossible=-1e06;
            cout<<"\n ";
            /*  outfile.open(writingFile,fstream::in | fstream::out | fstream::app);
             outfile<< ChromList->at(indice)->optimizedMathExpression;
             outfile<< "\n";
             outfile.close(); */
        }

        /*  outfile.open(writingFile,fstream::in | fstream::out | fstream::app);
        outfile<< "\n";
        outfile.close(); */
//throw 25;
        sort(ChromList->begin(),ChromList->begin()+indice,Cmp());
    }
    else
    {
        ChromList->at(0)->optimizedMathExpression= ChromList->at(indice)->mathExpression;
        ChromList->at(0)->simplifiedMathExpression=ChromList->at(indice)->mathExpression;
    }

    Chromosome *copiedChrom=0;
    copiedChrom=ChromList->at(0)->copyChromosome();//!! saving the fist chromosome accroding to its fitness beacuse Chrom list is sorted
    SelectedChromes->push_back(copiedChrom);//we save the best fitted chromosome in the first position of the newly selected generation
    cout<<"------------    "<<SelectedChromes->at(0)->optimizedFitness<<"     -----best fitness-------------"<<endl;
    cout<<"------------    "<<SelectedChromes->at(0)->optimizedMathExpression<<"     ----------------------"<<endl;
    /*
        for(int i=0; i<ChromList->at(k)->numberOfGenes; i++)
        {
            cout<<ChromList->at(k)->Genes->at(i)->toString()<<endl;
        }
    *///for testing purposes

//here we begin the random weighted selection

    for(unsigned int i=1; i<ChromList->size(); i++)
    {
        int sel=randSelectorTournamentOpt(ChromList);
        Chromosome *copiedChromNew=0;
        copiedChromNew=ChromList->at(sel)->copyChromosome();
        SelectedChromes->push_back(copiedChromNew);
    }

    delete vectorGol;
    return SelectedChromes;
}

//
//vector<Chromosome *>  *GEP::generateChromListFromStrings(int HeadLength,
//        int geneNumber, char linkingFunction ,  vector<string>* stringsList)
//{
//    string functii=functionSet;
//    string terminale=terminals+constants;
//    vector<Chromosome *> *ChromList=new vector<Chromosome *>();
//    int TotalGeneLength=HeadLength*2+1;
//    int nbOfStrings=stringsList->size();
//
//
//    while(nbOfStrings>0)
//    {
//        nbOfStrings--;
//          int thisGeneNb=0;
//        while(thisGeneNb<geneNumber)
//        {
//            head=stringsList
//
//
//            Gene *gene=new Gene(head,tail);
//        }
//
//
//
//
//    }
//}
