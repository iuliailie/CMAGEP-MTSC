#include "Global.h"
#include "GEP.h"


using namespace std;

//here we set the data size and number of explanatory variables

const string allPossibleCharsForTermials = "ABDGIJKMNOQRUVWZbdjkuvwz";


string writingFile;
string inputReadFile;
string targetReadFie;
string writingFileStats;
string currentResultFile;


int oneGEPRun(int k, clock_t currentTime, int finali);

int editDistance(string *s, string *t);

int myMin(int a, int b, int c);

void simplifyTestFunc();

void optimizedGep();

double **readInputFromFile(int number_of_fitness_cases, string file);

double *readTargetFromFile(int number_of_fitness_cases, string file);

long seedgen(void) {
    long s, seed, pid;
    pid = getpid();
    s = time(NULL);    /* get CPU seconds since 01/01/1970 */
    seed = abs(((s * 181) * ((pid - 83) * 359)) % 104729);
    cout << "seed " << seed << endl;
    return seed;
}


#include <string>
#include <sstream>
#include <vector>








//
//_x_2^_x_3


//int main()
//
//{
//    srand(time(NULL));
//    cout<<rand()%3+1;
//    cout<<rand()%3+1;
//    cout<<rand()%3+1;
//    cout<<rand()%3+1;
//}

















int main()
{
// {   cout << "You have entered " << argc
//          << " arguments:" << "\n";

    // for (int i = 0; i < argc; ++i)
    //     cout << argv[i] << "\n";


    number_of_fitness_cases=500;
    number_of_parameters=11;
    clock_t t1,t2;
    unsigned long int sec= time(NULL);

    Py_Initialize();//initialize python module
    PyRun_SimpleString(
            "import sys\n"
            "sys.path.append('/Users/iulia/multiclassCMAGEP/')\n"
            "sys.path.append('/Applications/anaconda3/lib/python3.6/site-packages')\n"
    );
    // int i=1;
    string resultRoot="/Users/iulia/multiclassCMAGEP/Results/VideoMTSC/";

    stringstream ss,ssStats;
    ss<<resultRoot+"Result_VideoMTSC_Cl0_0Lags.csv";
    ssStats<<resultRoot+"Stats_VideoMTSC_Cl0_0Lags.txt";
    std::ofstream out(resultRoot+"Out_VideoMTSC_Cl0_0Lags.txt");
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf


    currentResultFile=currentResultFile.substr(0,2);


   // std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    const char *c_name = ss.str().c_str();
    string asda=ss.str();
    writingFile=c_name;
    const char *c_nameStats = ssStats.str().c_str();
    string asZda=ssStats.str();
    writingFileStats=c_nameStats;
    fstream outfile;
    outfile.open(writingFile,ios::out | ios::app);
    outfile<<"fitness,max fitness possible,structure,simplified structure,optimized "
            "structure,tree size,optimizedFitness,type of fitness function,inputFile,"
            "number of generations,fitness evaluations,runtime"<<"\n";
    outfile.close();
    int k=0;
    int simplifyTest=0;

    if(!simplifyTest)
    {
        float diff;
        k=0 ;

        while(k>=0)
        {
            t1=clock();
            oneGEPRun(4,t1,k);
            diff=(double)(clock() - t1)/CLOCKS_PER_SEC;
            k--;
            outfile.open(writingFile,fstream::in | fstream::out | fstream::app);
            outfile<<(int)diff<<"\n";
            outfile.close();
            out.flush();
        }
    }

    Py_Finalize();
    //std::cout.rdbuf(coutbuf); //reset to standard output again
    return 0;
}


int oneGEPRun(int k, clock_t currentTime, int finali)
{

    // int i=1;
    string resultRoot="/Users/iulia/multiclassCMAGEP/Results/VideoMTSC/";

    string currentFolder="folderDeInlocuit";
    string numarDeInl="numarDeInlocuit";
//    cout<<numarDeInl<<endl;
//  cout<<numarDeInl.at(0)<<endl;
//   cout<<numarDeInl.at(0) - '0'<<endl;
    string fileRootLocation="/Users/iulia/UCD/Binh work/";


  
    inputReadFile=fileRootLocation+"x_XTrain_Cl0_0Lags.txt";
    targetReadFie=fileRootLocation+"x_YTrain_Cl0_0Lags.txt";
    //targetReadFie="/Users/iulia/UCD/testNewCMAGEPY.txt";

    cout<<inputReadFile;
    inputData=readInputFromFile(number_of_fitness_cases, inputReadFile);  //true for input file false for targetData file
    targetData=readTargetFromFile(number_of_fitness_cases, targetReadFie);


    //initiate the function alphabet and variables
    //set the parameters for the chromosomes


    terminalsToUse = "";
    for (int indCount = 1; indCount <= number_of_parameters; indCount++) {
        if (terminalsToUse.size() == 0)
            terminalsToUse = "x_" + to_string(indCount);
        else
            terminalsToUse = terminalsToUse + ' ' + "x_" + to_string(indCount);

    }

    int populationSize=100;
    int  headLength=6;
    int nrGene=4;
    char linkF='+';
    int maxNumberOfCycles=100000;
    double maxTimeToRun=600;
    double timeToStartCMA=1000;
    bool useCMA=false;
    double currentTimeRun;
    //terminalsToUse='_x';
    srand(time(NULL));
//    int headLength = 5;
//    int popSize = 10;
    //constants = "1 2 3 4 5 6 7 8 9";
    constants="1 2 3 4 5 6 7 8 9";
    string functii = "+ - * / ^ S L E";
    functionSet = functii;
    terminals= terminalsToUse;
    split(functii, functionsVector, ' ');
    split(terminals, terminalsVector, ' ');
    split(constants, constantsVector, ' ');

    vector<string> dimensionsToUse;
    //split("1 2 3 4 5 6 7 8", dimensionsToUse, ' ');

    if (dimensionsToUse.size()>0)
        for (int i=0;i<dimensionsToUse.size();i++)
        { int currDim=stoi(dimensionsToUse.at(i));
            terminalsToUseVector.push_back(terminalsVector.at(currDim));
        }
    else
        terminalsToUseVector=terminalsVector;


    allElementsVector.insert(allElementsVector.end(), functionsVector.begin(), functionsVector.end() );
    terminalsToUseVector.insert(terminalsToUseVector.end(), constantsVector.begin(), constantsVector.end() );
    allElementsVector.insert(allElementsVector.end(), terminalsToUseVector.begin(), terminalsToUseVector.end() );



//    functionSet="+-*/^SLE";//SELsC ^SLEsC;
//
//    terminals=allPossibleCharsForTermials.substr(0,number_of_parameters);
//    terminalsToUse="ABDGIJKMN";//fara P H F T S L E C g

    //constants="123456789";

    //generating the first population of chrmosomes randomly
    GEP * gep=new GEP();
    vector<Chromosome *> *ChromList= gep->generatePopulation(headLength,nrGene,linkF,populationSize);
    map<char,FunctionStructure> *localMap;
    Global *g=new Global;
    genesHeadLength=headLength;
    numberOfGenesInAChromosome=nrGene;
    g->mapCharactersWithFunctions();
    localMap=g->functionDictionary;
    //set the threshold accepted for difference between the predicted and actual data
    gep->fitnessType=(FitnessType)(k);
    gep->precisionForFitness=0.01;
    //set parameters for the mutation rates for the algorithm
    double mutationValue=(double)2/((2*headLength+1)*nrGene); //equivalent to 2 point mutation per chromosome
    gep->mutationRate=0.2;
    gep->inversionRate=0.05;
    gep->onePRRate=0.4;
    gep->twoPRRate=0.0;
    gep->geneRRate=0.1;
    gep->IS_transpRate=0.05;
    gep->RIS_transpRate=0.1;
    gep->geneTranspRate=0.1;
    //first training cycle where we  evaluate  the expresion trees and the fitnesses in the chromosomes of the first population
    /*saving the settings of the runs*/
    string currentDate="dataDeInlocuit";
    fstream settingsFile;
    string settingsWritingFile;
    string fisierSettings=resultRoot+"SettingsVideoMTSC_Cl0_0Lags.txt";
    //string fisierSettings=resultRoot+"Settings.txt";
    stringstream ss;
    ss<<fisierSettings;
    const char *c_n = ss.str().c_str();
    string asda=ss.str();
    settingsWritingFile=c_n;
    settingsFile.open(settingsWritingFile,fstream::in | fstream::out | fstream::trunc);
    settingsFile<<"samplesize "<<number_of_fitness_cases<<"\n";
    settingsFile<<"populationSize "<<populationSize<<"\n";
    settingsFile<<"functionSet "<<functionSet<<"\n";
    settingsFile<<"terminals "<<terminalsToUse<<"\n";
    settingsFile<<"constants "<<constants<<"\n";
    settingsFile<<"headLength "<<headLength<<"\n";
    settingsFile<<"nrGene "<<nrGene<<"\n";
    settingsFile<<"linkF "<<linkF<<"\n";
    settingsFile<<"maxNumberOfCycles "<<maxNumberOfCycles<<"\n";
    settingsFile<<"maxTimeToRun "<<maxTimeToRun<<"\n";
    settingsFile<<"timeToStartCMA "<<timeToStartCMA<<"\n";
    settingsFile<<"gep->fitnessType "<< gep->fitnessType<<"\n";
    settingsFile<<"gep->mutationRate "<<gep->mutationRate<<"\n";
    settingsFile<<"gep->geneTranspRate "<<gep->geneTranspRate<<"\n";
    settingsFile<<"gep->inversionRate "<<gep->inversionRate<<"\n";
    settingsFile<<"gep->onePRRate "<<gep->onePRRate<<"\n";
    settingsFile<<"gep->twoPRRate "<<gep->twoPRRate<<"\n";
    settingsFile<<"gep->geneRRate "<< gep->geneRRate<<"\n";
    settingsFile<<"gep->IS_transpRate "<<gep->IS_transpRate<<"\n";
    settingsFile<<"gep->RIS_transpRate"<<gep->RIS_transpRate<<"\n";
    settingsFile.close();
    /*ends here*/
    gep->doTrainingCycle(ChromList,localMap,false,true);
    vector<Chromosome *> *SelectedChromosomes;
    //if(currentTime>=timeToStartCMA)useCMA=true;
    SelectedChromosomes=gep->SelectChromosomesOnOptimizedFitness(ChromList,useCMA);//initial attribution of fitnesses, no optimization on this step useCMA
    delete ChromList;
    double currentMaxFitness=1e50;
    double MaxUnimproovedGenerationsAllowed=30;//we set how Few generations can have the previous fitness before stop
    int NumberOfUnimproovedGenerations=0;
    //!!---------------begin GEP
    int j;

    for(j=0; j<maxNumberOfCycles; j++)
    {
//        cout<<"cycle "<<j<<endl;
        clock_t initTimeCycle=clock();

        if(SelectedChromosomes->at(0)->optimizedFitness<=SelectedChromosomes->at(0)->maxFitnessPossible)
        {
            //cout<<"best fitness found at run "<<j<<endl;
            clock_t endTimeCycle=clock();
            gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
            cout<<"--------------------"<<j<<endl;
            break;
        }

        for(unsigned int i=1; i<SelectedChromosomes->size(); i++)
        {
            if(getARandom(0.0,1.0)<gep->mutationRate)
            {
                SelectedChromosomes->at(i)->Mutate();//
                break;
            }

            if(getARandom(0.0,1.0)<gep->inversionRate)
            {
                SelectedChromosomes->at(i)->doInversion();
                break;
            }

            int other=rand()%(SelectedChromosomes->size()-1)+1;

            if(getARandom(0.0,1.0)<gep->geneRRate)
            {
                SelectedChromosomes->at(i)->GeneSwap(SelectedChromosomes->at(other));
                break;
            }

            if(getARandom(0.0,1.0)<gep->IS_transpRate)
            {
                SelectedChromosomes->at(i)->doISTransposition();
                break;
            }

            if(getARandom(0.0,1.0)<gep->RIS_transpRate)
            {
                SelectedChromosomes->at(i)->doRootTransposition();
                break;
            }

            other=rand()%(SelectedChromosomes->size()-1)+1;

            if(getARandom(0.0,1.0)<gep->onePRRate)
            {
                SelectedChromosomes->at(i)->doOnePointRecombination(SelectedChromosomes->at(other));
                break;
            }

            other=rand()%(SelectedChromosomes->size()-1)+1;

//            if(getARandom(0.0,1.0)<gep->twoPRRate)
//            {
//                SelectedChromosomes->at(i)->doTwoPointRecombination(SelectedChromosomes->at(other));
//                break;
//            }
        }

        gep->doTrainingCycle(SelectedChromosomes,localMap,false,false);
        vector<Chromosome *> * SelectedChromosomesOld=SelectedChromosomes;
        currentTimeRun=(double)((clock() - currentTime)/CLOCKS_PER_SEC);

        if(currentTimeRun>=timeToStartCMA)
            useCMA=true;

        SelectedChromosomes=gep->SelectChromosomesOnOptimizedFitness(SelectedChromosomesOld,useCMA);

        if(SelectedChromosomesOld!=SelectedChromosomes)
        {
            for(int l=0; l<SelectedChromosomesOld->size(); l++)
            {
                delete SelectedChromosomesOld->at(l);
            }

            delete SelectedChromosomesOld;
            //free(SelectedChromosomesOld);
        }//we delete the previous generation from memory

        currentTimeRun=(double)((clock() - currentTime)/CLOCKS_PER_SEC);
        cout<<"currentTimeRun "<<currentTimeRun<<endl;
        if(maxTimeToRun<=currentTimeRun)
        {
            clock_t endTimeCycle=clock();
            gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
            break;
        }

        if(currentMaxFitness>SelectedChromosomes->at(0)->optimizedFitness)//we check how Few generations we have without improvement
        {
            currentMaxFitness=SelectedChromosomes->at(0)->optimizedFitness;
            NumberOfUnimproovedGenerations=0;
        }
        else
        {
            NumberOfUnimproovedGenerations++;

            if(NumberOfUnimproovedGenerations>MaxUnimproovedGenerationsAllowed)
            {
                //cout<<"best fitness found at run "<<j<<endl;
                cout<<"--------------------"<<j<<endl;
                clock_t endTimeCycle=clock();
                gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
                break;
            }
        }

        if(SelectedChromosomes->at(0)->optimizedFitness<=SelectedChromosomes->at(0)->maxFitnessPossible)
        {
            cout<<"best fitness found at run "<<j<<endl;
            cout<<"--------------------"<<j<<endl;
            clock_t endTimeCycle=clock();
            gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
            break;
        }

        if(j==maxNumberOfCycles)
        {
            cout<<"best fitness found at run "<<j<<endl;
            cout<<"--------------------"<<j<<endl;
        }

        //cout<<j<<endl;
        clock_t endTimeCycle=clock();
        gep->statistics->at(j).at(3)=(double)(endTimeCycle - initTimeCycle)/CLOCKS_PER_SEC;
        currentTimeRun=(double)((clock() - currentTime)/CLOCKS_PER_SEC);

        if(currentTimeRun>200)
        {
            gep->doTrainingCycle(SelectedChromosomes,localMap,false,false);
            stringstream ssHist;
            ssHist<<resultRoot+"History"+to_string(finali)+".csv";
            const char *c_nameHist = ssHist.str().c_str();
            fstream outfileHist ;
            outfileHist .open(c_nameHist,fstream::in | fstream::out | fstream::app);

            for(int indexH=0; indexH<10; indexH++)
            {
                outfileHist<<SelectedChromosomes->at(indexH)->mathExpression<<" , "<<SelectedChromosomes->at(indexH)->fitness;
                outfileHist<<"\n";
            }
        }
    }

    //for optimizing the last generation
   vector<Chromosome *> * SelectedChromosomesOld=SelectedChromosomes;
   SelectedChromosomes=gep->SelectChromosomesOnOptimizedFitness(SelectedChromosomesOld,true);

   if(SelectedChromosomesOld!=SelectedChromosomes)
   {
       for(int l=0; l<SelectedChromosomesOld->size(); l++)
       {
           delete SelectedChromosomesOld->at(l);
       }

       delete SelectedChromosomesOld;
    }//we delete the previous generation from memory

    //!!-----------end GEP
//    cout<<"main 264";
    fstream outfile;
    outfile.open(writingFile,fstream::in | fstream::out | fstream::app);
    gep->doTrainingCycle(SelectedChromosomes,localMap,true,false);
    outfile<<SelectedChromosomes->at(0)->fitness<<","<<SelectedChromosomes->at(0)->maxFitnessPossible<<",";
    outfile<<SelectedChromosomes->at(0)->mathExpression;
    outfile<<","<< SelectedChromosomes->at(0)->simplifiedMathExpression;
    outfile<<","<<SelectedChromosomes->at(0)->optimizedMathExpression;

   if(SelectedChromosomes->at(0)->optimizedParameters->size()>0)
            outfile<<","<<SelectedChromosomes->at(0)->optimizedParameters->size();//save tree size
        else
            outfile<<","<<SelectedChromosomes->at(0)->OrfSum;

    outfile<<","<<SelectedChromosomes->at(0)->optimizedFitness;
    outfile<<","<<gep->fitnessType;
    outfile<<","<<inputReadFile;
    outfile<<","<<(j+1);//saving number of generations in the results file
    outfile<<","<<gep->fitFuncCalls+NBCalls<<",";//saving number of fitness eval in the results file
    outfile.close();
    fstream outfileStats;
    outfileStats.open(writingFileStats,fstream::in | fstream::out | fstream::app);

    for(int indexstats=0; indexstats<=gep->statistics->size(); indexstats++)
    {
        outfileStats<<currentResultFile<<" "<<finali<<" "<<indexstats+1<<" ";

        if(indexstats<gep->statistics->size())
            for(int istats=0; istats<gep->statistics->at(indexstats).size(); istats++)
                outfileStats<<gep->statistics->at(indexstats).at(istats)<<" ";
        else
        {
            outfileStats<<SelectedChromosomes->at(0)->optimizedFitness<<" ";

            for(int istats=1; istats<gep->statistics->at(indexstats-1).size(); istats++)
                outfileStats<<gep->statistics->at(indexstats-1).at(istats)<<" ";
        }

        outfileStats<<endl;
    }

    outfileStats<<"\n"<<endl;
    outfileStats.close();

    /*//for testing purposes
        cout<<"-----------------"<<endl;
        double *pointer=SelectedChromosomes->at(0)->Values;
        for(int i=0;i<number_of_fitness_cases;i++)
        {
        cout<<endl<<*(pointer+i);
        }
        cout<<"-----------------"<<endl;
        */
    for(int i=0; i< SelectedChromosomes->size(); i++)
    {
        SelectedChromosomes-> at(i)->optimizedParameters->clear();

        // delete  SelectedChromosomes-> at(i)->optimizedParameters;
    }

    SelectedChromosomes->clear();
    gep->statistics->clear();
    delete gep, g;
    delete SelectedChromosomes;

    for (int i = 0; i < number_of_fitness_cases; ++i) {
        delete inputData[i];
        // each i-th pointer is now pointing to dynamic array (size number_of_fitness_cases) of actual double values
    }

    delete inputData;

    cout<<"GEP run done"<<endl;

    return 0;
}


double ** readInputFromFile(int  number_of_fitness_cases, string file)
{
    srand(seedgen());
    vector<string> STRING;
    ifstream infile;
    infile.open(file);
    vector<string> readStrings;
    string s;
    double ** inputData;

    inputData=new double *[number_of_fitness_cases];


    for (int i = 0; i < number_of_fitness_cases; ++i) {
        inputData[i] = new double[number_of_fitness_cases];
        // each i-th pointer is now pointing to dynamic array (size number_of_fitness_cases) of actual double values
    }



    int t=0;

    while(infile.good()) // To get you all the lines.
    {
        STRING.push_back("");
        getline(infile,STRING.at(t)); // Saves the line in STRING.
        t++;
    }

    infile.close();


    for(int i=0; i<number_of_fitness_cases; i++)
    {
        istringstream f(STRING.at(i));

           while(getline(f, s, '\r'))
            {
                readStrings=split(s.c_str(),'\t');
            }

            //cout<<readStrings.size()<<endl;

            for(int l=0; l<number_of_parameters; l++)
            {
                //cout<<readStrings.at(l)<<endl;
                inputData[i][l]=atof(readStrings.at(l).c_str());
            }
    }

    return inputData;
}


double * readTargetFromFile(int  number_of_fitness_cases, string file)
{
    srand(seedgen());
    vector<string> STRING;
    ifstream infile;
    infile.open(file);
    vector<string> readStrings;
    string s;
    double * targetData;
    targetData=new double[number_of_fitness_cases];
    int t=0;

    while(infile.good()) // To get you all the lines.
    {
        STRING.push_back("");
        getline(infile,STRING.at(t)); // Saves the line in STRING.
        t++;
    }

    infile.close();
    //int k= getARandom((int)0,187);
    int k=0;

    for(int i=k; i<k+number_of_fitness_cases; i++)
    {
        istringstream f(STRING.at(i));


        while(getline(f, s, '\r'))
        {
            readStrings.push_back(s);
        }

        //targetData[(i)]=log(atof(readStrings.at(i).c_str()));
        targetData[(i)]=atof(readStrings.at(i).c_str());
    }

    return targetData;
}






//int editDistance(string *s, string *t)
//{
//    int m,n;
//    m=s->length();
//    n=t->length();
//    int d[m+1][n+1];
//
//    for(int i=0; i<=m; i++)
//        d[i][0]=i;
//
//    for(int j=0; j<=n; j++)
//        d[0][j]=j;
//
//    for(int j=1; j<=n; j++)
//        for(int i=1; i<=m; i++)
//        {
//            if(s->at(i-1)==t->at(j-1))
//                d[i][j]=d[i-1][j-1];
//            else
//                d[i][j]=myMin(d[i-1][j] + 1,  // a deletion
//                              d[i][j-1] + 1,  // an insertion
//                              d[i-1][j-1] + 1 // a substitution
//                             );
//        }
//
////    cout<<endl;
////
////    for(int i=0; i<=m; i++)
////    {
////        for(int j=0; j<=n; j++)
////            cout<< d[i][j]<<" ";
////
////        cout<<endl;
////    }
//    cout<<endl<<"distanta este "<<d[m][n];
//    return d[m][n];
//}
//
//
//int myMin(int a, int b, int c)
//{
//    int temp=1e5;
//
//    if(temp>a)
//        temp=a;
//
//    if(temp>b)
//        temp=b;
//
//    if(temp>c)
//        temp=c;
//
////    cout<<endl<<a<<" "<<b<<" "<<c<<" "<<endl;
////    cout<<temp<<endl;
//    return temp;
//}
//
//
//
//
///*int argc, char **argv
//string inputFile,targetFile,currentPath,parameterSet,sympy_path;
//
//
//    int opt;
//    while ((opt = getopt(argc, argv, "X:Y:f:N:P:p:s")) != -1) {
//        switch(opt) {
//            case 'X':
//                inputFile = string (optarg);
//                break;
//            case 'Y':
//                targetFile = string(optarg);
//                break;
//            case 'f':
//                currentPath=string(optarg);
//                break;
//            case 'N':
//                number_of_fitness_cases = atoi(optarg);
//                break;
//            case 'P':
//                number_of_parameters = atoi(optarg);
//                break;
//            case 'p':
//                parameterSet = string(optarg);
//                break;
//            case 's':
//                sympy_path = string (optarg);
//                break;
//            default:
//                std::cout << "Usage: " << argv[0] << std::endl;
//                return -1;
//        }
//    }
//*/
//


//int main2() {
//
//
//    int noVars = 5;
//    terminalsToUse = "";
//    for (int indCount = 1; indCount <= noVars; indCount++) {
//        if (terminalsToUse.size() == 0)
//            terminalsToUse = "x_" + to_string(indCount);
//        else
//            terminalsToUse = terminalsToUse + ' ' + "x_" + to_string(indCount);
//
//    }
//
//    //terminalsToUse='_x';
//    srand(time(NULL));
//    int headLength = 5;
//    int popSize = 10;
//    constants = "2 3 5 6";
//    string functii = "+ - * / L E s C";
//    functionSet = functii;
//    string terminale = terminalsToUse + ' ' + constants;
//
//    split(functii, functionsVector, ' ');
//    split(terminale, terminalsVector, ' ');
//    split(functii + ' ' + terminale, allElementsVector, ' ');
//
//
//    vector<Chromosome *> *ChromList = new vector<Chromosome *>();
//
//    while (popSize > 0) {
//        int nrGene = 2;
//        Chromosome *crom = new Chromosome();
//        crom->numberOfGenes = nrGene;
//        crom->linkFunction = '+';
//
//        map<char, FunctionStructure> *localMap;
//        Global *g = new Global;
//        g->mapCharactersWithFunctions();
//        localMap = g->functionDictionary;
//
//        while (nrGene > 0) {  //while we still have genes to generate
//            string head = "", tail = "";//reinitialise working tail and head
//            int indexGeneConstruct = 0;//indexGeneConstruct for head and tail string
//
//            while (indexGeneConstruct < headLength) //while the head isnt't filled with values
//            {
//
//
//                if (head.size() <=
//                    0)//if it's the first position in head we must fill it with a function to make sure we get diversity
//                    // head.push_back((functii).at(rand()%functii.size()));
//                {
//                    head = head + functionsVector.at(rand() % functionsVector.size());
//                    indexGeneConstruct++;
//                } else {
//                    head = head + ' ' + (allElementsVector).at(rand() % allElementsVector.size());
//                    indexGeneConstruct++;
//
//                }
//            }
//
//
//            indexGeneConstruct = 0;
//
//            while (indexGeneConstruct <
//                   (headLength + 1))//the largest arity is 2 so t=h+1 from t=h*(n-1)+1 with n largest arity
//            {
//                if (tail.size() == 0)
//                    tail = tail + (terminalsVector).at(rand() % terminalsVector.size());
//                else tail = tail + ' ' + (terminalsVector).at(rand() % terminalsVector.size());
//                indexGeneConstruct++;
//
//            }
//
//
//            // cout << head << "---" << tail << endl;
//            Gene *gene = new Gene(head, tail);//we create a new gene with the head and tail recently generated
//
////            genesHeadLength=headLength;
////            numberOfGenesInAChromosome=nrGene;
//
//            gene->getOrfSize(localMap);
//            gene->GenerateNodeMap(localMap);
//
//            vector<string> elems = gene->geneComponentsVector;
//
//            crom->Genes->push_back(gene);//we add the new gene to the vector of genes of the Chromosome
//            nrGene--;
//        }
//
//
//        if (ChromList->size() > 2) {
//            crom->getChromosomeMathExpression();
//
//            string expressi1 = crom->mathExpression;
//            crom->doTwoPointRecombination(ChromList->at(0));
//
//            for (int i = 0; i < crom->numberOfGenes; i++) {
//                crom->Genes->at(i)->getOrfSize(localMap);
//                crom->Genes->at(i)->GenerateNodeMap(localMap);
//            }
//
//            crom->getChromosomeMathExpression();
//            string expressi2 = crom->mathExpression;
//            cout << expressi1 << " " << expressi2;
//
//        }
//
//        cout << endl;
//        ChromList->push_back(crom);
//        popSize--;
//    }
//
//    //return ChromList;
//
//}


