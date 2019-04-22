#include "FitnessFunction.h"

FitnessFunction::FitnessFunction()
{
    //ctor
}


double FitnessFunction::NumberOfHitsFitness()
{
    double *diff=new double[sampleSize];
    double sum=0;

    for(int i=0; i<sampleSize; i++)
    {
        if(!IsNumber(*(predicted+i))||!IsFiniteNumber(*(predicted+i)))
            return INFINITY;
        else
            *(diff+i)=abs(*(predicted+i)-*(actual+i));

        if(*(diff+i)<=precision)
            sum++;
    }

    maxFitness=1e-05;
    delete [] diff;
    return sampleSize-sum;
}

double FitnessFunction::sumOfSquares()
{
    double *diff=new double[sampleSize];
    double sum=0;

    for(int i=0; i<sampleSize; i++)
    {
        if(!IsNumber(*(predicted+i))||!IsFiniteNumber(*(predicted+i)))
            return INFINITY;
        else
            *(diff+i)=pow(*(predicted+i)-*(actual+i),2);

        sum=sum+(1-*(diff+i));
    }

    maxFitness=1e-05;
    delete [] diff;
    return sampleSize-sum;
}
double  FitnessFunction::relativeMeanSquaredError()
{
    double sum=0;
    double fitness=100000;

    for(int j=0; j<sampleSize; j++)
    {
        if(!IsNumber(*(predicted+j))||!IsFiniteNumber(*(predicted+j)))
            return INFINITY;

        sum+=pow((*(predicted+j)-*(actual+j)),2);
    }

    //double Ei=1/sampleSize*sum;
    //fitness=1000*(1/(1+Ei));
    //double fpp=fitness*(1+1/5000*(maxProgramSize-programSize)/(maxProgramSize-minProgramSize));
    // maxFitness=1.0002*1000;
    //return fpp;//with complexity weight
    fitness=sqrt(sum/sampleSize);
    maxFitness=1e-05;
    return fitness;
}
double FitnessFunction::rootRelativeSquaredError()

{
    double mean=0;
    double fitness;
    double sum1=0,sum2=0;

    for(int i=0; i<sampleSize; i++)
    {
        mean+=*(actual+i);
    }

    mean=mean/sampleSize;

    for(int i=0; i<sampleSize; i++)
    {
        if(!IsNumber(*(predicted+i))||!IsFiniteNumber(*(predicted+i)))
            return INFINITY;

        sum1+=pow((*(predicted+i)-*(actual+i)),2);
        sum2+=pow((*(actual+i)-mean),2);
    }

    double Ei=sqrt(sum1/sum2);
    fitness=Ei;
    double fpp=fitness*(1+1/500000000*(maxProgramSize-programSize)/(maxProgramSize-minProgramSize));
    //maxFitness=1.000000002*1000;
    maxFitness=1e-05;
    // return fpp;
//maxFitness=1000;
    return fitness;
}

double  FitnessFunction::RelativeAbsoluteError()
{
    double mean=0;
    double fitness;
    double sum1=0,sum2=0;

    for(int i=0; i<sampleSize; i++)
    {
        if(!IsNumber(*(predicted+i))||!IsFiniteNumber(*(predicted+i)))
            return INFINITY;

        mean+=*(actual+i);
    }

    mean=mean/sampleSize;

    for(int i=0; i<sampleSize; i++)
    {
        sum1+=abs((*(predicted+i)-*(actual+i))/(*(actual+i)));
        sum2+=abs((*(actual+i)-mean)/mean);
    }

    double Ei=sum1/sum2;
    //fitness=1000*(1/(1+Ei));
    fitness=Ei;
    double fpp=fitness*(1+1/5000*(maxProgramSize-programSize)/(maxProgramSize-minProgramSize));
    //maxFitness=1.0002*1000;
    // maxFitness=999.999;
    maxFitness=1e-05;
    //  return fpp;
    return fitness;
}
double FitnessFunction::ModellingEfficiency()
{
    double mean=0;
    double fitness;
    double sum1=0,sum2=0;

    for(int i=0; i<sampleSize; i++)
    {
        if(!IsNumber(*(predicted+i))||!IsFiniteNumber(*(predicted+i)))
            return INFINITY;

        mean+=*(actual+i);
    }

    mean=mean/sampleSize;

    for(int i=0; i<sampleSize; i++)
    {
        sum1+=pow((*(predicted+i)-*(actual+i)),2);
        sum2+=pow((*(actual+i)-mean),2);
    }

    double Ei=sum1/sum2;
    fitness=Ei;
    maxFitness=1e-05;
    // maxFitness=999.999;
    // return fpp;
    double fpp=fitness;
//    cout<<fitness<<endl;
    return fitness;
}
double FitnessFunction::MeanAbsoluteError()
{
    double mean=0;
    double fitness;
    double sum1=0;

    for(int i=0; i<sampleSize; i++)
    {
        if(!IsNumber(*(predicted+i))||!IsFiniteNumber(*(predicted+i)))
            return INFINITY;

        mean+=*(actual+i);
    }

    mean=mean/sampleSize;

    for(int i=0; i<sampleSize; i++)
    {
        sum1+=abs(*(predicted+i)-*(actual+i));
    }

    double mae=sum1/mean;
    fitness=mae;
    maxFitness=1e-05;
    // maxFitness=999.999;
    // return fpp;
    double fpp=fitness;
//    cout<<fitness<<endl;
    return fitness;
}

double FitnessFunction::AkaikeIndex()
{
    double mean=0;
    double fitness;
    double rss=0;

    try
    {
        for(int i=0; i<sampleSize; i++)
        {
            rss+=pow((*(predicted+i)-*(actual+i)),2);
        }

        if((sampleSize-programSize-1)!=0)
            fitness=sampleSize*log(rss/sampleSize)+2*programSize+(2*programSize*(programSize+1))/(sampleSize-programSize-1);
        else
            fitness=sampleSize*log(rss/sampleSize)+2*programSize+(2*programSize*(programSize+1))/(sampleSize-programSize);

        //fitness=-2*log(sum1)+programSize*(log(sampleSize)+1);
        //  fitness=2*programSize+sampleSize*log(sum1);
        //fitness=sampleSize*log(sum1/sampleSize)+(2*sampleSize*programSize)/(sampleSize-programSize-1);//aic
        // fitness=sampleSize*log(sum1/sampleSize)+programSize*log(sampleSize);//bic
    }
    catch
        (...)
    {
        fitness=1000;
    }

    return fitness;
}

double FitnessFunction::xAkaikeIndex()
{
    double mean=0;
    double fitness;
    double rss=0;
    double rssOfMean=0;
    double meanOfObserved=0;
    double programSizeTerm=0;
    double maxParamNumber;
    maxParamNumber=(2*genesHeadLength+1)*numberOfGenesInAChromosome;
    try
    {
        for(int i=0; i<sampleSize; i++)
        {
            rss+=pow((*(predicted+i)-*(actual+i)),2);
            meanOfObserved+=*(actual+i);
        }
        meanOfObserved=meanOfObserved/sampleSize;

         for(int i=0; i<sampleSize; i++)
        {
            rssOfMean+=pow((meanOfObserved-*(actual+i)),2);
        }


        if((sampleSize-maxParamNumber-1)!=0)
            programSizeTerm=(2*programSize+(2*programSize*(programSize+1))/(sampleSize-programSize))/
            ((2*maxParamNumber+(2*maxParamNumber*(maxParamNumber+1))/(sampleSize-maxParamNumber-1)));
             else
             programSizeTerm=(2*programSize+(2*programSize*(programSize+1))/(sampleSize-programSize))/
             ((2*maxParamNumber+(2*maxParamNumber*(maxParamNumber+1))/(sampleSize-maxParamNumber)));

            fitness=sampleSize*log((rss/rssOfMean)/sampleSize)+programSizeTerm;



        //fitness=-2*log(sum1)+programSize*(log(sampleSize)+1);
        //  fitness=2*programSize+sampleSize*log(sum1);
        //fitness=sampleSize*log(sum1/sampleSize)+(2*sampleSize*programSize)/(sampleSize-programSize-1);//aic
        // fitness=sampleSize*log(sum1/sampleSize)+programSize*log(sampleSize);//bic
    }
    catch
        (...)
    {
        fitness=1000;
    }

    return fitness;
}
double FitnessFunction::ComplexityCorrectedEM()
{
    double mean=0;
    double fitness;
    double rss=0;
    double rssOfMean=0;
    double meanOfObserved=0;
    double programSizeTerm=0;
    double maxParamNumber;
    maxParamNumber=(2*genesHeadLength+1)*numberOfGenesInAChromosome;
    try
    {



        if((sampleSize-maxParamNumber-1)!=0)
            programSizeTerm=(2*programSize+(2*programSize*(programSize+1))/(sampleSize-programSize))/
            ((2*maxParamNumber+(2*maxParamNumber*(maxParamNumber+1))/(sampleSize-maxParamNumber-1)));
             else
             programSizeTerm=(2*programSize+(2*programSize*(programSize+1))/(sampleSize-programSize))/
             ((2*maxParamNumber+(2*maxParamNumber*(maxParamNumber+1))/(sampleSize-maxParamNumber)));
            double me;
            me=ModellingEfficiency();
//            cout<<endl;
//            cout<<"meGP="<<me<<endl;
//            cout<<"GPprogramSizeTerm="<<programSizeTerm<<endl;
            double SE;
            SE=ShannonComplexity(ordinal_pattern_loop());
//            cout<<"GPSE="<<SE<<endl;
            fitness=sqrt(pow((me),2)+pow(programSizeTerm,2)+pow(SE,2));
//            cout<<"fitnessGP"<<fitness<<endl;

    }
    catch
        (...)
    {
        fitness=1000;
    }

    return fitness;
}

double FitnessFunction::ShannonComplexity(double *ODP)
{
    double fitness;
    int ndemb;
    ndemb=4;
    int  N= 24;//nfact
    double ODPprobs[N];
    double sumoODP=0;
    double ShanCompl=0;

    for(int i=0; i<N; i++)
    {
        sumoODP+=*(ODP+i);
//        cout<<*(ODP+i)<<"z"<<endl;
    }

    for(int i=0; i<N; i++)
    {
        ODPprobs[i]=*(ODP+i)/sumoODP;

        if(ODPprobs[i]>1e-30)
        {
            ShanCompl=ShanCompl+ODPprobs[i]*log(ODPprobs[i]);
        }

    }

    ShanCompl=((-1)*ShanCompl)/log(N);
    double diffCM=1-ShanCompl;
    fitness=(diffCM);
    delete ODP;
    return fitness;
}



double * FitnessFunction::ordinal_pattern_loop()//

{
    double * ifrec;
    double x[number_of_fitness_cases]; //residuals
    int ndemb,nx;
    ndemb=4;
    int ipa[ndemb][ndemb];
    int nd;
    int has_na;
    double xvec[ndemb];
    double epsilon=1e-10;
    //if (ISNA(x[0])) x[0]=0;
    //x[2]=10;
    //ifrec[0]=1;has
    ifrec=new double[24];//
    nx=number_of_fitness_cases;

    for(int i=0; i<number_of_fitness_cases; i++)
    {
        x[i]=(-*(predicted+i)+*(actual+i));
    }

    for(int i=0; i<24; i++)
    {
        ifrec[i]=0;
    }

    for(int nv=0; nv<nx-ndemb+1; nv++)
    {
        has_na=0;

        for(int ixvec=0; ixvec<ndemb; ixvec++)
        {
            if(std::isnan(x[nv+ixvec]))
            {
                has_na=1;
                break;
            }

            xvec[ixvec]=x[nv+ixvec]; /* Fill xvec */
        }

        if(has_na)
            continue;

        for(int i=0; i<ndemb; i++)
            for(int j=0; j<ndemb; j++)
                ipa[i][j]=0; /* Reset ipa matrix */

        for(int ilag=1; ilag<ndemb; ilag++)
        {
            for(int itime=ilag; itime<ndemb; itime++)
            {
                ipa[itime][ilag] = ipa[itime-1][ilag-1];

                if((xvec[itime] <= xvec[itime - ilag]) || (fabs(xvec[itime - ilag] - xvec[itime]) < epsilon))
                    ipa[itime][ilag] = ipa[itime][ilag] + 1;
            }
        }

        nd = ipa[ndemb-1][1];

        for(int ilag=2; ilag<ndemb; ilag++)
        {
            nd =(ilag+1) * nd + ipa[ndemb-1][ilag];
        }

        ifrec[nd] = ifrec[nd] + 1;
    }

//      for(int i=0; i<24; i++)
//    {
//        cout<<ifrec[i]<<"a"<<endl;
//    }

    return ifrec;
}


double FitnessFunction::ComputeFitness()
{
    double fitness;

    switch(type)
    {
    case NOH:
        return fitness=NumberOfHitsFitness();
        break;

    case rMSE:
        return fitness=relativeMeanSquaredError();
        break;

    case RSRE:
        return fitness=rootRelativeSquaredError();
        break;

    case rRAE:
        return fitness=RelativeAbsoluteError();
        break;

    case ME:
        return fitness=ModellingEfficiency();
        break;

    case sOS:
        return fitness=sumOfSquares();
        break;

    case AIC:
        return fitness=AkaikeIndex();
        break;

    case CM:
        return fitness=ShannonComplexity(ordinal_pattern_loop());
        break;
    case AICx:
        return fitness=xAkaikeIndex();
        break;
    case CEM:
        return fitness=ComplexityCorrectedEM();
        break;
    case MAE:
        return fitness=MeanAbsoluteError(); 
    default
            :
        break;
        return 0;
    }

    return 0;
}
