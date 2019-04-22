import re

def replaceConstants(strd):
    #print strd
    #strd=strd+'+P'
    tep1=re.findall(r'[!_\d\.]+(?:[eE][+-]?\d+)?',strd)
    iter=re.finditer(r'[!_\d\.]+(?:[eE][+-]?\d+)?',strd)
    indices = [m.start(0) for m in iter]
    tep=[]
    indT=[]
    for i in range(0, len(tep1)):
        if("_" not in tep1[i]):
            tep.append(tep1[i])
            indT.append(indices[i])

    try:
        tepp=[float(x) for x in tep]
        #tep=map(float, )
        l=len(tep)
        # print l
        for i in range(0,l):
                # print i
                sub="P";
                stringForReplacing=strd[0:indT[i]]+"P"+strd[indT[i]+len(tep[i]):]

                tempo=str(stringForReplacing)

            #print "in consts:  "+str(tempo)+"\n"
        return (str(tempo),tepp)
    except:
        return (str(strd),'')
    print(strd)
    return (str(strd),'')

import FitnessFunctionsForCma
import random
import cma

def CallFunctie(inputData, outputData, GEPmathExpression, nbParam, fitnessFunctionType, listOfTerminals,
                priorlyOptimisedConstants):

    print str(priorlyOptimisedConstants)

    PythonGEPmathExpression = GEPmathExpression.replace("^", "**")
    ParamReplaced_GEPmathExpression = PythonGEPmathExpression.replace("P", "paramToOptimize")
    print "cma function to optimize= " + str(ParamReplaced_GEPmathExpression)

    cmaESResults = [[0], ]
    cmaESResultsTemp = [[0], ]
    fitness = 1e50;
    conj = "conjugate"
    if conj in ParamReplaced_GEPmathExpression:  # make sure the simplified expression does not contain conjugate as it cannot be evaluated after
        print "a iesit din evs conj"
        return (cmaESResults, 2e50, ParamReplaced_GEPmathExpression, 0)

    optimizedParamList = []

    if len(priorlyOptimisedConstants) > 0:
        optimizedParamList = priorlyOptimisedConstants

    else:
        for jjin in range(0, len(priorlyOptimisedConstants)):
            optimizedParamList.append(random.uniform(0, 1))

    # print('a intrat la complexityMeasures \n')
    # print(str(optimizedParamList))
    if len(optimizedParamList) > 1:
        es = cma.CMAEvolutionStrategy(optimizedParamList, 0.05, {'verb_disp': 0, 'maxiter': 50, 'maxfevals': 500})
        # rezTemp=es.optimize(FitnessFunctionsForCma.ME,args=(X,Y,funcu,nbParam,terminalStr))

        cmaESResultsTemp = es.optimize(fitnessFunctionSwitch(fitnessFunctionType),
                                       args=(inputData, outputData, PythonGEPmathExpression, nbParam, listOfTerminals))
        # cmaESResultsTemp=cma.fmin(FitnessFunctionsForCma.ComplexityMeasures, optimizedParamList, 1,args=(inputData,outputData,PythonGEPmathExpression,nbParam,listOfTerminals))

        if cmaESResultsTemp[1] < fitness:
            fitness = cmaESResultsTemp[1]
            cmaESResults = cmaESResultsTemp

        paramToOptimize = [0]
        try:

            paramToOptimize = cmaESResults[0].tolist()
            me = FitnessFunctionsForCma.ComplexityMeasures
            call = me.called

            print str(call)

        except:
            return (optimizedParamList, 2e50, ParamReplaced_GEPmathExpression, 0)

        return (paramToOptimize, cmaESResults[1], ParamReplaced_GEPmathExpression, call)
    else:
        return (optimizedParamList, 2e50, ParamReplaced_GEPmathExpression, 0)

    return (cmaESResults, 2e50, ParamReplaced_GEPmathExpression, 0)


def fitnessFunctionSwitch(x):
    return {
        4: FitnessFunctionsForCma.ME,
        5: FitnessFunctionsForCma.AIC,
        6: FitnessFunctionsForCma.ComplexityMeasures,
        9: FitnessFunctionsForCma.CEM,
        10: FitnessFunctionsForCma.MAE,
    }[x]

# enum FitnessType {NOH=0,rMSE=1,RSRE=2,rRAE=3,ME=4,AIC=5,CM=6,AICx=7,sOS=8,CEM=9};


# def replaceConstantsOlder(strd):
#     #print strd
#     #strd=strd+'+P'
#     tep1=re.findall(r'[!^_\d\.]+(?:[eE][+-]?\d+)?',strd)
#     tep=[]
#     indT=[]
#     for i in range(0, len(tep1)):
#         if("_" not in tep1[i]):
#             tep.append(tep1[i])
#             indT.append(indices[i])
#
#     try:
#         tepp=[float(x) for x in tep]
#         #tep=map(float, )
#         l=len(tep)
#        # print l
#         if  l > 0:
#            sub="P";
#            stringForReplacing=strd.replace(tep[0],sub,1)
#            tempo=strd
#            for i in range(0,l):
#                # print i
#                 sub="P";
#                # stringForReplacing tep[i]
#                 stringForReplacing=strd.replace(tep[i],sub,1)
#                 #strd=stringForReplacing
#                 tempo=str(stringForReplacing)
#
#            #print "in consts:  "+str(tempo)+"\n"
#            return (str(tempo),tepp)
#     except:
#         return (str(strd),'')
#     print(strd)
#     return (str(strd),'')
#
#
# def replaceConstants2(strd):
#
#     final_list = []
#     for elem in strd.split():
#         try:
#             final_list.append(float(elem))
#         except ValueError:
#             pass
#     bla=strd.split()
#     print str(bla)
#     l=len(final_list)
#     print str(final_list)
#    # print l
#     if  l > 0:
#        sub="P";
#        tempString=strd.replace(final_list[0],sub,1)
#        StringToReturn=tempString
#        for i in range(0,l):
#            # print i
#             sub="P";
#            # print tep[i]
#             mghj=strd.replace(final_list[i],sub,1)
#             strd=mghj
#             StringToReturn=str(mghj)
#        return str(StringToReturn)
#     return str(strd)
