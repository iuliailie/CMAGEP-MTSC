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
        es = cma.CMAEvolutionStrategy(optimizedParamList, 0.05, {'verb_disp': 0, 'maxiter': 50, 'maxfevals': 500, 'popsize':10})
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
