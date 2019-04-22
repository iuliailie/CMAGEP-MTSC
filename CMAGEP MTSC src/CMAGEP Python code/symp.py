
from __future__ import division


import os

import subprocess
import sys
import random as ran

#
#sys.stdout = open('Results/outPy/outputPythonEU'+str(ran.randint(0,1000))+'.txt', 'w')#redirecting output to file
#sys.stderr = open('Results/outPy/erorLogEU'+str(ran.randint(0,1000))+'.txt', 'w')#redirecting output to file
'''py_function.py - Python source designed to '''
'''demonstrate the use of python embedding'''




def Teste(mathExpression):
    print "a intrat in simpy"


    mathExpression_pythonConverted=mathExpression.replace("^","**")
    print "mathExpression_pythonConverted="+str(mathExpression_pythonConverted)
    try:
#
        stringToRunPythonCode='python2.7 /Users/iulia/multiclassCMAGEP/timeout.py 10 \'python2.7 /Users/iulia/multiclassCMAGEP/interm.py -i "'+mathExpression_pythonConverted+'"\''#setting timeout for response of simplify

        simplifiedMathExpression= subprocess.check_output(stringToRunPythonCode, shell=True)

        mathExpression_NoGaps=str(simplifiedMathExpression[:-1]).replace("**","^")
        mathExpression_NoGaps=mathExpression_NoGaps.replace(" ","")
        print "simplifiedMathExpression="+str(mathExpression_NoGaps)
        con="conjugate"
        if con in mathExpression_NoGaps :
            print "a iesit din simpy conj"
            return '!'+mathExpression

        print "a iesit din simpy"
        print "\n"
        return mathExpression_NoGaps
    except:
        print "Unexpected error: ", sys.exc_info()[0]

        print "a iesit din simpy"
        print "\n"
        return '!'+mathExpression










