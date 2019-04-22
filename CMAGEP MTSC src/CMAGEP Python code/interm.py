


from __future__ import division
import sys, getopt

sys.path.append('/Users/iuliailie/multiclassCMAGEP')
sys.path.append('/Users/iulia/anaconda3/lib/python3.7/site-packages')
#sys.stderr = open('Results/outPy/erorLog'+str(ran.randint(0,1000))+'.txt', 'w')#redirecting output to file
from sympy import *
import math
def main(argv):

   #print ('inter')
   mathExpression_pythonConverted=''
   try:
      opts, args = getopt.getopt(argv,"hi:")
   except getopt.GetoptError:
      sys.exit(2)
   for opt, arg in opts:
       mathExpression_pythonConverted = arg
   #print(local_PyMathExpression)
   ls='';
   for i in range(1,1000):
       ls=ls+'x_'+str(i)+" ";
   #ls='x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9';#list of GEP symbols, order counts!!
   l=ls.split(' ')
   for i in range(0,len(l)-1):
              l[i]=symbols(l[i])
   local_PyMathExpression=mathExpression_pythonConverted
   #print(local_PyMathExpression)

   for i in range(0,len(l)-1):
       local_PyMathExpression=local_PyMathExpression.replace("("+str(l[i])+")","(l["+str(i)+"])")

   #print(local_PyMathExpression)
   try:
       symplifiedMathExpression=N(eval(local_PyMathExpression),2)#function gets simplified
       con="conjugate"
       if con in str(symplifiedMathExpression) :#remove strange Sympy expression with conjugate
           symplifiedMathExpression=simplify(eval(local_PyMathExpression))
   except:
        symplifiedMathExpression=nsimplify(eval(local_PyMathExpression),rational=True)



   print(str(symplifiedMathExpression))

if __name__ == "__main__":
   main(sys.argv[1:])
