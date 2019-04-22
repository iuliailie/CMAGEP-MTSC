#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <stdarg.h>
#include <algorithm>
#include <float.h>
#include<sstream>
#include<list>
#include <array>
#include <Python.h>


#include <stdlib.h>

#include <map>
#include <random>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include<sys/types.h>
#include <unistd.h>
#include <cstdio>
#include <ctime>
#include <stdio.h>


using namespace std;
extern string functionSet;
extern string terminals;
extern string terminalsToUse;

extern string constants;

extern double **inputData;
extern double *targetData;
extern int number_of_fitness_cases;
extern int number_of_parameters;
extern double NBCalls;
extern double genesHeadLength;
extern double numberOfGenesInAChromosome;
extern vector<string> functionsVector;
extern vector<string> terminalsVector;
extern vector<string> terminalsToUseVector;
extern vector<string> constantsVector;
extern vector<string> allElementsVector;

std::string ReplaceAll(std::string str, const std::string &from, const std::string &to);

string addMathToPythonString(string str);

string simplifySolution(string str);

tuple<string, double, vector<double> *>
CallCma2(tuple<string, int, vector<double> *> functAndParam, string initialstr, int fTm,
         vector<double> *constanteOptimizateInainte);

tuple<string, int, vector<double> *> addParametersInPosition(string strings);

string makeUniformStringsForEditDistance(string strings);

string addParanthesesInPosition(string str);

string addParanthesesInPositionDivide(string str);

string getStringFromVector(vector<string> stringVec);


class FunctionStructure {
public :
    double (*pFunction)(double a, double b);

    int Arity;
};

enum nodeType {
    Terminal = 0, UnOperator = 1, BinOperator = 2
};

class Node {
public :
    nodeType type;
    string value;//value when replacing terminals
    string funct;//function
    double param = 1;
    Node *left;
    Node *right;

    Node(string val) {
        type = Terminal;
        value = val;
        funct = '\0';
        //   this->left=new Node('0');
        // this->right=new Node('0');
    }

    Node(string op, Node *left, Node *right) {
        type = BinOperator;
        value = "\0";
        funct = op;
        this->left = left;
        this->right = right;
    }

    ~Node() {
        //delete left,right;
    }

};


class Global {
public:
    Global();

    static double randomOne();

    static char charFromInt(int i);

    static char charFromFunction(int f);

//        enum Function
//        {Plus, Minus, Times, Divide, Sqrt, Exp, LessThan, GreaterThan, Abs, Log, Sin};
    static const int Variable;
    static const int Constant;


    map<char, FunctionStructure> *functionDictionary;

    void mapCharactersWithFunctions();

    Node expTree(string gene, map<char, FunctionStructure> dic);

    //  static bool isFunction(char a);
    ~Global() {
        functionDictionary->clear();
        delete functionDictionary;
    }


private:
};

static map<char, string> *functionFromChar();

static bool isFunction(string a) {
    if (functionSet.find(a.at(0)) != std::string::npos)//we check if the character is in the functions string
        return true;
    else
        return false;
}

static bool IsNumber(double x) {
    // but it's false if x is a NaN.
    return (x == x);
}


static bool IsFiniteNumber(double x) {
    return (x <= DBL_MAX && x >= -DBL_MAX);
}

static bool isConstant(string c) {
    string list = constants;
    int k = 0;
    k = list.find_first_of(c.at(0));

    if (k >= 0)
        return true;

    return false;
}

//static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
//{
//    std::stringstream ss(s);
//    std::string item;
//
//    while(std::getline(ss, item, delim))
//    {
//        elems.push_back(item);
//    }
//
//    return elems;
//}

static size_t split(const std::string &txt, std::vector<std::string> &strs, char ch) {
    size_t pos = txt.find(ch);
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while (pos != std::string::npos) {
        strs.push_back(txt.substr(initialPos, pos - initialPos));
        initialPos = pos + 1;

        pos = txt.find(ch, initialPos);
    }

    // Add the last one
    strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

    return strs.size();
}

static double computeTreeValue(Node *node, double *valueArray, map<char, FunctionStructure> *localMap) {
    if (node->type == Terminal) //if the node  of terminal type
    {
        if (isConstant(node->value)) {
            //cout<<node->value<<endl;
            double x = (double) (node->value.at(0) - '0');
            //cout<<x<<" const"<<endl;
            return x;
        }
        else
        {

            std::vector<string>::iterator iter=find (terminalsVector.begin(), terminalsVector.end(), node->value);
            int indexTerminal  = iter - terminalsVector.begin();
            double valueOfTerminal=valueArray[indexTerminal];
            if (IsNumber(valueOfTerminal) && IsFiniteNumber(valueOfTerminal))
            {
            //cout<<valueArray[terminals.find(node->value)]<<"term"<<endl;
            return valueOfTerminal;//we return the current value}
        }
        }
    } else if ((node->type == UnOperator)) //if it's an unary function
    {
        vector<double> vector1, vector2;
        vector1.push_back(computeTreeValue(node->right, valueArray,localMap));//we add the value from the right node of the current node, if it's not a terminal we apply the function recursively until we reach a terminal
        vector2.push_back(0);//the vector in which we'll store the result of the transform
        transform(vector1.begin(), vector1.end(), vector1.begin(), vector2.begin(),
                  localMap->at(node->funct.at(0)).pFunction);

        //we apply the function in the current node to the value in vector 1 and store it in vector2
        if (IsNumber(vector2.at(0)) && IsFiniteNumber(vector2.at(0))) {
            double val = vector2.at(0);
            //cout<<localMap->at(node->funct.at(0)).pFunction<<"\t"<<vector1.at(0)<<" dau \t"<<vector2.at(0)<<endl;
            vector1.clear();
            vector2.clear();
            //cout<<endl<<"unary"<<val<<endl;
            //cout<<val;
            return val;
        }

        vector1.clear();
        vector2.clear();
    } else
        //else it's binary
    {
        vector<double> vector1, vector2, vector3;
        vector1.push_back(computeTreeValue(node->left, valueArray, localMap));
        vector2.push_back(computeTreeValue(node->right, valueArray, localMap));
        vector3.push_back(0);
        transform(vector1.begin(), vector1.end(), vector2.begin(), vector3.begin(),
                  localMap->at(node->funct.at(0)).pFunction);

        if (IsNumber(vector3.at(0)) && IsFiniteNumber(vector3.at(0))) {
            double val = vector3.at(0);
            //cout<<node->funct.at(0)<<"\t"<<vector1.at(0)<<"si \t"<<vector2.at(0)<<"dau \t"<<vector3.at(0)<<endl;
            vector1.clear();
            vector2.clear();
            vector3.clear();
            //cout<<endl<<"binary"<<val<<endl;
            //cout<<localMap->at(node->funct.at(0)).pFunction<<"\t"<<vector1.at(0)<<"si \t"<<vector2.at(0)<<"dau \t"<<vector3.at(0)<<endl;
            //cout<<"val="<<val<<endl;
            return val;
        }

        vector1.clear();
        vector2.clear();
        vector3.clear();
    }

    return INFINITY;
}

static string printExpresie(Node *node) {
    map<char, string> *ffc = functionFromChar();
    string expression = "";
    stringstream ss;

    if (node->type == Terminal) {
        if (node->param > 1) {
            ss << node->param << "*";
        }

        ss << node->value;
        ss >> expression;
        expression = "(" + expression + ")";
        ffc->clear();
        delete (ffc);
        return expression;
    } else if (node->type == UnOperator) //if it's an unary function
    {
        ss << ffc->at(node->funct.at(0));
        ss >> expression;
        expression += "(" + printExpresie(node->right) + ")";
        ffc->clear();
        delete (ffc);
        return expression;
    } else
        //else it's binary
    {
        ss << ffc->at(node->funct.at(0));
        ss >> expression;
        expression = "(" + printExpresie(node->left) + expression + printExpresie(node->right) + ")";
        ffc->clear();
        delete (ffc);
        return expression;
    }
}

static double getARandom(double m, double n)//return a random double between m & n
{
    //This will generate a number from some arbitrary m to some arbitrary n:
    double r3 = m + (double) rand() / ((double) RAND_MAX / (n - m));
    return r3;
}


static int getARandom(int m, int n)//return a random integer between m & n
{
    int randNumber = rand() % (n - m + 1) + m;
    return randNumber;
}


static void simplifyNodeMap(Node *workNode) {
    bool one, two, three;

    if (workNode->type == Terminal) {
        return;
    } else {
        one = (workNode->left->type == Terminal);
        two = (workNode->right->type == Terminal);
        three = (workNode->left->value == workNode->right->value);

        switch (workNode->funct.at(0)) {
            case '/' :
                if (one && two && three) {
                    workNode->value = '1';
                    workNode->funct = '\0';
                    workNode->type = Terminal;
                    return;
                }

                if (!one)
                    simplifyNodeMap(workNode->left);

                if (!two)
                    simplifyNodeMap(workNode->right);

                break;

            case '-':
                if (one && two && three) {
                    workNode->value = '0';
                    workNode->type = Terminal;
                    workNode->funct = '\0';
                    return;
                }

                if (!one)
                    simplifyNodeMap(workNode->left);

                if (!two)
                    simplifyNodeMap(workNode->right);

                break;

            case '+':
                if (one && two && three) {
                    workNode->value = workNode->right->value;
                    workNode->type = Terminal;
                    workNode->funct = '\0';
                    workNode->param++;
                    return;
                }

                if (!one)
                    simplifyNodeMap(workNode->left);

                if (!two)
                    simplifyNodeMap(workNode->right);

                if (three) {
                    if (workNode->left->left->type == Terminal && workNode->left->right->type == Terminal)
                        if (workNode->left->left->value == workNode->right->left->value &&
                            workNode->left->right->value == workNode->right->right->value)
                            workNode = workNode->left;

                    workNode->param++;
                    return;
                }

                break;
        }
    }
}

static string getTimeString() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    string s = (string) buf;
    return s;
}


static map<char, string> *functionFromChar() {
    map<char, string> *fFC = new map<char, string>();
    fFC->insert(std::pair<char, string>('+', "+"));
    fFC->insert(std::pair<char, string>('-', "-"));
    fFC->insert(std::pair<char, string>('*', "*"));
    fFC->insert(std::pair<char, string>('/', "/"));
    fFC->insert(std::pair<char, string>('S', "sqrt"));
    fFC->insert(std::pair<char, string>('E', "exp"));
    fFC->insert(std::pair<char, string>('L', "log"));
    fFC->insert(std::pair<char, string>('s', "sin"));
    fFC->insert(std::pair<char, string>('C', "cos"));
    fFC->insert(std::pair<char, string>('^', "^"));
    fFC->insert(std::pair<char, string>('Q', "Q10"));
    return fFC;
}

static vector<string> split(const char *str, char c) {
    vector<string> result;

    do {
        const char *begin = str;

        while (*str != c && *str)
            str++;

        result.push_back(string(begin, str));
    } while (0 != *str++);

    return result;
}


#endif // GLOBAL_H
