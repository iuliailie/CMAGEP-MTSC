#include "Global.h"

#include <stdio.h>
#include <cassert>
#include <algorithm>
#include <functional>
#include <cmath>
#include <sstream>


using namespace std;


string functionSet;
string terminals;
string constants;
string terminalsToUse;
double **inputData;
double *targetData;
int number_of_fitness_cases;
int number_of_parameters;
double NBCalls;
double genesHeadLength;
double numberOfGenesInAChromosome;
vector<string> functionsVector;
vector<string> terminalsVector;
vector<string> terminalsToUseVector;
vector<string> constantsVector;
vector<string> allElementsVector;

Global::Global() {
    //ctor
}


//definitions of functions we will use for evaluation of functions
double Plus(double a, double b) {
    return (a + b);
}

double Minus(double a, double b) {
    return (a - b);
}

double Divide(double a, double b) {
    return (a / b);
}

double Times(double a, double b) {
    return (a * b);
}

double Sqrt(double a, double b) {
    return sqrt(a);
}

double Exp(double a, double b) {
    return exp(a);
}

double Log(double a, double b) {
    return log(a);
}

double Sin(double a, double b) {
    return sin(a);
}

double Cos(double a, double b) {
    return cos(a);
}

double PoW(double a, double b) {
    double x;
    x = pow(a, b);
    return x;
}

double QTenHighF(double a,
                 double b)//if this UDF(user defined function) is used, it will take automaticaly as tau value the last element from the terminals vector
{
    double lnQTen = a;
    // double lnR=lnQTen*terminals.at(terminals.size()-1);
    double Resp = pow(a, (terminals.at(0) - 10) / 15);
    return Resp;
}


void Global::mapCharactersWithFunctions() {
    functionDictionary = new map<char, FunctionStructure>();
    FunctionStructure plusF, minusF, timesF, divideF, sqrtF, expF, logF, sinF, cosF, powF, QTenHighFFunct;
    plusF.Arity = 2;
    plusF.pFunction = Plus;
    functionDictionary->insert(std::pair<char, FunctionStructure>('+', plusF));
    minusF.Arity = 2;
    minusF.pFunction = Minus;
    functionDictionary->insert(std::pair<char, FunctionStructure>('-', minusF));
    timesF.Arity = 2;
    timesF.pFunction = Times;
    functionDictionary->insert(std::pair<char, FunctionStructure>('*', timesF));
    divideF.Arity = 2;
    divideF.pFunction = Divide;
    functionDictionary->insert(std::pair<char, FunctionStructure>('/', divideF));
    sqrtF.Arity = 1;
    sqrtF.pFunction = Sqrt;
    functionDictionary->insert(std::pair<char, FunctionStructure>('S', sqrtF));
    expF.Arity = 1;
    expF.pFunction = Exp;
    functionDictionary->insert(std::pair<char, FunctionStructure>('E', expF));
    logF.Arity = 1;
    logF.pFunction = Log;
    functionDictionary->insert(std::pair<char, FunctionStructure>('L', logF));
    sinF.Arity = 1;
    sinF.pFunction = Sin;
    functionDictionary->insert(std::pair<char, FunctionStructure>('s', sinF));
    cosF.Arity = 1;
    cosF.pFunction = Cos;
    functionDictionary->insert(std::pair<char, FunctionStructure>('C', cosF));
    powF.Arity = 2;
    powF.pFunction = PoW;
    functionDictionary->insert(std::pair<char, FunctionStructure>('^', powF));
    QTenHighFFunct.Arity = 1;
    QTenHighFFunct.pFunction = QTenHighF;
    functionDictionary->insert(std::pair<char, FunctionStructure>('Q', QTenHighFFunct));
}


string getStringFromVector(vector<string> stringVec) {
    string strToReturn;

    for (int i = 0; i < stringVec.size(); i++) {
        if (strToReturn.size() == 0)
            strToReturn = strToReturn + stringVec.at(i);
        else strToReturn = strToReturn + ' ' + stringVec.at(i);


    }
    return strToReturn;
}

Node Global::expTree(string gene, map<char, FunctionStructure> dic) {
    Node *node = new Node("!", NULL, NULL);
    return *node;
}


string addMathToPythonString(string str) {
    vector<string> *functions = new vector<string>();
    functions->push_back("sqrt");
    functions->push_back("exp");
    functions->push_back("log");
    functions->push_back("sin");
    functions->push_back("cos");

    functionsVector;
//    vector<char> *terminalVector = new vector<char>();
//
//    for (int i = 0; i < terminals.size(); i++) {
//        terminalVector->push_back(terminals.at(i));
//    }

    //  terminals->push_back("P");
    string sub; // str is string to search, sub is the substring to search for
    vector<size_t> positions; // holds all the positions that sub occurs within str

    for (int i = 0; i < functions->size(); i++) {
        sub = functions->at(i);
        size_t pos = str.find(sub, 0);
        int dif = 0;

        while (pos != string::npos) {
            positions.push_back(pos);
            pos = str.find(sub, pos + 1);
            dif++;
        }

        for (int k = positions.size() - dif; k < positions.size(); k++) {
            int j = positions.at(k);
            int found = 0;
            int t = j + (functions->at(i)).length();

            while (str.at(t) != ')') {
                if (str.at(t) == '(') {
                    found++;
                }

                t++;
            }

            while (found > 0) {
                if (str.at(t) == ')') {
                    found--;
                }

                t++;
            }

            int foundTerminal = 0;

            for (int s = 0; s < terminalsVector.size(); s++) {
                string substring = str.substr(j + (functions->at(i)).length(), t - (j + (functions->at(i).length())));
                int posT = substring.find(terminalsVector.at(s));

                if (posT >= 0)
                    foundTerminal++;
            }

            if (foundTerminal == 0) {
                str.replace(positions.at(k), (functions->at(i)).length(), "math." + functions->at(i));

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 5;
            }
        }
    }

    int a = 0;
    functions->clear();
    return str;
}


string addMathToPythonStringCMA(string str) {
    vector<string> *functions = new vector<string>();
    functions->push_back("sqrt");
    functions->push_back("exp");
    functions->push_back("log");
    functions->push_back("sin");
    functions->push_back("cos");
    vector<char> *terminalVector = new vector<char>();

    for (int i = 0; i < terminals.size(); i++) {
        terminalVector->push_back(terminals.at(i));
    }

    //  terminals->push_back("P");
    string sub; // str is string to search, sub is the substring to search for
    vector<size_t> positions; // holds all the positions that sub occurs within str

    for (int i = 0; i < functions->size(); i++) {
        sub = functions->at(i);
        size_t pos = str.find(sub, 0);
        int dif = 0;

        while (pos != string::npos) {
            positions.push_back(pos);
            pos = str.find(sub, pos + 1);
            dif++;
        }

        for (int k = positions.size() - dif; k < positions.size(); k++) {
            int j = positions.at(k);
            int found = 0;
            int t = j + (functions->at(i)).length();

            while (str.at(t) != ')') {
                if (str.at(t) == '(') {
                    found++;
                }

                t++;
            }

            while (found > 0) {
                if (str.at(t) == ')') {
                    found--;
                }

                t++;
            }

            int foundTerminal = 0;

            for (int s = 0; s < terminalVector->size(); s++) {
                string substring = str.substr(j + (functions->at(i)).length(), t - (j + (functions->at(i).length())));
                int posT = substring.find(terminalVector->at(s));

                if (posT >= 0)
                    foundTerminal++;
            }

            str.replace(positions.at(k), (functions->at(i)).length(), "np." + functions->at(i));

            for (int index = k; index < positions.size(); index++)
                positions.at(index) += 3;
        }
    }

    int a = 0;
    functions->clear();
    return str;
}

string simplifySolution(string str) {
    //str="log(((2^t+3*3)))*log((sin(cos(3))+2))*log((sin(cos(3))+2))";
    PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *pArgs, *prez;
//     Build the name object
    //pName = PyString_FromString("symp");
    // Load the module object
    pModule = PyImport_ImportModule("symp");
    if (pModule == nullptr) {
        PyErr_Print();
        std::exit(1);
    }
//     pDict is a borrowed reference
    pDict = PyModule_GetDict(pModule);
    // pFunc is also a borrowed reference
    pFunc = PyDict_GetItemString(pDict, "Teste");
    string result;
    string asd;

    if (PyCallable_Check(pFunc)) {
        pArgs = PyTuple_New(1);
        // asd=addMathToPythonString(str);
        str = addMathToPythonString(str);
        pValue = PyString_FromString(str.c_str());
        int k = PyTuple_SetItem(pArgs, 0, pValue);
        prez = PyObject_CallObject(pFunc, pArgs);

        if (prez != NULL) {
            result = PyString_AsString(prez);
//            Py_CLEAR(pValue);
//            Py_CLEAR(pArgs);
//            Py_CLEAR(pFunc);
//            Py_CLEAR(pDict);
//             Py_CLEAR(prez);
        } else {
            result = str;
            PyErr_Print();
//           Py_CLEAR(pValue);
//           Py_CLEAR(pArgs);
        }
    } else {
        PyErr_Print();
    }

//     Clean up
    Py_CLEAR(pModule);
    //Py_CLEAR(pName);
    // Finish the Python Interpreter
    return result;
}

PyObject *makelistFromVector(vector<double> *vect) {
    PyObject *l = PyList_New(vect->size());

    for (size_t i = 0; i != vect->size(); ++i) {
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(vect->at(i)));
    }

    return l;
}

PyObject *makelist(double *array, size_t size) {
    PyObject *l = PyList_New(size);

    for (size_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(array[i]));
    }

    return l;
}

PyObject *inputDataAsPythonList() {
    PyObject *l = PyList_New(number_of_fitness_cases);

    for (int i = 0; i < number_of_fitness_cases; i++) {
        PyObject *interList = PyList_New(number_of_parameters);

        for (int j = 0; j < number_of_parameters; j++) {
            PyList_SetItem(interList, j, PyFloat_FromDouble(inputData[i][j]));
        }

//      interList= makelist(*(inputData+number_of_parameters*i),number_of_fitness_cases);
        PyList_SetItem(l, i, interList);
    }

    return l;
}

tuple<string, double, vector<double> *>
CallCma2(tuple<string, int, vector<double> *> functAndParam, string initialstr, int fT,
         vector<double> *constanteOptimizateInainte) {
    vector<double> *constante;
    double fitness = 2e50;
    int size = number_of_fitness_cases;
    string strin = get<0>(functAndParam);

    int npar = get<1>(functAndParam);

    vector<double> *parametriDinGEP = get<2>(functAndParam);

    if (npar == 1) {
        strin = strin + "+P[1]";
        npar++;

    }

    cout << endl;
    string interStr = strin;

    strin = addMathToPythonStringCMA(strin);
    PyObject *pModule1 = PyImport_ImportModule("Consts");
    if (pModule1 == nullptr) {
        PyErr_Print();
        //std::exit(1);
    }
    if (pModule1 != NULL) {
        PyObject *pDict1 = PyModule_GetDict(pModule1);
        PyObject *myfunc1 = PyDict_GetItemString(pDict1, "CallFunctie");
        PyObject *mylist1 = PyList_New(size);
        PyObject *mylist2 = PyList_New(size);

        mylist1 = inputDataAsPythonList();
        mylist2 = makelist(targetData, size);


        PyObject *listaFostiParametri;
        if (constanteOptimizateInainte->size() > 0) {
            listaFostiParametri = PyList_New(constanteOptimizateInainte->size());
            listaFostiParametri = makelistFromVector(constanteOptimizateInainte);


        } else {

            listaFostiParametri = PyList_New(parametriDinGEP->size());
            listaFostiParametri = makelistFromVector(parametriDinGEP);


        }

        PyObject *arglist1 = PyTuple_New(7);
        PyTuple_SetItem(arglist1, 0, mylist1);
        PyTuple_SetItem(arglist1, 1, mylist2);
        PyTuple_SetItem(arglist1, 2, PyString_FromString(strin.c_str()));
        PyTuple_SetItem(arglist1, 3, PyInt_FromLong(npar));
        PyTuple_SetItem(arglist1, 4, PyInt_FromLong(fT));
        PyTuple_SetItem(arglist1, 5, PyString_FromString(terminals.c_str()));
        PyTuple_SetItem(arglist1, 6, listaFostiParametri);


        PyObject *result;
        string rez = "";


        constante = new vector<double>();
        if (npar < 100) {
            result = PyObject_CallObject(myfunc1, arglist1);
        } else {
            Py_CLEAR(pModule1);
            return make_tuple(initialstr, fitness, constante);
        }

        if (result != NULL && PyFloat_AsDouble(PyTuple_GetItem(result, 1)) < 1e50) {
            cout << "optim fitness= " << PyFloat_AsDouble(PyTuple_GetItem(result, 1)) << endl;
            PyObject *result1 = NULL;
            PyObject *inter = PyTuple_GetItem(result, 2);
            string sInter = interStr;
            NBCalls = PyFloat_AsDouble(PyTuple_GetItem(result, 3));
            sInter = addMathToPythonString(sInter);
            PyObject *tup = PyTuple_New(3);
            PyTuple_SetItem(tup, 0, PyTuple_GetItem(result, 0));
            PyTuple_SetItem(tup, 1, PyString_FromString(sInter.c_str()));
            PyTuple_SetItem(tup, 2, PyString_FromString(terminals.c_str()));
            PyObject *pnume = PyString_FromString("EvalFunc");
            PyObject *pModual = PyImport_Import(pnume);
            PyObject *pDictas = PyModule_GetDict(pModual);
            PyObject *mynewfunc = PyDict_GetItemString(pDictas, "Evs");
            result1 = PyObject_CallObject(mynewfunc, tup);
            if (pModual == nullptr) {
                PyErr_Print();
                std::exit(1);
            }
            if (result1 != NULL) {
                rez = PyString_AsString(PyTuple_GetItem(result1, 0));
                int k = (int) PyList_Size(PyTuple_GetItem(result1, 1));
                constante = new vector<double>();
                PyObject *listConst = PyTuple_GetItem(result1, 1);

                for (int i = 0; i < k; i++) {
                    constante->push_back(PyFloat_AsDouble(PyList_GetItem(listConst, i)));
                }

                fitness = PyFloat_AsDouble(PyTuple_GetItem(result, 1));
            } else {
                PyErr_Print();
                Py_CLEAR(pModual);
                Py_CLEAR(pnume);
                Py_CLEAR(pModule1);
                return make_tuple(initialstr, fitness, constante);
            }

            Py_CLEAR(pModual);
            Py_CLEAR(pnume);
        } else {
            PyErr_Print();
            Py_CLEAR(pModule1);
            return make_tuple(initialstr, fitness, constante);
        }

        Py_CLEAR(pModule1);
        return make_tuple(rez, fitness, constante);

    }
    return make_tuple(initialstr, fitness, constante);
}

bool isFloat(string myString) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail();
}

tuple<string, int, vector<double> *> addParametersInPosition(string strings) {

    strings = ReplaceAll(strings, "^0.5", "^(1/2)");
    //strings=strings+"+0";
    string str = strings;
    vector<double> *constante;
    vector<string> *functions = new vector<string>();
    functions->push_back("sqrt");
    functions->push_back("exp");
    functions->push_back("log");
    functions->push_back("sin");
    functions->push_back("cos");

    //string term=terminals;
    vector<char> *terminalVector = new vector<char>();

    for (int i = 0; i < terminals.size(); i++) {
        terminalVector->push_back(terminals.at(i));
    }

    string sub; // str is string to search, sub is the substring to search for
    vector<size_t> positions; // holds all the positions that sub occurs within str
    int nbOfParameters = 0;
    PyObject *pModule2 = PyImport_ImportModule("Consts");
    PyObject *pDict2 = PyModule_GetDict(pModule2);
    PyObject *myfunc2 = PyDict_GetItemString(pDict2, "replaceConstants");
    PyObject *arglist2 = PyTuple_New(1);
    PyTuple_SetItem(arglist2, 0, PyString_FromString(str.c_str()));

    if (PyCallable_Check(myfunc2)) {
        PyObject *result2 = PyObject_CallObject(myfunc2, arglist2);

        if (result2 != NULL) {
            constante = new vector<double>();
            str = PyString_AsString(PyTuple_GetItem(result2, 0));
            int k = (int) PyList_Size(PyTuple_GetItem(result2, 1));
            constante = new vector<double>();
            PyObject *listConst = PyTuple_GetItem(result2, 1);

            for (int i = 0; i < k; i++) {
                constante->push_back(PyFloat_AsDouble(PyList_GetItem(listConst, i)));
            }
            Py_CLEAR(result2);
        } else {
            PyErr_Print();
        }
    }

    Py_CLEAR(pModule2);

//add params for terminals
    for (int i = terminalsVector.size()-1; i >=0 ; i--) {
        sub = terminalsVector.at(i);
        size_t pos = str.find(sub, 0);
        int dif = 0;

        while (pos != string::npos) {
            positions.push_back(pos);
            pos = str.find(sub, pos + 1);
            dif++;
        }

        for (int k = positions.size() - dif; k < positions.size(); k++) {
            string s;
            s+=(terminalsVector.at(i));

            if (positions.at(k) < 2) {
                str.replace(positions.at(k), terminalsVector.at(i).size(), "(H*" + s + ")");

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 4;

                nbOfParameters++;
            } else if (str.substr((positions.at(k) - 2), 2) != "H*") {
                str.replace(positions.at(k), terminalsVector.at(i).size(), "(H*" + s + ")");

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 4;

                nbOfParameters++;
            }

            s.clear();
        }
    }


//
//    while

//add params in position for functions
    for (int i = 0; i < functions->size(); i++) {
        sub = functions->at(i);
        size_t pos = str.find(sub, 0);
        int dif = 0;

        while (pos != string::npos) {
            positions.push_back(pos);
            pos = str.find(sub, pos + 1);
            dif++;
        }

        for (int k = positions.size() - dif; k < positions.size(); k++) {
            if (positions.at(k) < 2) {
                str.replace(positions.at(k), (functions->at(i)).length(), "F*" + functions->at(i));

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 2;

                nbOfParameters++;
            } else if (str.at(positions.at(k) - 2) != 'F') {
                str.replace(positions.at(k), (functions->at(i)).length(), "F*" + functions->at(i));

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 2;

                nbOfParameters++;
            }
        }
    }

    string stringDeIntors = addParanthesesInPosition(str);
//stringDeIntors=addParanthesesInPositionDivide(stringDeIntors);
    str = stringDeIntors;

    vector<double> *vectorDeParametri = new vector<double>();

    int c1 = std::count(str.begin(), str.end(), 'P');
    int c2 = std::count(str.begin(), str.end(), 'F');

    int indexp = 0, indexf = 0;
    while (indexp < str.size()) {
        if (str.at(indexp) == 'P') {
            vectorDeParametri->push_back(constante->at(indexf));
            indexf++;
        }
        if (str.at(indexp) == 'F' || str.at(indexp) == 'H') {
            vectorDeParametri->push_back(1);

        }
        indexp++;
    }
    replace(str.begin(), str.end(), 'F', 'P');
    replace(str.begin(), str.end(), 'H', 'P');


    positions.clear();
    sub = "P";
    size_t pos = str.find(sub, 0);
    int dif = 0;

    while (pos != string::npos) {
        positions.push_back(pos);
        pos = str.find(sub, pos + 1);
        dif++;
    }

    for (int k = positions.size() - dif; k < positions.size(); k++) {
        string s = "P[" + to_string(k) + "]";
        string helper = str;
        str.replace(positions.at(k), 1, s);

        if (str == "P[0]")
            str = helper;

        for (int index = k; index < positions.size(); index++)
            positions.at(index) += s.length() - 1;
    }
    if (str == "P")
        str = "P[0]";

//        string stringDeIntors=addParanthesesInPosition(str);

    nbOfParameters = positions.size();//  return str;
    functions->clear();
    return make_tuple(str, nbOfParameters, vectorDeParametri);
}


std::string ReplaceAll(std::string str, const std::string &from, const std::string &to) {
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

int countSubstring(const std::string &str, const std::string &sub) {
    if (sub.length() == 0) return 0;
    int count = 0;
    for (size_t offset = str.find(sub); offset != std::string::npos;
         offset = str.find(sub, offset + sub.length())) {
        ++count;
    }
    return count;
}

string addParanthesesInPosition(string str) {
    int cntInit = countSubstring(str, "^(F");
    str = ReplaceAll(str, "^F", "^(F");

    str = ReplaceAll(str, "E", "exp(1)");
    bool incep = false;

    int lenStr = str.size();
    int r = 0, l = 0, index = 0, interIndex = 0;

    int cnt = countSubstring(str, "^(F") - cntInit;

    while (index < lenStr) {
        string sub = str.substr(0, index);
        string temp = str.substr(index, 3);
        if (incep == false && str.substr(index, 3) == "^(F") {
            incep = true;
            interIndex = index + 1;
        }

        if (incep && str.at(index) == '(') {
            r++;

        }
        if (incep && str.at(index) == ')') {
            l++;


        }

        if (cnt > 0 && (r - l) == cnt && l > cnt - 1 && (str.at(index - 1) == ')' || str.at(index - 1) == ']')) {

            str.insert(index, ")");
            index = interIndex + 1;
            lenStr++;
            l = 0;
            r = 0;
            incep = false;
            cnt--;
        }

        index++;


    }


    return str;

}


string addParanthesesInPositionDivide(string str) {
    int cntInit = countSubstring(str, "/(-");
    str = ReplaceAll(str, "/F", "/(F");

    str = ReplaceAll(str, "E", "exp(1)");
    bool incep = false;

    int lenStr = str.size();
    int r = 0, l = 0, index = 0, interIndex = 0;

    int cnt = countSubstring(str, "/(F");
    cnt = 1;
    while (index < lenStr) {
        string sub = str.substr(0, index);
        string temp = str.substr(index, 3);
        if (incep == false && str.substr(index, 3) == "/(F") {
            incep = true;
            interIndex = index + 1;
        }

        if (incep && str.at(index) == '(') {
            r++;

        }
        if (incep && str.at(index) == ')') {
            l++;


        }

        if ((r - l) == 1 && l > cnt - 1 && (str.at(index - 1) == ')' || str.at(index - 1) == ']')) {

            str.insert(index, ")");
            index = interIndex + 1;
            lenStr++;
            l = 0;
            r = 0;
            incep = false;
            //cnt++;
        }

        index++;



//(H*A)*(H*B)*F*exp(-(H*A)+(H*B)-P)/(H*D)-(H*A)+P*(H*B)+(H*B)/(F*cos(P*F*exp(F*exp((H*A)))/(H*B)))+(H*B)/(-(H*A)+F*log((H*A)+(H*B)))
    }


    return str;

    //"-(H*A)+(H*A)/(F*log((H*A)))+F*sin((H*A))+(H*B)*((H*B)+F*cos(F*exp((H*B)/(F*cos((H*A))))))/(H*A)+(P*(H*B)-(H*D)-(H*D)/(-(H*B)+(H*D))+F*sin((H*B)))/(H*A)"
//H*A)+(H*A)/(F*log((H*A)))+F*sin((H*A))+(H*B)*((H*B)+F*cos(F*exp((H*B)/(F*cos((H*A))))))/(H*A)+(P*(H*B)-(H*D)-(H*D)/(-(H*B)+(H*D))+F*sin((H*B)))/(H*A)
//"(H*A)*(H*B)*F*exp(-(H*A)+(H*B)-P)/(H*D)-(H*A)+P*(H*B)+(H*B)/(F*cos(P*F*exp(F*exp((H*A)))/(H*B)))+(H*B)/(-(H*A)+F*log((H*A)+(H*B)))
}


string makeUniformStringsForEditDistance(string strings) {
    string str = strings;
    vector<string> *functions = new vector<string>();
    functions->push_back("sqrt");
    functions->push_back("exp");
    functions->push_back("log");
    functions->push_back("sin");
    functions->push_back("cos");

    //string term=terminals;
    vector<char> *terminalVector = new vector<char>();

    for (int i = 0; i < terminals.size(); i++) {
        terminalVector->push_back(terminals.at(i));
    }

    string sub; // str is string to search, sub is the substring to search for
    vector<size_t> positions; // holds all the positions that sub occurs within str
    int nbOfParameters = 0;
    PyObject *pModule2 = PyImport_ImportModule("Consts");
    PyObject *pDict2 = PyModule_GetDict(pModule2);
    PyObject *myfunc2 = PyDict_GetItemString(pDict2, "replaceConstants");
    PyObject *arglist2 = PyTuple_New(1);
    PyTuple_SetItem(arglist2, 0, PyString_FromString(str.c_str()));

    if (PyCallable_Check(myfunc2)) {
        PyObject *result2 = PyObject_CallObject(myfunc2, arglist2);

        if (result2 != NULL) {
            str = PyString_AsString(PyTuple_GetItem(result2, 0));
            Py_CLEAR(result2);
        } else {
            PyErr_Print();
        }
    }

    Py_CLEAR(pModule2);

    for (int i = 0; i < functions->size(); i++) {
        sub = functions->at(i);
        size_t pos = str.find(sub, 0);
        int dif = 0;

        while (pos != string::npos) {
            positions.push_back(pos);
            pos = str.find(sub, pos + 1);
            dif++;
        }

        for (int k = positions.size() - dif; k < positions.size(); k++) {
            if (positions.at(k) < 2) {
                str.replace(positions.at(k), (functions->at(i)).length(), "P*" + functions->at(i));

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 2;

                nbOfParameters++;
            } else if (str.at(positions.at(k) - 2) != 'P') {
                str.replace(positions.at(k), (functions->at(i)).length(), "P*" + functions->at(i));

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 2;

                nbOfParameters++;
            }
        }
    }

    for (int i = 0; i < terminalVector->size(); i++) {
        sub = terminalVector->at(i);
        size_t pos = str.find(sub, 0);
        int dif = 0;

        while (pos != string::npos) {
            positions.push_back(pos);
            pos = str.find(sub, pos + 1);
            dif++;
        }

        for (int k = positions.size() - dif; k < positions.size(); k++) {
            string s;
            s.push_back(terminalVector->at(i));

            if (positions.at(k) < 2) {
                str.replace(positions.at(k), 1, "(P*" + s + ")");

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 4;

                nbOfParameters++;
            } else if (str.substr((positions.at(k) - 2), 2) != "P*") {
                str.replace(positions.at(k), 1, "(P*" + s + ")");

                for (int index = k; index < positions.size(); index++)
                    positions.at(index) += 4;

                nbOfParameters++;
            }

            s.clear();
        }
    }

    positions.clear();

    functions->clear();
    return str;
}







