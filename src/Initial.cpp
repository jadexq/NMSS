//
// Created by WANGJade on 4/3/2020.
//

// steps:
// 0. Initial initial(constant)

#include "Initial.h"

using namespace std;
using namespace Eigen;

Initial::Initial(Constant constant)
        : alpha(constant.vv),
          a1(constant.vv),
          a2(constant.vv),
          zetaA(constant.vv),
          beta(constant.vv),
          zetaB(constant.vv),
          b(constant.py),
          sigE(constant.vv)
{
    alpha = VectorXd::Constant(constant.vv, 0); // set initial value
    a1 = VectorXd::Constant(constant.vv, 0);
    a2 = VectorXd::Constant(constant.vv, 0);
    zetaA = VectorXd::Constant(constant.vv, 0);
    beta = VectorXd::Constant(constant.vv, 0);
    // beta = constant.beta; // true value as initial
    zetaB = VectorXd::Constant(constant.vv, 0);
    gamma = 0.0;
    b = VectorXd::Constant(constant.py, 0.0);
    sigE = VectorXd::Constant(constant.vv, 1);
    sigD = constant.sigDFix;
}

// set true value (read from .txt) as initial value for beta
int Initial::setTrueBeta(void)
{
    // read beta true value
    // two tmp variables
    string strB; // store a line as str
    vector<double> vecB; // store a line as vec
    // open file
    ifstream fileB;
    fileB.open("../data/betaVec_true.txt");
    // assign
    int i = 0;
    while (getline(fileB, strB))
    {
        strToVec(strB, vecB);
        beta(i) = vecB[0];
        i++;
    } // end while i
    fileB.close();

    //cout << "betaInit2248=" << beta(2247) << endl;

    return 0;
}

// set true value (read from .txt) as initial value for alpha and zetaA
int Initial::setTrueAlpha(void)
{
    // read alpha true (vec) value
    // two tmp variables
    string strA; // store a line as str
    vector<double> vecA; // store a line as vec
    // open file
    ifstream fileA;
    fileA.open("../data/alphaT.txt");
    // assign
    int i = 0;
    while (getline(fileA, strA))
    {
        strToVec(strA, vecA);
        alpha(i) = vecA[0];
        i++;
    } // end while i
    fileA.close();

    //cout << "alphaIn[729]=" << alpha(728) << endl;
    //cout << "alphaIn[730]=" << alpha(729) << endl;

    zetaA = alpha/0.8;
    //cout << "zetaAIn[729]=" << zetaA(728) << endl;
    //cout << "zetaAIn[730]=" << zetaA(729) << endl;

    return 0;
}

