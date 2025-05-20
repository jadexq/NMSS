//
// Created by WANGJade on 4/3/2020.
//

#ifndef SIMU1_INITIAL_H
#define SIMU1_INITIAL_H

#include "Constant.h"

#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include "strToVec.h"

using namespace std;
using namespace Eigen;

class Initial
{
public:
    // parameters
    VectorXd alpha;
    VectorXd a1;
    VectorXd a2;
    VectorXd zetaA;
    VectorXd beta;
    VectorXd zetaB;
    double gamma;
    VectorXd b;
    VectorXd sigE;
    double sigD;

    // set true value (read from .txt) as initial value for beta
    int setTrueBeta(void);
    int setTrueAlpha(void);

    Initial(Constant constant);
};


#endif //SIMU1_INITIAL_H
