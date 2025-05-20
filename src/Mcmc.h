//
// Created by WANGJade on 4/3/2020.
//

#ifndef SIMU1_MCMC_H
#define SIMU1_MCMC_H

#include "Constant.h"
#include "Generate.h"
#include "Initial.h"

#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include "rnorm.h"
#include "runif.h"
#include "rgamma.h"
#include "strToVec.h"

using namespace std;
using namespace Eigen;

class Mcmc
{
private:
    // parameters to be updated
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
    VectorXd sparseIdx; // idx for unimportant connections/locations
    VectorXd aIA_vec; // vector version of aIA
    VectorXd aIB_vec; // vector version of bIA

    // data
    VectorXd x;
    VectorXd z;
    MatrixXd m;
    MatrixXd mb;
    VectorXd y;
    VectorXd c11;
    VectorXd c12;
    MatrixXd c2;
    MatrixXd nbhi;
    MatrixXd nbhj;
    //VectorXd nNbh;
    VectorXd validloc; // idx of valid locations that need an udpate
    int v = 0; // the indexing

    // for fast computation
    // fixed
    double sumX2; // sum_{i=1}^n x_i^2, scalar
    VectorXd sumXM; // sum_{i=1}^n m_i(v) x_i, vector, length = V
    VectorXd sumXMb;
    VectorXd sumM2; // sum_{i=1}^n m_i(v)^2, vector, length = V
    VectorXd sumMb2;
    VectorXd sumYM; // sum_{i=1}^n y_i m_i(v), vector, length = V
    VectorXd sumYMb;
    // update together with beta
    VectorXd sumBetaMb; // sum_{v=1}^V beta(v) m_i(v), vector, length = n
    // related to a1, a2, b
    double sumc112; // sum_{i=1}^n c1_1_i^2, scalar
    double sumc122; // sum_{i=1}^n c1_2_i^2, scalar
    VectorXd sumc22; // sum_{i=1}^n c2[j]_i^2, len=py
    VectorXd sumc11M;
    VectorXd sumc12M;
    MatrixXd sumc2Mb;
    double sumc11X; // sum_{i=1}^n x_i c1_1_i, scalar
    double sumc12X; // sum_{i=1}^n x_i c1_2_i, scalar

public:
    // read data and nbh from .txt
    int readData(Constant constant, Generate generate, int rep);

    // set initial value
    int setInitial(Constant constant, Initial initial);

    // for fast computation
    int fastComp(Constant constant);

    // mcmc updates
    // update alpha, zeta_alpha, a1, a2, and sigma^2_epsilon
    int updateA(Constant constant);
    // update beta and zeta_beta
    int updateB(Constant constant);
    // update gamma
    int updateGamma(Constant constant);
    // update b
    int updateb(Constant constant);
    // update sigma^2_delta
    int updateSigD(Constant constant);

    // get values
    VectorXd getAlpha(void);
    VectorXd getZetaA(void);
    VectorXd getA1(void);
    VectorXd getA2(void);
    VectorXd getBeta(void);
    VectorXd getZetaB(void);
    double getGamma(void);
    VectorXd getB(void);
    VectorXd getSigE(void);
    double getSigD(void);
    VectorXd getX(void);
    VectorXd getZ(void);
    MatrixXd getM(void);
    VectorXd getY(void);

    // constructor
    Mcmc(Constant constant);
};


#endif //SIMU1_MCMC_H
