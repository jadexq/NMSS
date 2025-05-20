//
// Created by WANGJade on 3/3/2020.
//

#ifndef SIMU1_GENERATE_H
#define SIMU1_GENERATE_H

#include "Constant.h"

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include "rnorm.h"
#include "runif.h"
#include "runifVec.h"

using namespace std;
using namespace Eigen;

class Generate
{
private:
    // treatment x_i, len=n
    VectorXd x;
    // age z_i, len=n
    VectorXd z;
    // TBSS m_i, vv*n
    MatrixXd m;
    // cog y_i, len=n
    VectorXd y;
    // neighbourhood
    VectorXd nNbh; // number of nhb
    MatrixXd nbh; // index of nbh, =-1 if NA

public:
    // generate data
    // model:
    // M_i = Alpha * x_i + A * z_i + Epsilon_i
    // y_i = Beta * M_i + gamma * x_i + b * z_i + delta_i
    int generateX(Constant constant); // generate x_i, len=n
    int generateZ(Constant constant); // generate z_i, len=n
    int generateM(Constant constant); // generate M_i, vv*n
    int generateY(Constant constant); // generate y_i, len=n
    int generateNbh(Constant constant); // generate nNbh and nbh

    // neighbourhood


    // get data
    VectorXd getX(void);
    VectorXd getZ(void);
    MatrixXd getM(void);
    VectorXd getY(void);
    VectorXd getNnbh(void);
    MatrixXd getNbh(void);

    // write data
    int writeX(void);
    int writeZ(void);
    int writeM(void);
    int writeY(void);
    int writeNbh(void);

    Generate(Constant constant);
};


#endif //SIMU1_GENERATE_H
