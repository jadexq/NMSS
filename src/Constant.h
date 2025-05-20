//
// Created by WANGJade on 2/3/2020.
//

#ifndef SIMU1_CONSTANT_H
#define SIMU1_CONSTANT_H

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Constant
{
public:
    // *** setting: 0-simulation, 1-real data (or read existing simulated data)
    int setting = 1;
    // 0-do not apply Cv, 1-appply Cv=sqrt(V)
    int applyCv = 1;

    // *** mcmc
    int nIter = 2000; // #mcmc iterations
    int nBurn = 1000; // #mcmc burn-in
    int nRep = 1; // #replication << must be specified properly!!!
    int nPrint = 100; // print "iter"

    // *** data
    int n = 271; // N sample size
    int v1 = 1; // #row
    int v2 = 1; // #col
    int v3 = 0;
    int vv = 71631; // V1*(V2-1)/2; V1=V2=p=379
    double cv = 1/sqrt(25431); // vv << the const for sum beta(v)*m(v)
    int mNnbh = 377; // max #nbh = p-2
    double pi = 3.14159;
    int py = 3; // #confounders in y's/2nd regression
    int nvl = 25431;// number of valid locations that need an updated, change-25432, T2-21610, T1-22130

    // *** parameter true values
    VectorXd alpha; // alpha
    double xxx1 = 0.3; // value of alpha
    VectorXd a1; // lam1 for c_1i_1
    VectorXd a2; // lam2 for c_1i_2
    double yyy1 = 0.2; // value of a/lam1
    VectorXd beta; // beta
    double xxx2 = 0.3; // value of beta
    VectorXd b; // (lam2, lam3) for (c_2i, c_31)
    double gamma = 0.2; // gamma
    VectorXd sigE; // sig^2_epsilon, = 0.1 (defined in .cpp)
    double sigD = 0.1; // sig^2_delta

    // *** hyper-parameters
    // * Ising (zeta)
    // p(zeta) ~ exp[ a*zeta + sum_l sum_l' b*I(zeta_l = zeta_l') ]
    // zeta_alpha
    double aIA = -2; // tuning!! << -9.33333
    double aIA0 = -2; // for unimportant connections
    double bIA = 0.00095; // tuning!! << 0.00167333
    // zeta_beta
    double aIB = -1; // tuning!! << -3
    double aIB0 = -1; // for unimportant connections
    double bIB = 0.005; // tuning!! << 0.012
    // association between zeta_alpha and zeta_beta
    double cI = 0.8; // tuning!! <<
    // * spike-and-slab
    double sig1Alpha = 0.0165; // slab tuning!! << 0.366667
    // sig0Alpha = Inf; // spike, point mass
    double sig1Beta = 0.0001; // slab tuning!! << 0.00173333
    double sig1BetaA = 0.0001; // slab tuning!! (adaptive sig1Beta) <<
    // sig0Beta = Inf // spike, point mass
    double sigDFix = 0.00069; // fix the var in y_i's regression << 0.00103333
    // * gamma, normal prior
    double mu0Gamma = 0;
    double sig0Gamma = 100;
    // * sig^2_epsilon, IG prior
    double a0sigE = 3;
    double b0sigE = 5;
    // * sig^2_delta, IG prior
    double a0sigD = 3;
    double b0sigD = 5;
    // for a and b (effect of age)
    double mu0a = 0;
    double sig0a = 100;
    double mu0b = 0;
    double sig0b = 100;

    // change tuning parameters:
    // alpha: 1-sig1Alpha 2-aIA 3-bIA
    // beta:  4-sigDFix 5-sig1Beta 6-aIB 7-bIB
    // AB: 8-cI
    // #choices
    int nTun = 5;
    int nTun1 = nTun; // alpha
    int nTun2 = nTun;
    int nTun3 = nTun;
    int nTun4 = nTun; // beta
    int nTun5 = nTun;
    int nTun6 = nTun;
    int nTun7 = nTun;
    int nTun8 = nTun; // cI
    // tuning choices, specified in .cpp
    VectorXd tun1;
    VectorXd tun2;
    VectorXd tun3;
    VectorXd tun4;
    VectorXd tun5;
    VectorXd tun6;
    VectorXd tun7;
    VectorXd tun8;

    // update tunings in each replication
    int changeTun(int rep);

    Constant();
};


#endif //SIMU1_CONSTANT_H
