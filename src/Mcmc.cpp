//
// Created by WANGJade on 4/3/2020.
//

/*
 how to use this class:
 steps:
 0. Mcmc:mcmc(Constant constant)
 1. readData(Constant constant, Generate generate)
 2. setInitial(Constant constant)
 3. fastComp(Constant constant)
 3. update
 */

#include "Mcmc.h"

using namespace std;
using namespace Eigen;

Mcmc::Mcmc(Constant constant)
        : x(constant.n),
          z(constant.n),
          m(constant.vv, constant.n),
          mb(constant.vv, constant.n),
          y(constant.n),
          c11(constant.n),
          c12(constant.n),
          c2(constant.n, constant.py),
          alpha(constant.vv),
          a1(constant.vv),
          a2(constant.vv),
          zetaA(constant.vv),
          beta(constant.vv),
          zetaB(constant.vv),
          sparseIdx(constant.vv),
          aIA_vec(constant.vv),
          aIB_vec(constant.vv),
          sigE(constant.vv),
          nbhi(constant.vv, constant.mNnbh),
          nbhj(constant.vv, constant.mNnbh),
          //nNbh(constant.vv),
          validloc(constant.nvl),
          sumXM(constant.vv),
          sumXMb(constant.vv),
          sumM2(constant.vv),
          sumMb2(constant.vv),
          sumYM(constant.vv),
          sumYMb(constant.vv),
          sumBetaMb(constant.n),
          //sumZM(constant.vv),
          //sumZMb(constant.vv)
          sumc11M(constant.vv),
          sumc12M(constant.vv),
          sumc2Mb(constant.vv, constant.py)
{
    x = VectorXd::Constant(constant.n, 0);
    z = VectorXd::Constant(constant.n, 0);
    m = MatrixXd::Constant(constant.vv, constant.n, 0);
    mb = MatrixXd::Constant(constant.vv, constant.n, 0);
    c11 = VectorXd::Constant(constant.n, 0);
    c12 = VectorXd::Constant(constant.n, 0);
    c2 = MatrixXd::Constant(constant.n, constant.py, 0);
    y = VectorXd::Constant(constant.n, 0);
    alpha = VectorXd::Constant(constant.vv, 0);
    a1 = VectorXd::Constant(constant.vv, 0);
    a2 = VectorXd::Constant(constant.vv, 0);
    zetaA = VectorXd::Constant(constant.vv, 0);
    beta = VectorXd::Constant(constant.vv, 0);
    zetaB = VectorXd::Constant(constant.vv, 0);
    sparseIdx = VectorXd::Constant(constant.vv, 0);
    aIA_vec = VectorXd::Constant(constant.vv, 0);
    aIB_vec = VectorXd::Constant(constant.vv, 0);
    b = VectorXd::Constant(constant.py, 0);
    sigE = VectorXd::Constant(constant.vv, 0);
    nbhi = MatrixXd::Constant(constant.vv, constant.mNnbh, 0);
    nbhj = MatrixXd::Constant(constant.vv, constant.mNnbh, 0);
    //nNbh = VectorXd::Constant(constant.vv, 0);
    validloc = VectorXd::Constant(constant.nvl, 0);
    sumXM = VectorXd::Constant(constant.vv, 0);
    sumXMb = VectorXd::Constant(constant.vv, 0);
    sumM2 = VectorXd::Constant(constant.vv, 0);
    sumMb2 = VectorXd::Constant(constant.vv, 0);
    sumYM = VectorXd::Constant(constant.vv, 0);
    sumYMb = VectorXd::Constant(constant.vv, 0);
    sumBetaMb = VectorXd::Constant(constant.n, 0);
    //sumZM = VectorXd::Constant(constant.vv, 0);
    //sumZMb = VectorXd::Constant(constant.vv, 0);
    sumc11M = VectorXd::Constant(constant.vv, 0);
    sumc12M = VectorXd::Constant(constant.vv, 0);
    sumc2Mb = MatrixXd::Constant(constant.vv, constant.py, 0);
}

// read data (x, m, y) and nbh, nNbh from .txt
int Mcmc::readData(Constant constant, Generate generate, int rep)
{
    // *** read from "Generate" class
    // x = generate.getX();
    // z = generate.getZ();
    // m = generate.getM();
    // mb = m;
    // y = generate.getY();
    // nbh = generate.getNbh();
    // nNbh = generate.getNnbh();

    // *** read from file
    if (constant.setting == 0)
    {
        cout << "setting=0" << endl;
    }

    if (constant.setting == 1)
    {
        //cout << "read data..." << endl;
        // * x
        // two tmp variables
        string strX; // store a line as str
        vector<double> vecX; // store a line as vec
        // open file
        ifstream fileX;
        fileX.open("../data/it_std.txt");
        // assign
        int i = 0;
        while (getline(fileX, strX))
        {
            //cout << "line=" << i << endl;
            strToVec(strX, vecX);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                x(i-constant.n*rep) = vecX[0];
            }
            i++;
        } // end while i
        fileX.close();
        // cout << "x[1:5]=" << x.head(5) << endl;

        // * z
        // two tmp variables
        // string strZ; // store a line as str
        // vector<double> vecZ; // store a line as vec
        // open file
        // ifstream fileZ;
        // fileZ.open("../data/z_gen.txt");
        // assign
        // i = 0;
        // while (getline(fileZ, strZ))
        // {
        //    strToVec(strZ, vecZ);
        //    if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
        //     {
        //         z(i-constant.n*rep) = vecZ[0];
        //     }
        //     i++;
        // } // end while i
        // fileZ.close();
        //cout << "zRep[1:5]=" << z.head(5) << endl;

        // * c11 corresponding to para a1/lam1
        // two tmp variables
        string strC11; // store a line as str
        vector<double> vecC11; // store a line as vec
        // open file
        ifstream fileC11;
        fileC11.open("../data/age.txt");
        // assign
        i = 0;
        while (getline(fileC11, strC11))
        {
            strToVec(strC11, vecC11);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                c11(i-constant.n*rep) = vecC11[0];
            }
            i++;
        } // end while i
        fileC11.close();
        // cout << "c11[1:5]=" << c11.head(5) << endl;

        // * c12 corresponding to para a2
        // two tmp variables
        string strC12; // store a line as str
        vector<double> vecC12; // store a line as vec
        // open file
        ifstream fileC12;
        fileC12.open("../data/sex.txt");
        // assign
        i = 0;
        while (getline(fileC12, strC12))
        {
            strToVec(strC12, vecC12);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                c12(i-constant.n*rep) = vecC11[0];
            }
            i++;
        } // end while i
        fileC12.close();
        // cout << "c12[1:5]=" << c12.head(5) << endl;

        // * c2 corresponding to para b/(lam2, lam3, ...)
        // two tmp variables
        string strC2; // store a line as str
        vector<double> vecC2; // store a line as vec
        // open file
        ifstream fileC2;
        fileC2.open("../data/c3.txt");
        // assign
        i = 0;
        while (getline(fileC2, strC2))
        {
            strToVec(strC2, vecC2);
            for (int j=0; j<constant.py; j++)
            {
                c2(i,j) = vecC2[j];
            }

            i++;
        } // end while i
        fileC2.close();
        // cout << "c2[1:5,:]=" << c2.block(0,0,5,2) << endl;

        // * m/mb
        // two tmp variables
        string strM; // store a line as str
        vector<double> vecM; // store a line as vec
        // open file
        ifstream fileM;
        fileM.open("../data/mtC.txt");
        // assign
        i = 0;
        while (getline(fileM, strM))
        {
            // convert str to vec
            strToVec(strM, vecM);
            // assign row_i in file to col_i in matrix xi
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                for (int j=0; j<constant.vv; j++)
                {
                    m(j,i-constant.n*rep) = vecM[j];
                } // end j
            } // end if i
            i++;
        } // end while i
        fileM.close();
        // currently mb = m, No Cv
        mb = m;
        if (constant.applyCv == 1)
        {
            mb = m * constant.cv;
            cout << "Cv is applied to m" << endl;
        }
        // cout << "mb[1:5,1:5]=" << mb.block(0,0,5,5) << endl;

        // * y
        // two tmp variables
        string strY; // store a line as strr
        vector<double> vecY; // store a line as vec
        // open file
        ifstream fileY;
        fileY.open("../data/y.txt");
        // assign
        i = 0;
        while (getline(fileY, strY))
        {
            strToVec(strY, vecY);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                y(i-constant.n*rep) = vecY[0];
            }
            i++;
        } // end while i
        fileY.close();
        // cout << "y[1:5]=" << y.head(5) << endl;

        // * nbhi
        // two tmp variables
        string strNi; // store a line as str
        vector<double> vecNi; // store a line as vec
        // open file
        ifstream fileNi;
        fileNi.open("../data/nbhi.txt");
        // assign
        i = 0;
        while (getline(fileNi, strNi))
        {
            // convert str to vec
            strToVec(strNi, vecNi);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constant.mNnbh; j++)
            {
                nbhi(i,j) = vecNi[j];
            } // end j
            i++;
        } // end while i
        fileNi.close();
        // cout << "nbhi[1:5,1:5]=" << nbhi.block(0,0,5,5) << endl;

        // * nbhj
        // two tmp variables
        string strNj; // store a line as str
        vector<double> vecNj; // store a line as vec
        // open file
        ifstream fileNj;
        fileNj.open("../data/nbhj.txt");
        // assign
        i = 0;
        while (getline(fileNj, strNj))
        {
            // convert str to vec
            strToVec(strNj, vecNj);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constant.mNnbh; j++)
            {
                nbhj(i,j) = vecNj[j];
            } // end j
            i++;
        } // end while i
        fileNj.close();
        // cout << "nbhj[1:5,1:5]=" << nbhj.block(0,0,5,5) << endl;

        // * nNbh
        // two tmp variables
        //string strNn; // store a line as str
        //vector<double> vecNn; // store a line as vec
        // open file
        //ifstream fileNn;
        //fileNn.open("../data/nNbh.txt");
        // assign
        //i = 0;
        //while (getline(fileNn, strNn))
        //{
        //    strToVec(strNn, vecNn);
        //    nNbh(i) = vecNn[0];
        //    i++;
        //} // end while i
        //fileNn.close();

        // * sprarseIdx, for unimportant connections/locations
        // two tmp variables
        string strS; // store a line as str
        vector<double> vecS; // store a line as vec
        // open file
        ifstream fileS;
        fileS.open("../data/idx_dis80_T2.txt");
        // assign
        i = 0;
        while (getline(fileS, strS))
        {
            strToVec(strS, vecS);
            sparseIdx(i) = vecS[0];
            i++;
        } // end while i
        fileS.close();
        // cout << "idx_dis80[1:10]=" << sparseIdx.head(10) << endl;
        // aIA_vec and aIB_vec
        aIA_vec = constant.aIA0*sparseIdx.array() + constant.aIA*(1-sparseIdx.array());
        // cout << "aIA_vec[1:10]=" << aIA_vec.head(10) << endl;
        aIB_vec = constant.aIB0*sparseIdx.array() + constant.aIB*(1-sparseIdx.array());
        // cout << "aIB_vec[1:10]=" << aIB_vec.head(10) << endl;

        // * valid locations need update
        // two tmp variables
        string strV; // store a line as strr
        vector<double> vecV; // store a line as vec
        // open file
        ifstream fileV;
        fileV.open("../data/validlocC.txt");
        // assign
        i = 0;
        while (getline(fileV, strV))
        {
            strToVec(strV, vecV);
            validloc(i) = vecV[0];
            i++;
        } // end while i
        fileV.close();
        // cout << "y[1:5]=" << y.head(5) << endl;
    }

    return 0;
}

// set initial value for mcmc updates
int Mcmc::setInitial(Constant constant, Initial initial)
{
    alpha = initial.alpha;
    a1 = initial.a1;
    a2 = initial.a2;
    zetaA = initial.zetaA;
    beta = initial.beta;
    zetaB = initial.zetaB;
    gamma = initial.gamma;
    b = initial.b;
    sigE = initial.sigE;
    sigD = initial.sigD;

    //cout << "MCMC_ini_beta2147=" << beta(2146) << endl;

    return 0;
}

// for fast computation
int Mcmc::fastComp(Constant constant)
{
    // fixed
    sumX2 = pow(x.array(),2).sum(); // sum_{i=1}^n x_i^2, scalar
    sumXM = ((m.array().rowwise()) * (x.transpose().array())).rowwise().sum(); // sum_{i=1}^n m_i(v) x_i, vector, length = V
    sumXMb = ((mb.array().rowwise()) * (x.transpose().array())).rowwise().sum();
    sumM2 = pow(m.array(), 2).rowwise().sum(); // sum_{i=1}^n m_i(v)^2, vector, length = V
    sumMb2 = pow(mb.array(), 2).rowwise().sum();
    sumYM = ((m.array().rowwise()) * (y.transpose().array())).rowwise().sum();  // sum_{i=1}^n y_i m_i(v), vector, length = V
    sumYMb = ((mb.array().rowwise()) * (y.transpose().array())).rowwise().sum();
    // need to be updated along with beta
    sumBetaMb = beta.transpose() * mb; // sum_{v=1}^V beta(v) m_i(v), vector, length = n
    // related to a1, a2, b
    sumc112 = pow(c11.array(),2).sum(); // sum_{i=1}^n c11_i^2, scalar
    sumc122 = pow(c12.array(),2).sum(); // sum_{i=1}^n c12_i^2, scalar
    sumc22 = pow(c2.array(),2).colwise().sum(); // sum_{i=1}^n c2[j]_i^2, len=py
    sumc11M = ((m.array().rowwise()) * (c11.transpose().array())).rowwise().sum(); // sum_{i=1}^n m_i(v) c11_i, vector, length = V
    sumc12M = ((m.array().rowwise()) * (c12.transpose().array())).rowwise().sum(); // sum_{i=1}^n m_i(v) c12_i, vector, length = V
    //sumZMb = ((mb.array().rowwise()) * (z.transpose().array())).rowwise().sum();
    sumc2Mb = mb * c2; // dim=vv*py
    sumc11X = (x.array() * c11.array()).sum(); // sum_{i=1}^n x_i c11_i, scalar
    sumc12X = (x.array() * c12.array()).sum(); // sum_{i=1}^n x_i c12_i, scalar

    return 0;
}

// get values
VectorXd Mcmc::getAlpha(void)
{
    return alpha;
}

VectorXd Mcmc::getA1(void)
{
    return a1;
}

VectorXd Mcmc::getA2(void)
{
    return a2;
}

VectorXd Mcmc::getZetaA(void)
{
    return zetaA;
}

VectorXd Mcmc::getBeta(void)
{
    return beta;
}

VectorXd Mcmc::getZetaB(void)
{
    return zetaB;
}

double Mcmc::getGamma(void)
{
    return gamma;
}

VectorXd Mcmc::getB(void)
{
    return b;
}

VectorXd Mcmc::getSigE(void)
{
    return sigE;
}

double Mcmc::getSigD(void)
{
    return sigD;
}

VectorXd Mcmc::getX(void)
{
    return x;
}

VectorXd Mcmc::getZ(void)
{
    return z;
}

MatrixXd Mcmc::getM(void)
{
    return m;
}

VectorXd Mcmc::getY(void)
{
    return y;
}

// update alpha, zeta_alpha, a1, a2, and sigma^2_epsilon
int Mcmc::updateA(Constant constant)
{
    for (int i=0; i<constant.nvl; i++)
    {
        v = validloc(i);
        // *** update alpha and zeta_alpha
        // generate alpha(v) proposal from alpha(v)|zeta^alpha = 1, all others
        double sig2 = 1/1/(sumX2/sigE(v) + 1/constant.sig1Alpha);
        double barAlphaV = 0;
        double zb01 = 0; // sum((zetaAlpha[nbh[v,]] == 0)-(zetaAlpha[nbh[v,]] == 1), na.rm = T)
        int npv = 0; // number of nhb
        for (int i=0; i<constant.mNnbh; i++)
        {
            if ((nbhi(v,i) >= 0)&(nbhj(v,i) >= 0))
            {
                barAlphaV = barAlphaV + alpha(nbhi(v,i)) + alpha(nbhj(v,i));
                npv = npv + 1;
                //zb01 = zb01 + (zetaA(nbh(v,i)) == 0) - (zetaA(nbh(v,i)) == 1);
            }
            // Simple Ising prior
            zb01 = zb01 + ((zetaA(nbhi(v,i)) == 0)+(zetaA(nbhj(v,i)) == 0)) - ((zetaA(nbhi(v,i)) == 1)+(zetaA(nbhj(v,i)) == 1));
        }
        barAlphaV = 0.5*barAlphaV/(npv+0.001);
        double mu = sumXM(v)/sigE(v) - a1(v)*sumc11X/sigE(v) - a2(v)*sumc12X/sigE(v) + barAlphaV/constant.sig1Alpha;
        mu = sig2 * mu;
        double alphaVP = mu + rnorm(1)(0)*sqrt(sig2);
        //cout << "alphaVP=" << alphaVP << endl;
        // compute g_v
        // where g_v = P(alpha(v)=0, zeta(v)=0|all) / P(alpha(v)=alpha(v)P, zeta(v)=1|all)
        double gv = 0;
        gv = gv -0.5/sigE(v)*(2*alphaVP*sumXM(v)-a1(v)*sumc11X-a2(v)*sumc12X-pow(alphaVP,2)*sumX2) + 0.5/constant.sig1Alpha*pow(alphaVP - barAlphaV,2);
        //gv = gv - constant.aIA + 2*constant.bIA*zb01 + constant.cI*((zetaB(v) == 0) - (zetaB(v) == 1));
        //gv = gv - constant.aIA + 2*constant.bIA*zb01 - constant.cI*zetaB(v); // [aIA as scalar]
        gv = gv - aIA_vec(v) + 2*constant.bIA*zb01 - 2*constant.cI*zetaB(v); // [aIA vector]
        //cout << "gv_logll=" << (-0.5/sigE(v)*(2*alphaVP*sumXM(v)-a(v)*sumc11X-pow(alphaVP,2)*sumX2)) << endl;
        //cout << "gv_log_GMRF=" << (0.5/constant.sig1Alpha*pow(alphaVP - barAlphaV,2)) << endl;
        //cout << "gv_log_BI=" << (- aIA_vec(v) + 2*constant.bIA*zb01 - constant.cI*zetaB(v)) << endl;
        gv = sqrt(2*constant.pi*constant.sig1Alpha)*exp(gv);
        // prob of alpha(v)=alphaVP, zeta(v)=1| all = 1/(g_v + 1)
        //cout << "gv=" << gv << endl;
        double prob = 1/(gv+1);
        //cout << "probV=" << prob << endl;
        // ac
        double rand = runif();
        if (rand < prob)
        {
            alpha(v) = alphaVP;
            zetaA(v) = 1;
        }
        else
        {
            alpha(v) = 0;
            zetaA(v) = 0;
        }

        // fix alpha=0 for debugging
        //alpha(v) = 0; 
        //zetaA(v) = 0;

        // update a1(v), a1(v)
        double barAV1 = 0; // GMRF prior mean
        double barAV2 = 0;
        for (int i=0; i<constant.mNnbh; i++)
        {
            barAV1 = barAV1 + a1(nbhi(v,i)) + a1(nbhj(v,i));
            barAV2 = barAV2 + a2(nbhi(v,i)) + a2(nbhj(v,i));
        }
        barAV1 = barAV1/constant.mNnbh;
        barAV2 = barAV2/constant.mNnbh;
        double sigStar1 = 1/(sumc112/sigE(v) + 1/constant.sig0a);
        double sigStar2 = 1/(sumc122/sigE(v) + 1/constant.sig0a);
        //double muStar = (sumc11M(v) - alpha(v)*sumc11X)/sigE(v) + constant.mu0a/constant.sig0a; // simple normal prior
        double muStar1 = (sumc11M(v) - alpha(v)*sumc11X - a2(v)*sumc11X)/sigE(v) + barAV1/constant.sig0a; // update a1(v)
        muStar1 = sigStar1 * muStar1;
        a1(v) = muStar1 + rnorm(1)(0)*sqrt(sigStar1); 
        double muStar2 = (sumc12M(v) - alpha(v)*sumc12X - a1(v)*sumc12X)/sigE(v) + barAV2/constant.sig0a; // update a2(v)
        muStar2 = sigStar2 * muStar2;
        a2(v) = muStar2 + rnorm(1)(0)*sqrt(sigStar2); 

        // fix a = 0 for debugging
        //a1(v) = 0; 
        //a1(v) = 0; 

        // *** update sigE(v)
        double aS = constant.a0sigE + 0.5*constant.n;
        double bS = constant.b0sigE + 0.5*pow(m.row(v).array()-alpha(v)*x.transpose().array()-a1(v)*c11.transpose().array()-a2(v)*c12.transpose().array(), 2).sum();
        sigE(v) = 1/rgamma(1, aS, 1/bS)(0); // NOTE rgamma parameterization
    }

    return 0;
}


// update beta and zeta_beta
int Mcmc::updateB(Constant constant)
{
    for (int i=0; i<constant.nvl; i++)
    {
        v = validloc(i);
        // sum_{v1!=v} beta(v1) m_i(v1)
        if (beta(v) != 0)
        {
            sumBetaMb = sumBetaMb - beta(v) * mb.row(v).transpose();
        }
        //cout << "v=" << v << endl;
        //cout << "nbhi dim=" << nbhi.rows() << "*" << nbhi.cols() << endl;
        //cout << "mnbh=" << constant.mNnbh << endl;

        // adjust sig1Beta with zetaA(v)
        //double sig1Beta = constant.sig1Beta * (1-zetaA(v)) + constant.sig1BetaA*zetaA(v);
        double sig1Beta = constant.sig1Beta; // no adjustment
        // generate beta(v) proposal from beta(v)|zeta^beta = 1, all others
        double sig2 = 1/(sumMb2(v)/sigD + 1/sig1Beta);
        double barBetaV = 0;
        double zb01 = 0; // sum((zetaBeta[nbh[v,]] == 0)-(zetaBeta[nbh[v,]] == 1), na.rm = T)
        int npv = 0; // number of nhb
        for (int i=0; i<constant.mNnbh; i++)
        {
            if ((nbhi(v,i)>=0)&(nbhj(v,i)>=0))
            {
                barBetaV = barBetaV + beta(nbhi(v,i)) + beta(nbhj(v,i));
                npv = npv + 1;
                // with zeta iteration only
                //zb01 = zb01 + ((zetaB(nbhi(v,i)) == 0)|(zetaB(nbhj(v,i)) == 0)) - ((zetaB(nbhi(v,i)) == 1)&(zetaB(nbhj(v,i)) == 1));
            }
            // Simple Ising prior
            zb01 = zb01 + ((zetaB(nbhi(v,i)) == 0)+(zetaB(nbhj(v,i)) == 0)) - ((zetaB(nbhi(v,i)) == 1)+(zetaB(nbhj(v,i)) == 1));
        }
        barBetaV = 0.5*barBetaV/(npv+0.001);
        //cout << "barBetaV=" << barBetaV << endl;
        double mu = 0;
        mu = (sumYMb(v) - gamma*sumXMb(v) - sumc2Mb.row(v)*b)/sigD - (mb.row(v).transpose().array()*sumBetaMb.array()).sum() - barBetaV/sig1Beta;
        mu = sig2 * mu;
        //cout << "sumYMb(v)=" << sumYMb(v) << endl;
        //cout << "mu=" << mu << endl;
        double betaVP = mu + rnorm(1)(0)*sqrt(sig2);
        // compute g_v
        // where g_v = P(beta(v)=0, zeta(v)=0|all) / P(beta(v)=beta(v)P, zeta(v)=1|all)
        double gv = 0;
        gv = -0.5/sigD * (2*betaVP*(sumYMb(v) - (mb.row(v).transpose().array()*sumBetaMb.array()).sum() - gamma*sumXMb(v) - sumc2Mb.row(v)*b) - pow(betaVP,2)*sumMb2(v));
        gv = gv + 0.5/sig1Beta * pow(betaVP - barBetaV,2);
        //gv = gv - constant.aIB + 2*constant.bIB*zb01 + constant.cI*((zetaA(v)==0) - (zetaA(v)==1)); [wrong]
        //gv = gv - constant.aIB + 2*constant.bIB*zb01 - constant.cI*zetaA(v); [aIB scalar]
        gv = gv - aIB_vec(v) + 2*constant.bIB*zb01 - 2*constant.cI*zetaA(v); // [aIB vector]
        //gv = gv - aIB_vec(v) + 2*(constant.bIB + constant.bIA*zetaA(v))*zb01 - constant.cI*zetaA(v); // [aIB vector and adjust bIB = bIB + bIA*zetaA(v)]
        //cout << "v=" << v << " bIA*zetaA(v)=" << constant.bIA*zetaA(v) << endl;
        //cout << "gv_logll=" << (-0.5/sigD * (2*betaVP*(sumYMb(v) - (mb.row(v).transpose().array()*sumBetaMb.array()).sum() - gamma*sumXMb(v) - sumc21Mb.row(v)*b) - pow(betaVP,2)*sumMb2(v))) << endl;
        //cout << "gv_log_GMRF=" << (0.5/sig1Beta * pow(betaVP - barBetaV,2)) << endl;
        //cout << "gv_log_BI=" << (- aIB_vec(v) + 2*constant.bIB*zb01 - constant.cI*zetaA(v)) << endl;
        gv = sqrt(2*constant.pi*sig1Beta) * exp(gv);
        //cout << "gv=" << gv << endl;
        // prob of beta(v)=betaVP, zeta(v)=1| all = 1/(g_v + 1)
        double prob = 1/(gv+1);
        //cout << "probV=" << prob << endl;
        // ac
        double rand = runif();
        if (rand < prob)
        {
            beta(v) = betaVP;
            zetaB(v) = 1;
            // update sumBetaM = sum_{v=1}^V beta(v) m_i(v)
            sumBetaMb = sumBetaMb + beta(v) * mb.row(v).transpose();
        }
        else
        {
            beta(v) = 0;
            zetaB(v) = 0;
        }
    }

    return 0;
}

// update gamma
int Mcmc::updateGamma(Constant constant)
{
    double sigStar = 1/(sumX2/sigD + 1/constant.sig0Gamma);
    double muStar = ((y - sumBetaMb - c2*b).array()*x.array()).sum()/sigD + constant.mu0Gamma/constant.sig0Gamma;
    muStar = sigStar * muStar;
    gamma = muStar + rnorm(1)(0)*sqrt(sigStar);

    return 0;
}

// update b
int Mcmc::updateb(Constant constant)
{
    for (int j=0; j<constant.py; j++)
    {
        double sigStar = 1/(sumc22(j)/sigD + 1/constant.sig0b);
        double muStar = ((y - sumBetaMb - gamma*x - c2*b + c2.col(j)*b(j)).array()*c2.col(j).array()).sum()/sigD + constant.mu0b/constant.sig0b;
        muStar = sigStar * muStar;
        b(j) = muStar + rnorm(1)(0)*sqrt(sigStar);
        //cout << "muStar=" << muStar << " sigStar=" << sigStar << " b=" << b << endl;
        //cout << "sumZ2=" << sumZ2 << endl;
    }

    // fix the effect of bertain b at 0
    //b(2) = 0;
    //b(6) = 0;

    return 0;
}

// update sigD
int Mcmc::updateSigD(Constant constant)
{
    double aS = constant.a0sigD + 0.5*constant.n;
    double bS = constant.b0sigD + pow((y - sumBetaMb - gamma*x - b*c2).array(),2).sum();
    sigD = 1/rgamma(1, aS, 1/bS)(0); // NOTE rgamma parameterization

    return 0;
}

