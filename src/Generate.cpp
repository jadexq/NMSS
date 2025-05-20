//
// Created by WANGJade on 3/3/2020.
//

// steps:
// 0. Generate generate(Constant constant)
// 1. generate.generateX(constant), treatment
// 2. generate.generateZ(constant), age
// 3. generate.generateM(constant), TBSS
// 4. generate.generateY(constants), cog
// 5. generate.generateNbh(constants), nbh and nNbh

#include "Generate.h"

using namespace std;
using namespace Eigen;

Generate::Generate(Constant constant)
        : x(constant.n),
          z(constant.n),
          m(constant.vv, constant.n),
          y(constant.n),
          nbh(constant.vv, constant.mNnbh),
          nNbh(constant.vv)
{
    x = VectorXd::Constant(constant.n, 0);
    z = VectorXd::Constant(constant.n, 0);
    m = MatrixXd::Constant(constant.vv, constant.n, 0);
    y = VectorXd::Constant(constant.n, 0);
    nbh = MatrixXd::Constant(constant.vv, constant.mNnbh, -1);
    nNbh = VectorXd::Constant(constant.vv, 0);
}

// model:
// M_i = Alpha * x_i + A * z_i + Epsilon_i
// y_i = Beta * M_i + gamma * x_i + b * z_i + delta_i

// generate x
int Generate::generateX(Constant constant)
{
    // prob of treatment/x_i=1
    double p1=0.5;
    VectorXd runif_vec = runifVec(constant.n);
    x = (runif_vec.array() < p1).cast<double>();

    return 0;
}

// generate z (age)
int Generate::generateZ(Constant constant)
{
    // N[0,1]
    z = rnorm(constant.n);
    // binary
    //VectorXd runif_vec = runifVec(constant.n);
    //z = (runif_vec.array() < 0.5).cast<double>();

    return 0;
}

// generate m, vv*n
// M_i = Alpha * x_i + A * z_i + Epsilon_i
int Generate::generateM(Constant constant)
{
    // Alpha * x_i
    MatrixXd AA = constant.alpha.rowwise().replicate(constant.n); // vv*n
    //cout << "#row=" << constant.alpha.rows() << "  #col=" << constant.alpha.cols() << endl;
    //cout << "#row=" << AA.rows() << "  #col=" << AA.cols() << endl;
    AA = (AA.array().rowwise() * x.transpose().array()).matrix();
    // A * z_i
    MatrixXd AZ = constant.a1.rowwise().replicate(constant.n); // vv*n
    AZ = (AZ.array().rowwise() * z.transpose().array()).matrix();
    // Epsilon_i
    VectorXd rnorm_vec = rnorm(constant.n * constant.vv);
    Map<MatrixXd> E(rnorm_vec.data(), constant.vv, constant.n);
    E = E.array().colwise() * sqrt(constant.sigE.array());
    // sum up
    m = AA + AZ + E;

    return 0;
}

// generate y
// y_i = Beta * M_i + gamma * x_i + b * z_i + delta_i
int Generate::generateY(Constant constant)
{
    // Beta * M_i
    VectorXd BM = (constant.beta.transpose() * m).transpose();
    // gamma * x_i
    VectorXd GX = constant.gamma * x;
    // b * z_i
    VectorXd BZ = constant.b * z;
    // delta_i
    VectorXd D = rnorm(constant.n) * sqrt(constant.sigD);
    // sum up
    y = BM + GX + BZ + D;

    return 0;
}

// index of neighbourhood and number of neighbourhood
int Generate::generateNbh(Constant constant)
{
    // 2D
    if (constant.mNnbh == 4)
    {
        // nbh
        for (int v=0; v<constant.vv; v++)
        {
            VectorXd vNbh = VectorXd::Constant(constant.mNnbh, -1);
            VectorXd tmp(4);
            tmp << v-1, v+1, v-constant.v1, v+constant.v1;
            ArrayXd l0 = (tmp.array() >= 0).cast<double>().array();
            ArrayXd sv = (tmp.array() < (constant.vv-1)).cast<double>().array();
            ArrayXd effective = l0 * sv;
            vNbh = vNbh.array()*(1-effective) + tmp.array()*effective;
            // assign
            nbh.row(v) = vNbh;
        }
        // nNbh
        MatrixXd tmp = (nbh.array() > 0).cast<double>();
        nNbh = tmp.rowwise().sum();
    }

    // 3D
    if (constant.mNnbh == 6)
    {
        // nbh
        for (int v=0; v<constant.vv; v++)
        {
            VectorXd vNbh = VectorXd::Constant(6, -1);
            VectorXd tmp(6);
            tmp << v-1, v+1, v-constant.v1, v+constant.v1, v-constant.v1*constant.v2, v+constant.v1*constant.v2;
            ArrayXd l0 = (tmp.array() >= 0).cast<double>().array();
            ArrayXd sv = (tmp.array() < (constant.vv-1)).cast<double>().array();
            ArrayXd effective = l0 * sv;
            vNbh = vNbh.array()*(1-effective) + tmp.array()*effective;
            // assign
            nbh.row(v) = vNbh;
        }

        int i=0;
        int j=0;
        int k=0;
        for (k=0; k<constant.v3; k++)
        {
            for (j=0; j<constant.v1; j++)
            {
                // left face
                nbh(k*(constant.v1*constant.v2)+j,2) = -1;
                //cout << "idx=" << k*(v1*v2)+j+1 << endl;
                // right face
                nbh((k+1)*(constant.v1*constant.v2)-j-1,3) = -1;
                //cout << "idx=" << (k+1)*(v1*v2)-j << endl;
            }
            for (i=0; i<constant.v2;i++)
            {
                // back face
                nbh(k*(constant.v1*constant.v2)+i*constant.v1,0) = -1;
                //cout << "idx=" << k*(v1*v2)+i*v1+1 << endl;
                // front face
                nbh(k*(constant.v1*constant.v2)+(i+1)*constant.v1-1,1) = -1;
                //cout << "idx=" << k*(v1*v2)+(i+1)*v1 << endl;
            }
        }

        // nNbh
        MatrixXd tmp = (nbh.array() > 0).cast<double>();
        nNbh = tmp.rowwise().sum();
    }

    return 0;
}



// get data
VectorXd Generate::getX(void)
{
    return x;
}

VectorXd Generate::getZ(void)
{
    return z;
}

MatrixXd Generate::getM(void)
{
    return m;
}

VectorXd Generate::getY(void)
{
    return y;
}

MatrixXd Generate::getNbh(void)
{
    return nbh;
}

VectorXd Generate::getNnbh(void)
{
    return nNbh;
}

// write data to .txt files
int Generate::writeX(void)
{
    ofstream file;
    file.open("../data/x_gen.txt");
    file << x << endl;
    file.close();

    return 0;
}

int Generate::writeZ(void)
{
    ofstream file;
    file.open("../data/z_gen.txt");
    file << z << endl;
    file.close();

    return 0;
}

int Generate::writeM(void)
{
    ofstream file;
    file.open("../data/m_gen.txt");
    file << m.transpose() << endl;
    file.close();

    return 0;
}

int Generate::writeY(void)
{
    ofstream file;
    file.open("../data/y_gen.txt");
    file << y << endl;
    file.close();

    return 0;
}

int Generate::writeNbh(void)
{
    ofstream file;
    file.open("../data/nbh.txt");
    file << nbh << endl;
    file.close();

    ofstream file2;
    file2.open("../data/nNbh.txt");
    file2 << nNbh << endl;
    file2.close();

    return 0;
}
