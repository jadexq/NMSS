//
// Created by WANGJade on 2/3/2020.
//

#include "Constant.h"
Constant::Constant() // use initializer list to construct the matrices
        : alpha(vv),
          //alphaM(v1*v2, v3),
          beta(vv),
          //betaM(v1*v2, v3),
          sigE(vv),
          tun1(nTun1),
          tun2(nTun2),
          tun3(nTun3),
          tun4(nTun4),
          tun5(nTun5),
          tun6(nTun6),
          tun7(nTun7),
          tun8(nTun8)
{
    // alpha
    //alphaAxialM = MatrixXd::Constant(v1, v2, 0);
    //alphaAxial = VectorXd::Constant(v1*v2, 0);
    //alphaM = MatrixXd::Constant(v1*v2, v3, 0);
    alpha = VectorXd::Constant(vv,0);
    //alphaAxialM.block(5,8,10,4) = MatrixXd::Constant(10,4,xxx1);
    //Map<VectorXd> tmpAA(alphaAxialM.data(), alphaAxialM.size());
    //alphaAxial = tmpAA;
    //for (int ia=8; ia<12; ia++)
    //{
    //    alphaM.col(ia) = alphaAxial;
    //}
    //Map<VectorXd> tmpA(alphaM.data(), alphaM.size());
    //alpha = tmpA;
    // a
    a1 = VectorXd::Constant(vv,0);
    a2 = VectorXd::Constant(vv,0);
    // beta
    //betaAxialM = MatrixXd::Constant(v1, v2, 0);
    //betaAxial = VectorXd::Constant(v1*v2, 0);
    //betaM = MatrixXd::Constant(v1*v2, v3, 0);
    beta = VectorXd::Constant(vv,0);
    //betaAxialM.block(8,5,4,10) = MatrixXd::Constant(4,10,xxx1);
    //Map<VectorXd> tmpBB(betaAxialM.data(), betaAxialM.size());
    //betaAxial = tmpBB;
    //for (int ia=8; ia<12; ia++)
    //{
    //        betaM.col(ia) = betaAxial;
    //}
    //Map<VectorXd> tmpB(betaM.data(), betaM.size());
    //beta = tmpB;
    b = VectorXd::Constant(py,0);
    // sigE
    sigE = VectorXd::Constant(vv, 0);
    // choice of tunings
    tun1 = VectorXd::LinSpaced(nTun1, 0.1, 1); // sig1Alpha
    tun2 = VectorXd::LinSpaced(nTun2, -4, 0); // aIA
    tun3 = VectorXd::LinSpaced(nTun3, 1, 2); // bIA
    tun4 = VectorXd::LinSpaced(nTun4, 0.55, 0.6); // sigDFix, for beta
    tun5 = VectorXd::LinSpaced(nTun5, 0.005, 0.0015); // sig1Beta
    tun6 = VectorXd::LinSpaced(nTun6, -1.5, 0); // aIB
    tun7 = VectorXd::LinSpaced(nTun7, 1, 2); // bIB
    tun8 = VectorXd::LinSpaced(nTun8, 0, 2); // cI
}

// change tunings
// alpha: 1-sig1Alpha 2-aIA 3-bIA
// beta:  4-sigDFix 5-sig1Beta 6-aIB 7-bIB
// AB: 8-cI
int Constant::changeTun(int rep)
{
    // for alpha ONLY, other tunings using default values
    if (0)
    {
        int idx1 = 0;
        int idx2 = 0;
        int idx3 = 0;
        idx3 = rep/(nTun1 * nTun2);
        idx2 = (rep-idx3*nTun1*nTun2)/nTun1;
        idx1 = rep % nTun1;
        sig1Alpha = tun1(idx1);
        aIA = tun2(idx2);
        bIA = tun3(idx3);
        cout << "idx1=" << idx1 << " " << "idx2=" << idx2 << " " << "idx3=" << idx3 << endl;
    }

    // for beta ONLY
    if (0)
    {
        int idx4 = 0;
        int idx5 = 0;
        int idx6 = 0;
        int idx7 = 0;
        idx7 = rep/(nTun4 * nTun5 * nTun6);
        idx6 = (rep-idx7*nTun4*nTun5*nTun6)/(nTun4 * nTun5);
        idx5 = (rep-idx7*nTun4*nTun5*nTun6-idx6*nTun4*nTun5)/nTun4;
        idx4 = rep % nTun4;
        sigDFix = tun4(idx4);
        sig1Beta = tun5(idx5);
        aIB = tun6(idx6);
        bIB = tun7(idx7);
        cout << "idx4=" << idx4 << " " << "idx5=" << idx5 << " " << "idx6=" << idx6 << " " << "idx7=" << idx7 << endl;
    }

    // for cI ONLY
    if (0)
    {
        int idx8 = rep;
        cI = tun8(idx8);
        cout << "idx8=" << idx8 << endl;
    }

    return 0;
}
