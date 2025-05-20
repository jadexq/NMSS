#include <iostream>
#include <ctime>
#include <Eigen/Dense>

#include "Constant.h"
#include "Generate.h"
#include "Initial.h"
#include "Mcmc.h"

using namespace std;
using namespace Eigen;

int main()
{
    // clock
    clock_t start_t, end_t;
    double total_t;
    // starting time
    start_t = clock();

    // *** instantiation "Constant"
    int rep=0;
    Constant constant;
    cout << "setting: " << constant.setting << endl;

    // *** store rep values, and the quality for cross validation
    MatrixXd rep_alpha = MatrixXd::Zero(constant.nRep, constant.vv);
    MatrixXd rep_zetaA = MatrixXd::Zero(constant.nRep, constant.vv);
    MatrixXd rep_a1 = MatrixXd::Zero(constant.nRep, constant.vv);
    MatrixXd rep_a2 = MatrixXd::Zero(constant.nRep, constant.vv);
    MatrixXd rep_beta = MatrixXd::Zero(constant.nRep, constant.vv);
    MatrixXd rep_zetaB = MatrixXd::Zero(constant.nRep, constant.vv);
    VectorXd rep_gamma = VectorXd::Zero(constant.nRep);
    MatrixXd rep_b = MatrixXd::Zero(constant.nRep, constant.py);
    MatrixXd rep_sigE = MatrixXd::Zero(constant.nRep, constant.vv);
    // VectorXd rep_sigD = VectorXd::Zero(constant.nRep);

    // *** replication
    for (rep=0; rep<constant.nRep; rep++)
    {
        cout << "replication " << rep+1 << endl;
        ofstream file_rep;
        file_rep.open("../output/rep.txt");
        file_rep << rep+1 << endl;
        file_rep.close();

        // *** change tuning choice
        //constant.changeTun(rep);
        cout << "sig1Alpha=" << constant.sig1Alpha << endl;
        cout << "aIA=" << constant.aIA << endl;
        cout << "bIA=" << constant.bIA << endl;
        cout << "sigDFix=" << constant.sigDFix << endl;
        cout << "sig1Beta=" << constant.sig1Beta << endl;
        cout << "aIB=" << constant.aIB << endl;
        cout << "bIB=" << constant.bIB << endl;
        cout << "cI=" << constant.cI << endl;

        // *** store posterior samples
        MatrixXd samples_alpha =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        MatrixXd samples_zetaA =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        MatrixXd samples_a1 =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        MatrixXd samples_a2 =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        MatrixXd samples_beta =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        MatrixXd samples_zetaB =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        VectorXd samples_gamma =  VectorXd::Constant(constant.nIter-constant.nBurn, 0);
        MatrixXd samples_b =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.py, 0);
        MatrixXd samples_sigE =  MatrixXd::Constant(constant.nIter-constant.nBurn, constant.vv, 0);
        // VectorXd samples_sigD =  VectorXd::Constant(constant.nIter-constant.nBurn, 0);

        // *** instantiation "Generate"
        Generate generate(constant);
        if (constant.setting == 0)
        {
            // generate data
            generate.generateX(constant);
            generate.generateZ(constant);
            generate.generateM(constant);
            generate.generateY(constant);
            generate.generateNbh(constant);
            // write *.txt in ./out
            generate.writeX();
            generate.writeZ();
            generate.writeM();
            generate.writeY();
            generate.writeNbh();
        }

        // *** instantiation "Initial"
        Initial initial(constant);
        // set beta init at true value
        //initial.setTrueBeta();
        // set alpha and zetaA init at true value
        //initial.setTrueAlpha();

        // *** instantiation "Mcmc"
        Mcmc mcmc(constant);
        // read data and nbh
        mcmc.readData(constant, generate, rep);
        // set initial value
        mcmc.setInitial(constant, initial);
        // for fast computation
        mcmc.fastComp(constant);
        // update
        for (int g=0; g<constant.nIter; g++)
        {
            // print the iterator
            if ((g+1) % constant.nPrint == 0)
            {
                cout << "     | iteration: " << g+1 << endl;
            }

            // update alpha, zeta_alpha, a1, a2, sigma^2_epsilon
            mcmc.updateA(constant);
            // update beta, zeta_beta
            mcmc.updateB(constant);
            // update gamma
            mcmc.updateGamma(constant);
            // update b
            mcmc.updateb(constant);
            // update sigma^wettto
            // 2_delta, fix as a tuning if necessary
            //mcmc.updateSigD(constant);


            // store the posterior samples
            if (g > (constant.nBurn-1))
            {
                samples_alpha.row(g-constant.nBurn) = mcmc.getAlpha().transpose();
                samples_zetaA.row(g-constant.nBurn) = mcmc.getZetaA().transpose();
                samples_a1.row(g-constant.nBurn) = mcmc.getA1().transpose();
                samples_a2.row(g-constant.nBurn) = mcmc.getA2().transpose();
                samples_beta.row(g-constant.nBurn) = mcmc.getBeta().transpose();
                samples_zetaB.row(g-constant.nBurn) = mcmc.getZetaB().transpose();
                samples_gamma(g-constant.nBurn) = mcmc.getGamma();
                samples_b.row(g-constant.nBurn) = mcmc.getB();
                samples_sigE.row(g-constant.nBurn) = mcmc.getSigE().transpose();
                // samples_sigD(g-constant.nBurn) = mcmc.getSigD();
            }
        } // end loop g

        // estimation
        VectorXd zetaA_est = samples_zetaA.colwise().sum();
        VectorXd alpha_est = (samples_alpha.array()*samples_zetaA.array()).colwise().sum();
        alpha_est = alpha_est.array()/(constant.nIter - constant.nBurn);
        zetaA_est = zetaA_est/(constant.nIter-constant.nBurn);
        VectorXd a1_est = samples_a1.colwise().sum()/(constant.nIter-constant.nBurn);
        VectorXd a2_est = samples_a2.colwise().sum()/(constant.nIter-constant.nBurn);
        VectorXd zetaB_est = samples_zetaB.colwise().sum();
        VectorXd beta_est = (samples_beta.array()*samples_zetaB.array()).colwise().sum();
        beta_est = beta_est.array()/(constant.nIter - constant.nBurn);
        zetaB_est = zetaB_est/(constant.nIter-constant.nBurn);
        double gamma_est = samples_gamma.mean();
        VectorXd b_est = samples_b.colwise().sum()/(constant.nIter-constant.nBurn);
        VectorXd sigE_est = samples_sigE.colwise().mean();
        // double sigD_est = samples_sigD.mean();

        // record est, and quality for cross-validation
        rep_alpha.row(rep) = alpha_est.transpose();
        rep_zetaA.row(rep) = zetaA_est.transpose();
        rep_a1.row(rep) = a1_est.transpose();
        rep_a2.row(rep) = a2_est.transpose();
        rep_beta.row(rep) = beta_est.transpose();
        rep_zetaB.row(rep) = zetaB_est.transpose();
        rep_gamma(rep) = gamma_est;
        rep_b.row(rep) = b_est.transpose();
        rep_sigE.row(rep) = sigE_est.transpose();
        // rep_sigD(rep) = sigD_est;

        // write the samples, and cv quality of current replication
        ofstream file_alpha;
        file_alpha.open("../output/samples_alpha.txt");
        file_alpha << samples_alpha << endl;
        file_alpha.close();
        ofstream file_zetaA;
        file_zetaA.open("../output/samples_zetaA.txt");
        file_zetaA << samples_zetaA << endl;
        file_zetaA.close();
        ofstream file_a1;
        file_a1.open("../output/samples_a1.txt");
        file_a1 << samples_a1 << endl;
        file_a1.close();
        ofstream file_a2;
        file_a2.open("../output/samples_a2.txt");
        file_a2 << samples_a2 << endl;
        file_a2.close();
        ofstream file_beta;
        file_beta.open("../output/samples_beta.txt");
        file_beta << samples_beta << endl;
        file_beta.close();
        ofstream file_zetaB;
        file_zetaB.open("../output/samples_zetaB.txt");
        file_zetaB << samples_zetaB << endl;
        file_zetaB.close();
        ofstream file_gamma;
        file_gamma.open("../output/samples_gamma.txt");
        file_gamma << samples_gamma << endl;
        file_gamma.close();
        ofstream file_b;
        file_b.open("../output/samples_b.txt");
        file_b << samples_b << endl;
        file_b.close();
        ofstream file_sigE;
        file_sigE.open("../output/samples_sigE.txt");
        file_sigE << samples_sigE << endl;
        file_sigE.close();
        // ofstream file_sigD;
        // file_sigD.open("../output/samples_sigD.txt");
        // file_sigD << samples_sigD << endl;
        // file_sigD.close();
    } // end loop rep

    // write the est of all replications
    ofstream file_repAlpha;
    file_repAlpha.open("../output/rep_alpha.txt");
    file_repAlpha << rep_alpha << endl;
    file_repAlpha.close();
    ofstream file_repZetaA;
    file_repZetaA.open("../output/rep_zetaA.txt");
    file_repZetaA << rep_zetaA << endl;
    file_repZetaA.close();
    ofstream file_repA1;
    file_repA1.open("../output/rep_a1.txt");
    file_repA1 << rep_a1 << endl;
    file_repA1.close();
    ofstream file_repA2;
    file_repA2.open("../output/rep_a2.txt");
    file_repA2 << rep_a2 << endl;
    file_repA2.close();
    ofstream file_repBeta;
    file_repBeta.open("../output/rep_beta.txt");
    file_repBeta << rep_beta << endl;
    file_repBeta.close();
    ofstream file_repZetaB;
    file_repZetaB.open("../output/rep_zetaB.txt");
    file_repZetaB << rep_zetaB << endl;
    file_repZetaB.close();
    ofstream file_repGamma;
    file_repGamma.open("../output/rep_gamma.txt");
    file_repGamma << rep_gamma << endl;
    file_repGamma.close();
    ofstream file_repB;
    file_repB.open("../output/rep_b.txt");
    file_repB << rep_b << endl;
    file_repB.close();
    ofstream file_repSigE;
    file_repSigE.open("../output/rep_sigE.txt");
    file_repSigE << rep_sigE << endl;
    file_repSigE.close();
    // ofstream file_repSigD;
    // file_repSigD.open("../output/rep_sigD.txt");
    // file_repSigD << rep_sigD << endl;
    // file_repSigD.close();

    // ending time
    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    cout << "total time:" << total_t << endl;

    return 0;
}
