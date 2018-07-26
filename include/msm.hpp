//
// Created by maojrs on 7/4/18.
//

#pragma once
#include <array>
#include <algorithm>
#include "particle.hpp"
#include "randomgen.hpp"


/**
 * Abstract base class for Markov state models
 */
class msmbase {
protected:
    const long double tolerance = 1*pow(10.0L,-10);
    long seed;
    randomgen randg;
public:
    int msmid;
    double lagtime;
    std::vector<std::vector<double>> tmatrix;
    unsigned int nstates;
    std::vector<double> Dlist;
    std::vector<double> Drotlist;
    /**
     * Constructors
     * @param msmid ID of the msm, corresponds to the particle type
     * @param lagtime msm lagtime (in ctmsm it is calculated after each propagation step)
     * @param tmatrix transition matrix (for ctmsm transition rate matrix)
     * @param nstates number of states in the msm (obtained directly from matrix size)
     * Additional variables
     * @param Dlist list of diffusion coefficients for each state
     * @param Drotlist list of rotational diffusion coefficients for each state
     * @param tolerance tolerance limit for MSM integrity check
     * @param seed variable for random number generation; seed = -1 corresponds to random_device;
     * @param randg random number generator class based on mt19937
     */

    // Base constructor, can receive std::vectior matrix or numpy array matrix (through pybind)
    msmbase(int msmid,  std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed)
            : msmid(msmid), lagtime(lagtime), seed(seed){

        // Resize vectors by input matrix size and set seed of random number generator
        nstates = static_cast<int>(tempmatrix.size());
        Dlist.resize(nstates);
        Drotlist.resize(nstates);
        tmatrix.resize(nstates);
        randg.setSeed(seed);

        // Verify input 2D vector is a square matrix and fill tmatrix
        for (const auto& row : tempmatrix) {
            if (tempmatrix.size() != row.size()) {
                throw std::range_error("MSM matrix must be a square matrix");
            }
        }
        for (int i=0; i<nstates; i++) {
            tmatrix[i].resize(nstates);
            std::copy_n(tempmatrix[i].begin(), nstates, tmatrix[i].begin());
        }
    };

    // Main functions definitions (=0 for abstract class)
    virtual void propagate(particle &part, int ksteps) = 0;
    // Get and set functions (**some needed for pybindinng)
    int getID() { return  msmid; }
    int getNstates() { return  nstates; }
    double getLagtime() {return lagtime; }
    std::vector<std::vector<double>> getTmatrix() { return  tmatrix; }
    void setD(std::vector<double> &D){
        Dlist.resize(nstates);
        Dlist=D;
    }
    void setDrot(std::vector<double> &Drot){
        Drotlist.resize(nstates);
        Drotlist=Drot;
    }
};


/**
 * Child classes of msmbase definitions: discrete-time (msm) and continuous-time (ctmsm)
 */

class msm: public msmbase {
public:
    msm(int msmid,  std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed);
    void propagate(particle &part, int ksteps) override;
};


class ctmsm: public msmbase {
private:
    std::vector<double> lambda0;
    std::vector<std::vector<double>> ratescumsum;
    void calculateParameters();
public:
    ctmsm(int msmid,  std::vector<std::vector<double>> &tempmatrix, long seed);
    void propagate(particle &part, int ksteps) override;
};






