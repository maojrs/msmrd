//
// Created by maojrs on 7/4/18.
//

#pragma once
#include <array>
#include <algorithm>
#include "particle.hpp"
#include "randomgen.hpp"

/**
 * Abstract base class for Markov state models of particles
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
     * @param tolerance tolerance limit for MSM integrity check
     * @param seed variable for random number generation; seed = -1 corresponds to random_device;
     * @param randg random number generator class based on mt19937
     * @param msmid ID of the msm, corresponds to the particle type
     * @param lagtime msm lagtime (in ctmsm it is calculated after each propagation step)
     * @param tmatrix transition matrix (for ctmsm transition rate matrix)
     * @param nstates number of states in the msm (obtained directly from matrix size)
     * @param Dlist list of diffusion coefficients for each state
     * @param Drotlist list of rotational diffusion coefficients for each state
     */

    // Base constructor, can receive std::vectior matrix or numpy array matrix (through pybind)
    msmbase(int msmid,  std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed);

    // Main functions definitions (=0 for abstract class)
    virtual void propagate(particleMS &part, int ksteps) = 0;
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
 * Discrete-time (msm) class declaration
 */
class msm: public msmbase {
public:
    msm(int msmid,  std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed);
    void propagate(particleMS &part, int ksteps) override;
};


/**
 * Continuous time Markov state model (ctmsm) class declaration
 */
class ctmsm: public msmbase {
private:
    std::vector<double> lambda0;
    std::vector<std::vector<double>> ratescumsum;
    void calculateParameters();
public:
    ctmsm(int msmid,  std::vector<std::vector<double>> &tempmatrix, long seed);
    void propagate(particleMS &part, int ksteps) override;
};






