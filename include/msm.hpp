//
// Created by maojrs on 7/4/18.
//

#pragma once
#include <array>
#include <algorithm>
#include "particle.hpp"
#include "randomgen.hpp"

namespace msmrd {
    /**
     * Abstract base class for Markov state models of particles
     */
    class baseMarkovStateModel {
    protected:
        int msmid;
        const long double tolerance = 1 * pow(10.0L, -10);
        long seed;
        randomgen randg;
    public:
        double lagtime;
        std::vector<std::vector<double>> tmatrix;
        unsigned int nstates;
        std::vector<double> Dlist;
        std::vector<double> Drotlist;
        /**
         * @param msmid ID of the msm, corresponds to the particle type
         * @param tolerance tolerance limit for MSM integrity check
         * @param seed variable for random number generation; seed = -1 corresponds to random_device;
         * @param randg random number generator class based on mt19937
         * @param lagtime msm lagtime (in ctmsm it is calculated after each propagation step)
         * @param tmatrix transition matrix (for ctmsm transition rate matrix)
         * @param nstates number of states in the msm (obtained directly from matrix size)
         * @param Dlist list of diffusion coefficients for each state
         * @param Drotlist list of rotational diffusion coefficients for each state
         */

        // Base constructor, can receive std::vectior matrix or numpy array matrix (through pybind)
        baseMarkovStateModel(int msmid, std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed);

        // Main functions definitions (=0 for abstract class)
        virtual void propagate(particleMS &part, int ksteps) = 0;

        // Get and set functions (**some needed for pybindinng)
        int getID() const { return msmid; }

        int getNstates() const { return nstates; }

        double getLagtime() const { return lagtime; }

        std::vector<std::vector<double>> getTmatrix() const { return tmatrix; }

        void setD(std::vector<double> &D) {
            Dlist.resize(nstates);
            Dlist = D;
        }

        void setDrot(std::vector<double> &Drot) {
            Drotlist.resize(nstates);
            Drotlist = Drot;
        }
    };


    /**
     * Discrete-time Markov state model class declaration
     */
    class discreteTimeMarkovStateModel : public baseMarkovStateModel {
    public:
        discreteTimeMarkovStateModel(int msmid, std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed);

        void propagate(particleMS &part, int ksteps) override;
    };


    /**
     * Continuous time Markov state model (ctmsm) class declaration
     */
    class continuousTimeMarkovStateModel : public baseMarkovStateModel {
    private:
        std::vector<double> lambda0;
        std::vector<std::vector<double>> ratescumsum;

        void calculateParameters();

    public:
        continuousTimeMarkovStateModel(int msmid, std::vector<std::vector<double>> &tempmatrix, long seed);

        void propagate(particleMS &part, int ksteps) override;

        void propagateNoUpdate(particleMS &part, int ksteps);

        // Getters for testing purposes
        std::vector<double> getLambda0() const { return lambda0; };

        std::vector<std::vector<double>> getRatescumsum() const { return ratescumsum; };
    };

}






