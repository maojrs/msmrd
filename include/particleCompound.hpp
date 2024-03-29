//
// Created by maojrs on 10/23/19.
//

#pragma once
#include <map>
#include "quaternion.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     * Declaration of particle Compound class. This will easily keep track which particles
     * are bound together in multiparticle MSM/RD, and their position with respect to the compound.
     * WORK IN PROGRESS
     */
    class particleCompound {
    public:
        double D;
        double Drot;
        vec3<double> position;
        quaternion<double> orientation = quaternion<double>(1.0, 0.0, 0.0, 0.0);
        std::map<std::tuple<int,int>, int> boundPairsDictionary = {};
        std::map<int, vec3<double>> relativePositions = {};
        std::map<int, quaternion<double>> relativeOrientations = {};
        int referenceParticleIndex = -1;
        bool active = true;
        std::vector<double> Dlist;
        std::vector<double> Drotlist;
        /**
         * @param D diffusion constant
         * @param Drot rotational diffusion constant
         * @param position position of particle compound (on chosen center, e.g. a center of mass or a geometric center)
         * @param positonReference position of reference particle within compound. The other
         * particles postion can be built from this value, the compound position and the bound states.
         * @param orientation orientation of particle compound. It will always be created with identity
         * orientation. Since the particle orientations can be reconstructed from reference particle.
         * @param orientationReference orientation of reference particle within compund. The other
         * particles orientation can be built from this value, the compound orientation and the bound states.
         * @param boundPairsDictionary dictionary, whose key corresponds to the tuple of pairs of indexes of the
         * particles bound with each other within the compound. The value is the state in which
         * the corresponding pair of particles is bound in.
         * @param relativePositions dictionary with the particle index as key and its relative position
         * with respect to the center of the particle compound.
         * @param relativeOrientations dictionary with the particle index as key and its relative orientation
         * with respect to the orientation of the main particle in the particle compound.
         * @param referenceParticleIndex index of reference particle in particleList. The reference particle will be chosen
         * as the one with smallest index when the compund is created. Once the compound is created it will not
         * change reference particle, unless the reference particle unbounds from complex. If two compounds join
         * the reference particle will remain the one with the smallest index.
         * compounds. If 0, it means compound is inactive.
         * @param active if true the compound is active, if false it is no longer active and can be deleted
         * by integrator to free up memory.
         * @param Dlist list of diffusion coefficients. Compounds of different sizes will diffuse
         * with different diffusion coefficients.
         * @param Drotlist same as Dlist, but for rotational diffusion.
         */

        /* Constructors for particle complex. */

        particleCompound() {};

        particleCompound(vec3<double> position);

        particleCompound(std::vector<double> &position);

        particleCompound(std::map<std::tuple<int,int>, int> boundPairsDictionary);

        particleCompound(vec3<double> position, std::map<std::tuple<int,int>, int> boundPairsDictionary);

        particleCompound(std::vector<double> &position, std::map<std::tuple<int,int>, int> boundPairsDictionary);

        // Utility functions

        void deactivateCompound();

        void joinCompound(particleCompound partComplex);

        void updateDiffusionCoefficients();

        //std::vector<particleCompound> splitCompound(std::tuple<int,int> pairToBeRemoved);

        bool isActive() { return active; }

        int getSizeOfCompound() {return relativePositions.size(); }

        int getNumberOfbindings() {return boundPairsDictionary.size(); }

        // Setters

        void setD(double Dnew) { D = Dnew; }

        void setDrot(double Drotnew) { Drot = Drotnew; }

        void setDs(double Dnew, double Drotnew) { D = Dnew, Drot = Drotnew; }

        void setDlist(std::vector<double> newDlist) { Dlist = newDlist; }

        void setDrotlist(std::vector<double> newDrotlist) { Drotlist = newDrotlist; }
    };



}