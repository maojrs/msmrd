//
// Created by maojrs on 10/23/19.
//

#include "particleCompound.hpp"

namespace msmrd {
/**
     * Implementation of particleCompound class
     */

    /* Constructors */
    particleCompound::particleCompound(vec3<double> position) : position(position) { };

    
    particleCompound::particleCompound(std::vector<double> &position) : position(position) { };

    
    particleCompound::particleCompound(std::map<std::tuple<int,int>, int> boundPairsDictionary):
            boundPairsDictionary(boundPairsDictionary) { };

    particleCompound::particleCompound(vec3<double> position, std::map<std::tuple<int,int>, int> boundPairsDictionary):
            position(position), boundPairsDictionary(boundPairsDictionary) { };

    particleCompound::particleCompound(std::vector<double> &position,
                                       std::map<std::tuple<int,int>, int> boundPairsDictionary) :
            position(position), boundPairsDictionary(boundPairsDictionary) { };

    /* Deactivates compound, clears all dictionaries and sets active to false */
    void particleCompound::deactivateCompound() {
        relativePositions.clear();
        relativeOrientations.clear();
        boundPairsDictionary.clear();
        active = false;
    };

    /* Joins another particle complex into this particle complex. The local dictionary has preference if
     * equal keys. However, that should never happen. */
    void particleCompound::joinCompound(particleCompound partComplex) {
        boundPairsDictionary.insert(partComplex.boundPairsDictionary.begin(),
                                    partComplex.boundPairsDictionary.end());
    };

    /* Updates diffusion coefficients fo compound depending on the size of the compound. It
     * set D and Drot from the correct value of Dlist and Drot list.
     */
    void particleCompound::updateDiffusionCoefficients() {
        auto sizeOfCompound = getSizeOfCompound();
        if (sizeOfCompound >= 2) {
            D = Dlist[sizeOfCompound - 2];
            Drot = Drotlist[sizeOfCompound - 2];
        } else{
            throw std::domain_error("Size of active particle compound cannot be smaller than two.");
        }
    }


    // TODO: FUNCTION BELOW IN PROGESS

//    /* Splits compound if a pair is removed and if possible, returns resulting compounds */
//    std::vector<particleCompound> particleCompound::splitCompound(std::tuple<int,int> pairToBeRemoved) {
//        std::vector<particleCompound> particleCompounds;
//
//        // Remove pairToBeRemoved from boundsPairDictionary
//        boundPairsDictionary.erase (pairToBeRemoved);
//
//        /* Check if compound was broken (If we can still reach from one node to other of the pair removed,
//         * even after removing the pair, the compound was not broken.) */
//        int startingNode = std::get<0>(pairToBeRemoved);
//        std::tuple<int,int> key;
//        for (const auto &entry : boundPairsDictionary) {
//            key = entry.first;
//            if ( (std::get<0>(key) == startingNode) or (std::get<1>(key) == startingNode) ) {
//            }
//        }
//        }
//    };


}