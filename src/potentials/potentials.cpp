//
// Created by maojrs on 8/21/18.
//
#include <tuple>
#include "potentials/potentials.hpp"
#include "particle.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     *  All these functions are needed for PyBinding evaluate and forceTorque functions, since pyBind only accepts
     *  vectors as input. The different versions of the same functions depend on the AUXVARIABLES template that can
     *  be used for particles with no orientation, rod-like particles with their orientation described by one vectors,
     *  and by rigidbody particles with no symmetry and their orientation described by quaternions. AUXVARIABLES
     *  can also incorporate particle types to incorporate type-specific behavior in the potential function.
     */



}