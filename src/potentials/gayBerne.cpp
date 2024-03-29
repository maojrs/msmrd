//
// Created by dibakma on 17.08.18.
//

#include <vec3.hpp>
#include "potentials/gayBerne.hpp"

namespace msmrd {
    /**
     * Class implementation for Gay Berne pair potential
     * @param a steric anisotropy of the particles: Length to breadth ratio
     * @param d ratio of binding energy in the side-by-side to end-to-end configuration
     * @param eps0 depth of the potential well
     * @param sig0 length of the short axis of the ellipsoid
     */
    gayBerne::gayBerne(double a, double d, double eps0, double sig0) : eps0(eps0), sig0(sig0) {
        chi = (a * a - 1) / (a * a + 1);
        chip = (std::sqrt(d) - 1) / (std::sqrt(d) + 1);
    }

    // Evaluate Gay Berne potential
    double gayBerne::evaluate(particle &part1, particle &part2) {
        vec3<double> u1 = part1.orientvector;
        vec3<double> u2 = part2.orientvector;
        vec3<double> r = relativePosition(part2.position, part1.position); //part1.position - part2.position;
        double rabs = r.norm();
        vec3<double> rhat = r / rabs;
        double ru1 = rhat * u1;
        double ru2 = rhat * u2;
        double u1u2 = u1 * u2;
        double rup = ru1 + ru2; //ru plus
        double rup2 = rup * rup;
        double rum = ru1 - ru2; //ru minus
        double rum2 = rum * rum;
        double denp = 1 + chi * u1u2; // denominator plus
        double denm = 1 - chi * u1u2; // denominator minus
        double denpp = 1 + chip * u1u2; // denominator prime plus
        double denpm = 1 - chip * u1u2; // denominator prime minus
        double epsbrak = (1 - chip / 2 * (rup2 / denpp + rum2 / denpm)); // bracket in the epsilon term
        double epsden = std::sqrt(denp * denm); //denominator in the epsilon term
        double eps =
                eps0 * epsbrak * epsbrak /
                epsden;  // epsilon function of the Gay-Berne potential, spelled out the square
        double sig = 1.0 / std::sqrt(1 - chi / 2 * (rup2 / denp + rum2 / denm));
        double frac = 1 / (rabs/sig0 - sig + 1); // bracket in the potential term
        double frac2 = frac * frac;
        double frac6 = frac2 * frac2 * frac2;
        double en_pot = eps * frac6 * (frac6 - 1);
        return en_pot;
    }


    // Returns force and torque exerted on the first particle by the second one.
    std::array<vec3<double>, 4> gayBerne::forceTorque(particle &part1, particle &part2) {
        vec3<double> u1 = part1.orientvector;
        vec3<double> u2 = part2.orientvector;
        vec3<double> r = relativePosition(part2.position, part1.position); //part1.position - part2.position;
        double rabs = r.norm();
        vec3<double> rhat = r / rabs;
        double ru1 = rhat * u1;
        double ru2 = rhat * u2;
        double u1u2 = u1 * u2;
        double rup = ru1 + ru2; //ru plus
        double rup2 = rup * rup;
        double rum = ru1 - ru2; //ru minus
        double rum2 = rum * rum;
        double denp = 1 + chi * u1u2; // denominator plus
        double denm = 1 - chi * u1u2; // denominator minus
        double denpp = 1 + chip * u1u2; // denominator prime plus
        double denpm = 1 - chip * u1u2; // denominator prime minus
        double epsbrak = (1 - chip / 2 * (rup2 / denpp + rum2 / denpm)); // bracket in the epsilon term
        double epsden = std::sqrt(denp * denm); //denominator in the epsilon term
        double eps = eps0 * epsbrak * epsbrak /
                epsden;  // epsilon function of the Gay-Berne potential, spelled out the square
        double sig = 1 / std::sqrt(1 - chi / 2 * (rup2 / denp + rum2 / denm));
        double sig3 = sig * sig * sig;
        double frac = 1 / (rabs / sig0 - sig + 1); // bracket in the potential term
        double frac2 = frac * frac;
        double frac6 = frac2 * frac2 * frac2;
        double dU_drabs = -6 * eps * frac * frac6 * (2 * frac6 - 1) / sig0;
        double dU_dru1 = -2 * eps0 * frac6 * (frac6 - 1) * epsbrak * chip * (rup / denpp + rum / denpm) / epsden -
                         dU_drabs / 2 * sig3 * chi * (rup / denp + rum / denm);
        double dU_dru2 = -2 * eps0 * frac6 * (frac6 - 1) * epsbrak * chip * (rup / denpp - rum / denpm) / epsden -
                         dU_drabs / 2 * sig3 * chi * (rup / denp - rum / denm);
        double dU_du12 = (chi * chi * u1u2 / (epsden * epsden * epsden) * epsbrak * epsbrak +
                          chip * chip / epsden * epsbrak * (rup2 / (denpp * denpp) - rum2 / (denpm * denpm))) * frac6 *
                         (frac6 - 1) + dU_drabs * chi * chi * sig3 * (rup2 / (denp * denp) - rum2 / (denm * denm)) / 4;

        vec3<double> force = -dU_drabs * rhat - (dU_dru1 * (u1 - ru1 * rhat) + dU_dru2 * (u2 - ru2 * rhat)) / rabs;
        vec3<double> torque1 = dU_dru1 * rhat.cross(u1) - dU_du12 * u1.cross(u2);
        vec3<double> torque2 = dU_dru2 * rhat.cross(u2) - dU_du12 * u2.cross(u1);
        return {force, torque1, -1*force, torque2};
    }
}