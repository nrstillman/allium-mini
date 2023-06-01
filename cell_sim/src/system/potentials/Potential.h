//
// Created by Namid Stillman on 10/25/20.
//
#ifndef CAPMD_POTENTIAL_H
#define CAPMD_POTENTIAL_H

#include <array>
#include <memory>
#include <vector>

#include "Parameters.h"
#include "Particle.h"

class Potential {

    public:
        Parameters p;
        // Constructor
        Potential(Parameters);

        int ntypes; //! Number of different types of particles that exist in the simulation
        bool periodic; //!< Indicating whether boundary is periodic or not
        double Lx, Ly; //!< Boundary edges
        double cutoffZ;
        int btype;

        std::vector<std::vector<double>> pairstiff; //! Particle stiffnesses ... (k_ij)
        std::vector<std::vector<double>> pairatt; //! Particle attraction strengths (epsilon_ij)
        std::vector<double> delta; //! size of attractive well (Del_i)

        double fade; //! Fade-in (or out) time for particle interactions

        //! Calculates vector between two particles
        std::array<double,2> calc_dr(std::array<double,2>, std::array<double,2> );

        //! Calculates distance between two particles (or doubles)
        double dist(std::array<double,2>, std::array<double,2> );

        //! Computes the mechanical force between particles
        void computeForce(std::shared_ptr<Particle>, std::shared_ptr<Particle>);
};

#endif //CAPMD_POTENTIAL_H
