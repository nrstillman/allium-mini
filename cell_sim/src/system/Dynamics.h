// Created by N.R. Stillman & S. Henkes 2020

#ifndef CAPMD_DYNAMICS_H
#define CAPMD_DYNAMICS_H

#include "Particle.h"
#include "Parameters.h"

#include <memory>
#include <random>
/**
    \file Dynamics.h
    Handles active propulsion and overdamped motion
*/

/*!
    Used for all dynamic motion, including updating particle position. Population dynamics is treated in Population.h
*/
class Dynamics {

    public:
        //! Dynamics Constructor
        Dynamics(Parameters);
        //! Moves a particle according to the force law and adds a active motion
        void step(std::shared_ptr<Particle>);

    private:
        //! Parameter object (specifically for assigning...)
        std::vector<double> factive; //!< Magnitude of the active force
        std::vector<double> zeta; //!< Substrate friction
        std::vector<double> tau; //!< Correlation time of the active motion
        std::vector<double> alignmentTorque; //!< Alignment Torque
        bool periodic = false; //!< Bool condition for periodic boundary conditions
        double Lx, Ly; //!< Boundary edges
        int NCells; //!< Number of cells for computing neighbourhood
        double dt;

        typedef std::mt19937 Engine; //!< Mersenne-Twister Random Number Engine Template
        typedef std::normal_distribution<double> Distribution; //!< Normal distribution template

        Engine gen;  //!< Engine for random num gen
        Distribution dist; //!< Normal distribution used for updating theta
};

#endif //CAPMD_DYNAMICS_H
