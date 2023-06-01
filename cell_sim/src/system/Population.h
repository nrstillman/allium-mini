// Created by N.R. Stillman & S. Henkes 2020
//
#ifndef CAPMD_POPULATION_H
#define CAPMD_POPULATION_H

#include <vector>
#include <random>

#include "Parameters.h"
#include "Particle.h"

/**
    \file Population.h
    Handles all population dynamics including division and death
*/
class Population{
    public:
            //! Population Dynamics constructor
            Population(Parameters);

            /*! Test for particle division
            \param type - Particle type (eg boundary particles don't divide)
            \param z - Number of contact neighbours
            \param timeint - Time since last division
            */
            bool testDivide(int, double, double);

            /*! Test for particle division
            \param type - Particle type (eg boundary particles don't divide)
            \param z - Number of contact neighbours
            \param timeint - Time since last division
            */
            bool testDeath(int, double);

    private:

        double alpha;//!< division constant from data
        int ntypes;//!< Number of types
        std::vector<double> divrates; //!< Division rate
        std::vector<double> deathrates; //!< death rate
        int maxZ; //!< maximum density at division (in units of particle neighbours)

        std::random_device rd;
        typedef std::mt19937 Engine;
        typedef std::uniform_real_distribution<double> Distribution;

        Engine gen;
        Distribution dist;
};


#endif //CAPMD_POPULATION_H