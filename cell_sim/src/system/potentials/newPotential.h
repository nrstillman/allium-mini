// Created by N.R. Stillman & S. Henkes 2020
//
#include "Potential.h"
#include "Parameters.h"
#include <cmath>

#ifndef CAPMD_SOFT_H
#define CAPMD_SOFT_H

class newPotential : public Potential {

    public:
        newPotential(Parameters);
        //! Computes the mechanical force between particles
        void computeForce(std::shared_ptr<Particle>, std::shared_ptr<Particle>);
};

#endif //CAPMD_SOFT_H
