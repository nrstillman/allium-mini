//
// Created by Namid Stillman on 10/25/20.
//
#ifndef CAPMD_INTERACTION_H
#define CAPMD_INTERACTION_H

#include "Parameters.h"
#include "Potential.h"

class Interaction {

    public:
        enum Potentials {
            soft, adh, hard
        };

        Interaction();

        static std::shared_ptr<Potential> createPotential(Parameters);
};

#endif //CAPMD_INTERACTION_H
