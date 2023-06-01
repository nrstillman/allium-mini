//
// Created by Namid Stillman on 10/25/20.
//
#include "Interaction.h"
#include "adhPotential.h"
//#include "newPotential.h"

Interaction::Interaction(){};

std::shared_ptr<Potential> Interaction::createPotential(Parameters params) {
    switch (params.potential) {

        case 1:
            // std::cout << "Using original potential" << std::endl;
            return std::make_shared<Potential>(params);

        case 2:
            std::cout << "Chosen potential w large adhesive well" << std::endl;
            return std::make_shared<adhPotential>(params);

        default:
            std::cout << "Chosen original potentials" << std::endl;
            return std::make_shared<Potential>(params);

    }
}
