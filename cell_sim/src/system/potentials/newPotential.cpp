//
// Created by Namid Stillman on 10/25/20.
//
#include "newPotential.h"

newPotential::newPotential(Parameters params) : Potential(params){
}

void newPotential::computeForce(std::shared_ptr<Particle> i, std::shared_ptr<Particle> j) {
    std::cout << "Chosen new potential" << std::endl;

    // get pair parameters
    double kij = pairstiff[i->getType()][j->getType()];
    double eps = pairatt[i->getType()][j->getType()];

    // compute vector distance between particles
    std::array<double,2> dr = calc_dr(i->getPosition(),j->getPosition());
    // compute distance
    double dx = dist(i->getPosition(), j->getPosition());

    // actual force computation according to our potentials
    // several lines since piecewise defined
    std::array<double,2> force = {0,0};

    double bij = i->getRadius()+ j->getRadius();
    if (dx < bij*(1 + eps)) {
        i ->addZ(1);
        force = {-kij*(bij - dx)*(dr[0]*dr[0])/dx,-kij*(bij - dx)*(dr[1]*dr[1])/dx};
    }
    else if (dx < bij*(1 + 2*eps)){
        force = { kij*(bij - dx - 2*eps)*(dr[0]*dr[0])/dx, kij*(bij - dx - 2*eps)*(dr[1]*dr[1])/dx};
    }
    // multiply resulting force by amount of fade-in required. Cumulative if both particles are fading
    // used for particle fade-in post division.
    if (j->getType() != btype){

        double multi = 1.0;
        if (i->getAge()<fade) {multi = i->getAge()/fade;}
        double multj = 1.0;
        if (j->getAge()<fade) {multj = j->getAge()/fade;}

        double multiplier = multi*multj;

        std::transform(force.begin(), force.end(), force.begin(), [&multiplier](auto& c){return c*multiplier;});
    }
    // add force
    i->addForce(force);
}