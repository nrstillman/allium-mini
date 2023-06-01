//
// Created by Namid Stillman on 10/25/20.
//
#include "adhPotential.h"

adhPotential::adhPotential(Parameters params) : Potential(params){
}

void adhPotential::computeForce(std::shared_ptr<Particle> i, std::shared_ptr<Particle> j) {
    std::cout << "Chosen potential w large well" << std::endl;

    // get pair parameters
    double kij = pairstiff[i->getType()][j->getType()];
    double eps = pairatt[i->getType()][j->getType()];
    double del = delta[i->getType()];

    // compute vector distance between particles
    std::array<double,2> dr = calc_dr(i->getPosition(),j->getPosition());
    // compute distance
    double dx = dist(i->getPosition(), j->getPosition());

    // actual force computation according to our potentials
    // several lines since piecewise defined
    std::array<double,2> force = {0,0};

    double bij = i->getRadius()+ j->getRadius();
    if (dx/bij - 1 <= eps) {
        force = {kij*(dx - bij)*dr[0]/dx,kij*(dx - bij)*dr[1]/dx};
    }
    else if (dx/bij - 1 > eps + del*eps) {
        force = {eps * dr[0] / dx, eps * dr[1] / dx};
    }
    else if ((dx/bij - 1 > eps) && (dx/bij - 1 <= 2*eps)){
        force = {-kij*(dx - bij - 2*bij*eps)*dr[0]/dx, -kij*(dx - bij - 2*bij*eps)*dr[1]/dx};
    }
    if (dx/bij < 1 + 2*eps){
        i ->addZ(1);
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