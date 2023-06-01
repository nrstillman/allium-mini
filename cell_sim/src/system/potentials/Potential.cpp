//
// Created by Namid Stillman on 10/25/20.
//

#include "Potential.h"
#include <cmath>

#define _USE_MATH_DEFINES

double TOL = 1E-6;

Potential::Potential(Parameters params){
    // read out number of types, pair stiffnesses and pair attractions
    ntypes = params.ntypes;
    pairstiff = params.pairstiff;
    pairatt= params.pairatt;
    cutoffZ= params.cutoffZ;
    btype = params.btype;
    fade = params.fade;
    Lx = params.Lx;
    Ly = params.Ly;
    if (params.bc_opt == "periodic"){periodic = true;}
    // std::cout << "Initialised Potential" << std::endl;
}

// vector between two particles
std::array<double,2> Potential::calc_dr(std::array<double,2> xi, std::array<double,2> xj)
{
    double x = xj[0] - xi[0];
    double y = xj[1] - xi[1];
    if (periodic){
        if (x>=Lx/2){x-=Lx; }// std::cout << "new x is " << x << std::endl;}
        else if (x<=-Lx/2){x+=Lx; }// std::cout << "new x is " << x << std::endl;}
        if (y>=Ly/2){y-=Ly; }//std::cout << "new y is " << y << std::endl;}
        else if (y<=-Ly/2){y+=Ly; }// std::cout << "new y is " << y << std::endl;}
    }
    return {x, y};
}

// distance between two particles
double Potential::dist(std::array<double,2> i, std::array<double,2> j)
{
    std::array<double,2> dr = calc_dr(i,j);

    if (i == j) return TOL; else return sqrt(dr[0]*dr[0] + dr[1]*dr[1]);
}

void Potential::computeForce(std::shared_ptr<Particle> i, std::shared_ptr<Particle> j) {
//    std::cout << "Chosen old potential" << std::endl;
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

    double bij = i->getRadius() + j->getRadius();

    if (dx/bij - 1 <= eps) {
        force = {kij*(dx - bij)*dr[0]/dx,kij*(dx - bij)*dr[1]/dx};
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