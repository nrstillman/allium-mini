// Created by N.R. Stillman & S. Henkes 2020

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <random>

#include "Population.h"

Population::Population(Parameters params) {
// read out number of types, division and death rates
    alpha = params.alpha;
    ntypes = params.ntypes;
// type-specific division rates
    divrates = params.divrate;
// type-specific death rates
    deathrates = params.deathrate;
    // maximum density at division (in units of particle neighbours)
    maxZ = params.maxZ;

// it has its own random number generator
    gen = Engine(params.popseed);
    dist = Distribution(0,1);
}


// check if the particle will divide
// takes a particle, a number of neighbours for it, and a random number generator
// the time interval is the amount of real time that has elapsed since the last check
bool Population::testDivide(int type, double z, double timeint) {
    // actual density-dependent division rate
    double divreal = divrates[type]*(1-z/maxZ);
    // new density-dependent division rate calculated from data
    //    double divreal = divrates[type]*std::exp(alpha*z);

    // include timeinterval to make this a division probability.
    // Note that the resulting number needs to remain << 1 for accuracy
    divreal = divreal*timeint;
	bool div;
    // final check on division: If divreal is larger than a randomly chosen number in (0,1)
    // divid, else don't
    if (divreal > dist(gen)) div = true; else div= false;
             return div;
}

// check if a particle will die. Similar to division, just without the z-dependence (which can be added if desired)
bool Population::testDeath(int type, double timeint) {
    // include timeinterval to make this a division probability.
    // Note that the resulting number needs to remain << 1 for accuracy
    double deathreal = deathrates[type]*timeint;
    // final check on death: If deathreal is larger than a randomly chosen number in (0,1)
    // divid, else don't
    if (deathreal > dist(gen)) return true; else return false;
}
