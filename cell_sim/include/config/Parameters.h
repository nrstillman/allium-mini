// Created by N.R. Stillman & S. Henkes 2020

#ifndef CAPMD_PARAMS_H
#define CAPMD_PARAMS_H

#include <string>
#include <vector>
#include <array>

#include <iostream>
#include <cmath>

struct Parameters {

    // everything relevant, including
    // information about system size
    // dynamical parameters
    // division death parameters
    // interaction parameters

    // io
    std::string configfile = "parameters.json";

    // io
    std::string filename = "testing";
    std::string outputfolder = "output/";
    std::string output_type = "both"; //< options are vtp, text and both
    int output_time = 100;
    bool neighbour_print = false;
    bool log = false;

    // Random
    int initseed = 2022;
    int angseed = 2001;
    int popseed = 2021;

    //Particle
    int flag = 0;
    int type = 1;
    double R = 1;
    double circle_r = 5; // To initialise in circle
    double poly = 0.3; // Matching [Silke, H, Soft Matt. 2016]

    // Simulation
    int t_final = 100000;
    int popdynfreq = 101;
    int zaptime = 15000000;
    double dt  = 0.01;  // fixed (upper bound)
    int dim = 2;
    double maxmove = 0.25; // rescale based on r (including polydispersity)
    double cutoff = 3; // rescale based on r (including polydispersity)
    double eps = 1E-6; 

    //Population
    int N = 1000; // Population size
    std::string init_opt = "random_unif"; //"random_unif";
    std::string bc_opt = "periodic";//"bounded";
    int Ntracer = N*0.05;
    int NTA = N - Ntracer;
    int Nstem = 0;
    int btype = 0;  // type for boundary cells

    //Domain
    // x and y length
    double phi = 1.;
    double Lx = 50;
    double Ly = 50;

    int NCells = 20; //discretisation of domain into Ncells*Ncells for neighbourlist calc

    //Interaction
    int potential = 1;
    std::vector<double> delta = {0, 0, 0};
    int ntypes = 4;
    double fade = 3; //also run with 1

    std::vector<std::vector<double>> pairstiff =
        {
            {10, 10, 10,10},
            {10, 1, 1,1},
            {10, 1, 1,1},
            {10, 1, 1,1},
        };
    std::vector<std::vector<double>> pairatt =
        {
            { 0, 0, 0, 0},
            { 0, 0.2, 0.2,0.2}, //max 0.2
            { 0, 0.2, 0.2,0.2}, //max 0.2
            { 0, 0.2, 0.2,0.2}, //max 0.2
        };

    //Dynamics (active parameters)
    std::vector<double> factive = {0, 0.03, 0.03,0.03}; //active force
    std::vector<double> alignmentTorque = {0,0,0,0}; //alignment torque
    std::vector<double> zeta = {0, 1, 1, 1}; //friction
    std::vector<double> tau = {0, 40, 100,100}; //angular noise term

    // Division/death rates
    double alpha = 1; // parameter for exponential density-dependence of division rate (unused)
    std::vector<double> deathrate = {0.01,0.01,0.01,0.01};
    std::vector<double> divrate = {0.1,0.1,0.1,0.1};

    double cutoffZ = 1 + pairatt[type][type]; // (1 + 2*eps) - matching [Silke, H, Soft Matt. 2016]
    double maxZ = 6;

};

#endif //CAPMD_PARAMS_H
