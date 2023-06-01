// Created by N.R. Stillman & S. Henkes 2020
//
#include "Population.h"
#include "Domain.h"
#include "Dynamics.h"
#include "Particle.h"
#include "Output.h"

#include "potentials/Interaction.h"
#include "potentials/Potential.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <algorithm> // fpr count
#include <memory> // For std::shared_ptr

#ifndef CAPMD_SIMULATION_H
#define CAPMD_SIMULATION_H
/**
    \file Simulation.h
    Main functioning and activity of the simulator which contains other classes
*/

class Simulation : virtual public Particle{

    private:
        Parameters params;

        // As container, has the different sub-pieces
        std::shared_ptr<Domain> domain;
        std::shared_ptr<Dynamics> dynamics;
        std::shared_ptr<Population> population;
//        std::shared_ptr<Interaction> interaction;
        std::shared_ptr<Potential> interaction;

        typedef std::mt19937 Engine;
        typedef std::uniform_real_distribution<double> Distribution;

        Engine gen;
        Distribution disx;
        Distribution disy;
        Distribution disr;
        Distribution distheta;

    public:
        std::shared_ptr<Output> output;

        // Access through interface (should make these protected)
        std::vector<std::shared_ptr<Particle>> particles;
        int cells;

        // Setting parameters
        void setParams(Parameters);

        // Methods for particle dynamics
        void move();
//        void move(int);
        void populationDynamics(int);

        // Methods for getting sim data
        int popSize(void){ return particles.size() - boundarysize;}
        int totalSize(void){ return particles.size();}
        int getBoundarySize(void){ return boundarysize;}
        std::string getFileName(){return params.filename;};

        // Methods for particle data
        std::shared_ptr<Particle> getParticle(int);
        std::shared_ptr<Particle> getParticlebyId(int);

        void removeParticle(int);
        void removeParticles(std::vector<int>);

        void changeParticles(std::vector<int>, int);


        std::vector<std::shared_ptr<Particle>> getAllParticles(void){return particles;};
        std::vector<std::array<double,2>> getPopulationPosition(std::vector<int>);
        std::vector<std::array<double,2>> getPopulationVelocity(std::vector<int>);
        std::vector<int> getPopulationId(std::vector<int>);
        std::vector<int> getPopulationZ(std::vector<int>);
        std::vector<int> getPopulationType(std::vector<int>);
        std::vector<double> getPopulationTheta(std::vector<int>);

        std::vector<std::array<double,2>> getBoundaryPosition();
        std::vector<int> getBoundaryId();

        std::vector<double> getPopulationRadius(std::list<int> &index);

        //  I/O for sim data
        void loadPopulation(std::string);
		void saveData(std::string outtype);

        // Access through other classes

        // Members (taken from parameters but included incase changing through python script)
        // Number of particles
        int N;

        // time step
        double dt;
		// current time step
		int timestep;
        // neighbour list total move threshold (should be of the order of 0.5)
        double maxmove;
        // running counter on particle flags
        int currentflag;

        int boundarysize = 0;

        // Constructors:
        Simulation();
        Simulation(Parameters);

        // Methods to initialise the system, including creating the particle vector and the first NeighbourList
        void initialise();
        void initPopulation(); //  To initialise a population
        void initBoundary(); //  To create the boundary

        void updateOutput(); //  To create the boundary
};


#endif //CAPMD_SIMULATION_H
