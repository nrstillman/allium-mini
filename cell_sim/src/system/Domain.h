// Created by N.R. Stillman & S. Henkes 2020

#ifndef CAPMD_DOMAIN_H
#define CAPMD_DOMAIN_H

#include <algorithm>    // std::transform
#include <functional>   // std::minus
#include <iostream>
#include <list>
#include <memory>
#include <vector>
#include <map> // for idx to id map

#include "Particle.h"
#include "Parameters.h"
/**
    \file Domain.h
    Contains boundary conditions and creates neighbour list
*/
/*!
    Essentially a container for NeighbourList, PrevPositions, etc.
*/
class Domain
    {
    public:
        double cutoff; //!< Outer cutoff for the particle
        double cutoffZ; //!< Cutoff for the particles contact
        int boundarysize;  //!< Maximum boundary size for the domain (where is this used?)
        double maxmove;  //!< Maximum distance of a particle before neighbour list is recalculated
        double Lx, Ly; //!< Boundary edges
        int NCells;
        bool periodic; //!< Indicating whether boundary is periodic or not
        int neighbour_print;
        //! Domain Constructor
        Domain(Parameters);

        //! Calculates vector between two particles
        std::array<double,2> calc_dr(std::array<double,2>, std::array<double,2> );

        //! Calculate mod (inc negative)
        int mod(int , int );

        //! Calculates distance between two particles (or doubles)
        double dist(std::array<double,2>, std::array<double,2> );

        //! Count number of particles in contact with a particle
        int countZ(std::vector<std::shared_ptr<Particle>>, int);

        //! Calculates a particles neighbour list - This is one of the most intensive aspects of code
        void makeNeighbourList(std::vector<std::shared_ptr<Particle>>);
        //! Calculates the cell list to assist neighbour list calculation
        void makeCellList(std::vector<std::shared_ptr<Particle>>);
        //! Check to see whether a particle's neighbour list needs to be recalculated based on distance travelled
        bool checkRebuild(std::vector<std::shared_ptr<Particle>>);
        //! Returns list of neighbours by index
        std::list<std::shared_ptr<Particle>> getNeighbours(int);
        //! Gets index of a particle using the position of the particle in the vector of all particles
        int getIdx(int x){ return idxmap[x];}

        //! Adjusts the boundary size (why is this used?)
        void setBoundarySize( int x) { boundarysize = x;}


private:
        //! Vector of all particles neighbours
        std::vector<std::list<std::shared_ptr<Particle>>> NeighbourList;
        //! Vector of all particles within a cell
        std::vector<std::vector<std::shared_ptr<Particle>>> CellList;
        //! Previous positions of all particle (used to calculate checkRebuild)
        std::vector<std::vector<double>> PrevPositions;
        //! Map from position in vectors of all particles to particle id
        std::map<int,int> idxmap;
};

#endif //CAPMD_DOMAIN_H
