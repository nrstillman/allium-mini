// Created by N.R. Stillman & S. Henkes 2020
//
#ifndef CAPMD_PARTICLE_H
#define CAPMD_PARTICLE_H

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <iterator>     // for reading txt
#include <sstream>      // for reading txt
#include <algorithm>	// for 'transform' to generate iterators
#include <memory> 		// for shared pointers
/**
    \file Particle.h
    Class that mostly gives and returns individual particle information
*/
/*!
   Note, we do not include a destructor. This may need to be revisited early on
*/
class Particle
    {
    public:
        double radius; //!< Particle radius
        double age;  //!< Particle age (since creation) which is used in death calculation
        int type;  //!< Particle type (such as boundary cell)
        int  cellidx;
        std::array<double,2> position; //!< Particle position
        double theta;  //!< Particle angle
        std::array<double,2> force;  //!< Particle force
        std::array<double,2> vel;  //!< Particle velocity
        std::array<double,2> activeforce;  //!< Particle active force (which is output)
        std::vector<int> neighbours; //
        //! Particle constructor
        Particle(int pid = 0, int ptype = 0, std::array<double,2> px = {0,0},double ptheta = 0.,double  pr = 1.);
        //! Particle copy constructor
        Particle(const Particle &);
        //! Particle constructor overloaded with string (why?)
        Particle(std::string);
//
//        ~Particle(); // destructor

    /*!
      Various methods for accessing and altering particle properties
    */

        int getId() const { return id; };
        void setId( int x) {std::cout << "attempted id change!"<< std::endl;}

        int getCid() const { return cellidx; };
        void setCid(int x) {cellidx = x;}

        int getIndex() const { return idx; };
        void setIndex( int x) {idx = x;}

        int getType() const { return type; };
        void setType( int x) { type = x;}

        double getAge() const { return age;}
        void setAge( double x) { age = x;}

        double getRadius() const { return radius;}
        void setRadius( double x) { radius = x;}

        double getTheta() const { return theta;}
        void setTheta( double x) { theta = x;}

        int getNumNeigh() { return numneigh;}
        void setNumNeigh( int x) { numneigh= x;}

        int getZ() {return z;}
        void setZ(int x) { z= x;}
        void addZ(int x) { z+= x;};

        std::array<double,2> getPosition() { return position;}
        void setPosition(std::array<double,2> x) { position = x;}

        std::array<double,2> getPrevPosition() { return prevposition;}
        void setPrevPosition() { prevposition = position;}

        std::array<double,2> getForce() { return force;}
        void setForce(std::array<double,2> x)  {force = x;}
        void addForce(std::array<double,2>);

        std::array<double,2> getVel() { return vel;}
        void setVel(std::array<double,2> x)  { vel = x;}

        std::vector<int> getNeighbours() { return neighbours;}
        void setNeighbours(std::vector<int> x)  { neighbours = x;}

        void setActiveForce(std::array<double,2> x)  { activeforce = x;}

        friend std::ostream& operator<< (std::ostream &, const Particle &); //!< Overload the out stream for saving particle data
        friend std::istream& operator>> (std::istream &, const Particle &); //!< Overload the in stream for loading particle data

        static void split(const std::string &s, char delim, std::vector<double> &elems);  //!< Split function for ??

    private:
        int id; //!< Particle id which should not change
        int idx; //!< Particle position within the vector of all particles
        std::array<double,2> prevposition;  //!< Position of the particle for the previous timestep (used in neighbourlist calculation)
        int numneigh; //!< Number of neighbours
        double z; //!< Number of particles in contact with the particle
};

#endif //CAPMD_PARTICLE_H
