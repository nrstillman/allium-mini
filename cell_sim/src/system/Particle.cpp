// Created by N.R. Stillman & S. Henkes 2020
//
#include <Particle.h>
#include <cmath>

// Particle constructor
Particle::Particle(int pid, int ptype, std::array<double,2> px, double ptheta, double pr)
{
    age = 0;
    force = {0,0};
    vel = {0,0};
    prevposition = {0,0};
    id = pid;
    cellidx = 0;
    type = ptype;
    position = px;
    theta = ptheta;
    radius = pr;
    numneigh = 0;
    z = 0;

    idx = 0;

//    std::cout << "Particle " << pid <<" Initialised" << std::endl;
}

Particle::Particle(const Particle & rhs)
{
    this->setId(rhs.getId()); // check removing this doesn't break code
    std::cout << "particle copied " << std::endl;
    cellidx = rhs.cellidx;
    age = rhs.age;
    force = rhs.force;
    vel= rhs.vel;
    type = rhs.type;
    position = rhs.position;
    prevposition = rhs.prevposition;
    theta = rhs.theta;
    radius = rhs.radius;
    numneigh = rhs.numneigh;
    z = rhs.z;
}

// Particle constructor
Particle::Particle(std::string line)
{
    std::vector<double> properties;
    Particle::split(line, '\t', properties);

    if (properties.size() == 0){std::cout << "Particle made with empty string";}

    id = properties.at(0);
    type = properties.at(1);
    age = properties.at(2);
    position = {properties.at(3), properties.at(4)};
    theta = properties.at(5);
    radius = properties.at(6);
    neighbours = {};
    numneigh = 0;
    z = 0;
}

void Particle::addForce(std::array<double,2> f){
    force[0] += f[0];
    force[1] += f[1];
}
std::ostream& operator<<(std::ostream& out,const Particle& p)
{
    return out << p.getId() << ',' << p.type << ','
        << std::setprecision(8) << p.radius << ','
        << std::setprecision(8) << p.position[0] << ','
        << std::setprecision(8) << p.position[1] << ','
        << std::setprecision(8) << p.vel[0] << ','
        << std::setprecision(8) << p.vel[1] << ','
        << std::setprecision(8) << cos(p.theta) << ','
        << std::setprecision(8) << sin(p.theta) << std::endl ;
}

// taken from http://stackoverflow.com/a/236803/248823
void Particle::split(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(stod(item));
    }
}
