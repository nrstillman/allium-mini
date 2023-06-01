// Created by N.R. Stillman & S. Henkes 2020
//
#include "Domain.h"
#include <cmath>
#include <iterator>

#define _USE_MATH_DEFINES

// Domain constructor
Domain::Domain(Parameters params){
    neighbour_print = params.neighbour_print;
    cutoff = params.cutoff;
    cutoffZ = params.cutoffZ;
    maxmove = params.maxmove;
    Lx = params.Lx;
    Ly = params.Ly;
    NCells  = params.NCells;
    std::vector<std::vector<int>> CellList(params.NCells*params.NCells);
    if (params.bc_opt == "periodic"){periodic = true;}
    // std::cout << "Initialised Domain" << std::endl;
}
//included as mod negative is machine depndent
int Domain::mod(int k, int n) {
    return ((k %= n) < 0) ? k+n : k;
}

// vector between two particles
std::array<double,2> Domain::calc_dr(std::array<double,2> xi, std::array<double,2> xj)
{
    double x = xj[0] - xi[0];
    double y = xj[1] - xi[1];
    if (periodic){
        if (x>=Lx/2){x-=Lx; }
        else if (x<=-Lx/2){x+=Lx; }
        if (y>=Ly/2){y-=Ly; }
        else if (y<=-Ly/2){y+=Ly; }
    }
    return {x, y};
}

// distance between two particles
double Domain::dist(std::array<double,2> i, std::array<double,2> j)
{
    std::array<double,2> dr = calc_dr(i,j);
    return sqrt(dr[0]*dr[0] + dr[1]*dr[1]);
}

// compute the actual number of neighbours here, based on the interaction range
// The cutoffZ here is *mandatorily* smaller than the neighbour list cutoff
// This is now done with each interaction calculation
int Domain::countZ(std::vector<std::shared_ptr<Particle>> particles, int i) {
    int z = 0;
    // get the particles which are in the local neighbour list
    std::list <std::shared_ptr<Particle>> neighs = NeighbourList[i];
    for (auto n : neighs) {
            double dist_ip = dist(particles[i]->getPosition(), n->getPosition());
            if (dist_ip <cutoffZ) {z+=1;}
    }
    return z;
}

// create the neighbour list
// gets passed all the currently existing particles
// and a suitable cutoff, which is *larger* than the maximum existing interaction range,
// optimal value is in the range of the first maximum of g(r), about 1.4 interaction ranges

void Domain::makeCellList(std::vector<std::shared_ptr<Particle>> particles){
    //std::cout << "Remade cell list" <<std::endl;
    //vector of previous positions of particle (used in rebuild)
    std::vector<std::vector<std::shared_ptr<Particle>>> tmpList(NCells*NCells);
    auto p = particles.begin();
    std::advance(p, boundarysize);
    int idx = 0;
    if (NCells < 9){
        std::cout << "domain is too small" << std::endl;
        throw;
    }
    while (p != particles.end()) {
        // construct grid and assign particle into cell
        double dx = Lx/NCells;
        double dy = Ly/NCells;
        int cid;
        for (int n = 0; n<=NCells;n++){
            for (int m = 0; m<=NCells;m++){
                bool xcheck = ((*p)->getPosition()[0] >= (-Lx/2 + m*dx) && ((*p)->getPosition()[0]< (-Lx/2 + dx*(m+1))));
                bool ycheck = ((*p)->getPosition()[1] >= (-Ly/2 + n*dy) && ((*p)->getPosition()[1]< (-Ly/2 + dy*(n+1))));
                if (xcheck && ycheck){
                    cid = NCells*m + n;
//                    std::cout << cid << std::endl;
                    tmpList[cid].push_back((*p));
                    (*p)->setCid(cid);
                }
            }
        }
        ++p;
    }
    CellList = tmpList;
}
void Domain::makeNeighbourList(std::vector<std::shared_ptr<Particle>> particles){
    NeighbourList.clear();
    //vector of previous positions of particle (used in rebuild)
    auto p = particles.begin();
    std::advance(p, boundarysize);
    int idx = 0;

    while (p != particles.end()) {
        /// Gets neighbouring cells for particle
        std::array<int, 9> cellNeighbours;
        int a = (int)(*p)->getCid()/NCells;
        int b = mod((*p)->getCid(),NCells);
        int cidx = 0;
        for (int i = a-1; i < a+2; i++) {
            for (int j = b-1; j < b+2; j++) {
                cellNeighbours[cidx] = NCells*(mod(i,NCells)) + mod(j,NCells);
                ++cidx;
            }
        }
        (*p)->setIndex(idx);
        (*p)->setPrevPosition();

        std::list<std::shared_ptr<Particle>> pneighs;
        std::vector<int> ids;
        int numneighs = 0;
        //TODO: Can we combine both of these into a single loop? - for now, have two
        for (cidx = 0; cidx < 9; cidx++) {
            std::vector<std::shared_ptr<Particle>>::iterator neigh, end;
            for (neigh = CellList[cellNeighbours[cidx]].begin(), end = CellList[cellNeighbours[cidx]].end();
                neigh != end; ++neigh) {
                if ((*p)->getId() != (*neigh)->getId()) {
                    double dist_pq = dist((*p)->getPosition(), (*neigh)->getPosition());
                    if (dist_pq < cutoff) {
                        pneighs.push_back((*neigh));
                        ids.push_back((*neigh)->getId());
                        numneighs += 1;
                    }
                }
            }
        }
        NeighbourList.push_back(pneighs);
        (*p)->setNumNeigh(numneighs);
        if (neighbour_print > 0){
            (*p)->setNeighbours(ids);
        }
		/* std::cout << "particle " << idx << " with id " << (*p)->getId() << " has " << numneighs << " in the neighbour list " << std::endl;
		std::cout << "neighbourcomputeForces are:" << std::endl;
		for (int n = 0; n<numneighs; n++){
		    std::cout << " " << ids[n];
		}
		std::cout << "." <<std::endl;
		*/
        idxmap[(*p)->getId()] = idx + boundarysize;

        pneighs.clear();
        ++p;
        idx += 1;
    }
}

// check for a neighbour list rebuild based on max motion of particles
bool Domain::checkRebuild(std::vector<std::shared_ptr<Particle>> particles) {
    for  (int i = boundarysize; i< particles.size(); ++i)  {
        // this is not pretty
//        std::array<double,2> drmove = calc_dr(particles[i]->getPosition(), particles[i]->getPrevPosition());
//        double distmove = sqrt(drmove[0]*drmove[0] + drmove[1]*drmove[1]);
        double distmove = dist(particles[i]->getPosition(), particles[i]->getPrevPosition());
        if (distmove > maxmove){
//            std::cout << "Rebuild neighbour list triggered!" << std::endl;
            return true;
        }
    }
    return false;
}

// return the list of neighbours of particle i
// cannot be used to get boundary cell neighbours (which aren't stored)
std::list<std::shared_ptr<Particle>> Domain::getNeighbours(int i) {
    return NeighbourList[i];
}
