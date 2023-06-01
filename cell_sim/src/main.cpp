// Created by N.R. Stillman & S. Henkes 2020

#include "Parameters.h"
#include "Interface.h"
#include "Simulation.h"
#include <chrono>
#include <vector>

int main(int argc, char* argv[]) {
//    cout << "Number of args = " << argc << endl;
    auto start = std::chrono::high_resolution_clock::now();

    Parameters params;

    if (argc > 1){
        params =  Interface::loadJSON(argv[1]);
    }
    Interface sim = Interface(params);

    std::array<double, 2> x ;
    double Rlength = params.Lx/4;
    std::cout << Rlength << std::endl;
    std::array<double, 2> maxR = {Rlength/2, double(params.Ly)};
    std::array<double, 2> minR = {-Rlength/2, -1*double(params.Ly)};

    for (int t = 0; t<= params.t_final; t++){
        sim.move();
        // Test for output
        if (t % params.output_time == 0) {
            //sim.updateOutput();
            sim.output->log(t);
            if (params.output_type.compare("all") == 0) {
                sim.saveData("text");
                sim.saveData("vtp");
            } else {
                sim.saveData(params.output_type);
            }
        }
        // Test for zap
        if (t == params.zaptime){
            int pop = sim.popSize();
            std::vector<int> popidx;
            for (int i = sim.boundarysize; i < sim.boundarysize + pop; i++) {popidx.push_back(i);}
            std::vector<int> popId = sim.getPopulationId(popidx);
            std::vector<std::array<double,2>> popPosn = sim.getPopulationPosition(popidx);
            std::vector<int> zapList;
            for (int i = 0; i < pop; i++){
                std::array<double,2> x = popPosn[i];
                if ((x[0] < maxR[0]) && (x[0] > minR[0])){
                    if ((x[1] < maxR[1]) && (x[1] > minR[1])){
                        zapList.push_back(popId[i]);
                    }
                }
            }
            sim.killCells(zapList);
            std::cout << "Cell zapping stage completed" << std::endl;
        }
        // Test for population dynamicss
        if (t % params.popdynfreq == 0){sim.populationDynamics(params.popdynfreq);}
    }
    return 0;
}

