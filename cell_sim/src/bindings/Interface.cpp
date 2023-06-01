//
// Created by Namid Stillman on 7/21/20.
//

#include "Interface.h"
// for convenience
using json = nlohmann::json;

Interface::Interface() : Simulation(){
};

Interface::Interface(Parameters params) : Simulation(params){
};

void Interface::trackCell(int cell_idx){
    std::shared_ptr<Particle> p = Simulation::getParticle(cell_idx);
    p->setType(99);
}

void Interface::killCell(int cell_id){
    Simulation::removeParticle(cell_id);
    std::cout << "Cell with id " << cell_id << " killed" << std::endl;
}

void Interface::killCells(std::vector<int> cell_ids){
    Simulation::removeParticles(cell_ids);
}

void Interface::setCellType(int cell_idx, int new_cell_type){
    std::shared_ptr<Particle> p = Simulation::getParticle(cell_idx);
    p->setType(new_cell_type);
};

void Interface::setCellTypes(std::vector<int> cell_ids, int new_cell_type){
    Simulation::changeParticles(cell_ids, new_cell_type);
};

Parameters Interface::loadJSON(std::string filename){
    // read a JSON file
    std::ifstream f("../include/config/" + filename);

    std::string line;
    try {
        f.is_open();
        std::cout << "File found @ " << "'../include/config/" + filename + "'." << std::endl;
        std::cout << "---------" << std::endl;
    }
    catch (int e){
        std::cout << "Error " << e << ". Could not find file. ";
    }

    try
    {
        // parsing input with a syntax error
        json j = json::parse(f);
        // range-based for
        // special iterator member functions for objects
//        TODO: Make this case insenstive
        Parameters params;
        for (json::iterator it = j.begin(); it != j.end(); ++it) {
            if (it.key() == "filename"){
                params.filename = it.value();
            }
            else if (it.key() == "outputfolder"){
                params.outputfolder = it.value();
            }
            else if (it.key() == "neighbour_print"){
                params.neighbour_print = it.value();
            }
            else if (it.key() == "output_time"){
                params.output_time = it.value();
            }
            else if (it.key() == "output_type"){
                params.output_type = it.value();
            }
            else if (it.key() == "init_opt"){
                params.init_opt = it.value();
            }
            else if (it.key() == "bc_opt"){
                params.bc_opt = it.value();
            }
            else if (it.key() == "popdynfreq"){
                params.popdynfreq = it.value();
            }
            else if (it.key() == "zaptime"){
                params.zaptime = it.value();
            }
            else if (it.key() == "t_final"){
                params.t_final = it.value();
            }
            else if (it.key() == "dt") {
                params.dt = it.value();
            }
            else if (it.key() == "N") {
                params.N = it.value();
            }
            else if (it.key() == "NTA") {
                params.NTA = it.value();
            }
            else if (it.key() == "Ntracer") {
                params.Ntracer = it.value();
            }
            else if (it.key() == "ntypes") {
                params.ntypes = it.value();
            }
            else if (it.key() == "phi"){
                params.phi = it.value();
            }
            else if (it.key() == "R"){
                params.R = it.value();
            }
            else if (it.key() == "Lx"){
                params.Lx = it.value();
            }
            else if (it.key() == "Ly"){
                params.Ly = it.value();
            }
            else if (it.key() == "circle_r"){
                params.circle_r = it.value();
            }
            else if (it.key() == "initseed"){
                params.initseed = it.value();
            }
            else if (it.key() == "fade"){
                params.fade = it.value();
            }
            else if (it.key() == "maxZ"){
                params.maxZ = it.value();
            }
            else if (it.key() == "cutoffZ"){
                params.cutoffZ = it.value();
            }
            else if (it.key() == "deathrate"){
                params.deathrate[0] = it.value()[0];
                params.deathrate[1] = it.value()[1];
                params.deathrate[2] = it.value()[2];
            }
            else if (it.key() == "divrate"){
                params.divrate[0] = it.value()[0];
                params.divrate[1] = it.value()[1];
                params.divrate[2] = it.value()[2];
            }
            else if (it.key() == "factive"){
                params.factive[0] = it.value()[0];
                params.factive[1] = it.value()[1];
                params.factive[2] = it.value()[2];
            }
            else if (it.key() == "tau"){
                params.tau[0] = it.value()[0];
                params.tau[1] = it.value()[1];
                params.tau[2] = it.value()[2];
            }
            else if (it.key() == "alignmentTorque"){
                params.alignmentTorque[0] = it.value()[0];
                params.alignmentTorque[1] = it.value()[1];
                params.alignmentTorque[2] = it.value()[2];
            }
            else if (it.key() == "pairatt"){
                params.pairatt[0][0] = it.value()[0][0];
                params.pairatt[0][1] = it.value()[0][1];
                params.pairatt[0][2] = it.value()[0][2];
                params.pairatt[1][0] = it.value()[1][0];
                params.pairatt[1][1] = it.value()[1][1];
                params.pairatt[1][2] = it.value()[1][2];
                params.pairatt[2][0] = it.value()[2][0];
                params.pairatt[2][1] = it.value()[2][1];
                params.pairatt[2][2] = it.value()[2][2];
            }
            else if (it.key() == "pairstiff"){
                params.pairstiff[0][0] = it.value()[0][0];
                params.pairstiff[0][1] = it.value()[0][1];
                params.pairstiff[0][2] = it.value()[0][2];
                params.pairstiff[1][0] = it.value()[1][0];
                params.pairstiff[1][1] = it.value()[1][1];
                params.pairstiff[1][2] = it.value()[1][2];
                params.pairstiff[2][0] = it.value()[2][0];
                params.pairstiff[2][1] = it.value()[2][1];
                params.pairstiff[2][2] = it.value()[2][2];
            }
            else{
                std::cout << "\tNote:" << it.key() << " not updated." << std::endl;
            }
        }
        std::cout << "---------" << std::endl;
        return params;
    }
    catch (json::parse_error& e)
    {
        // output exception information
        std::cout << "\nmessage: " << e.what() << '\n'
                  << "\nexception id: " << e.id << '\n'
                  << "\nbyte position of error: " << e.byte << std::endl;
    }
}
//void setCellAttrib(int cell_idx, int attribute_value, std::string attribute){};
