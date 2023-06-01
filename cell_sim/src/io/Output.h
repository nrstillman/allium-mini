// Created by N.R. Stillman & S. Henkes 2020
//
#ifndef CAPMD_OUTPUT_H
#define CAPMD_OUTPUT_H

// #include <vtkCellArray.h>
// #include <vtkDoubleArray.h>
// #include <vtkXMLPolyDataWriter.h>
// #include <vtkPolyData.h>
// #include <vtkSmartPointer.h>
// #include <vtkPoints.h>
// #include <vtkPointData.h>
#include <chrono>
#include <cmath>

#include "Parameters.h"
#include "Particle.h"

class Output{

    typedef std::chrono::steady_clock Clock;

private:
    double elapsed;

public:
    /// Constructor
    std::chrono::time_point<Clock>  last_log = std::chrono::steady_clock::now();

    Output(Parameters, int, std::vector<std::shared_ptr<Particle>>);

    std::string file_name;
    std::string output_folder;
    int t_final;
    std::vector<double> factive;

    int N;
    int NB;
    std::vector<std::shared_ptr<Particle>> particles;

    void update(Parameters, int, std::vector<std::shared_ptr<Particle>>);
    void log(int);
    // void vtp(int);
    void savePopulation(int);

    void setParticles(std::vector<std::shared_ptr<Particle>> x) { particles = x;}
};
#endif //CAPMD_OUTPUT_H
