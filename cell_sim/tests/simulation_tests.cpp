#define CATCH_CONFIG_MAIN
#include "catch.hh"
#include "Simulation.h"
#include "Parameters.h"
#include "Output.h"

TEST_CASE( "Check particle class", "[ParticleTests]" ) {


    Particle p =  Particle();

    SECTION("0 -> ParticleConstructors") {
        // initialise a particle
        Particle p0;
        CHECK(p0.getId() == 0);
        // initialise a particle with certain properties
        Particle p1(2, 2, std::array<double,2> {1,1}, 2, 2);
        CHECK(p1.getId() == 2);
        // copy a particle
        Particle p2;
        p2 = p1;
        REQUIRE(p2.getId() == 2);
    }

    SECTION("1 -> ParticleProperties") {

        // initialise with position
        Particle p (2, 2, std::array<double,2> {12,4}, 2, 2);
        std::array<double,2> x = p.getPosition();
        CHECK(x[0] == 12);
        CHECK(x[1] == 4);

        // set particle position
        p.setPosition(std::array<double,2>{5, 9});

        x = p.getPosition();
        CHECK(x[0] == 5);
        CHECK(x[1] == 9);
        }
    }


TEST_CASE( "Test Simulation Basics", "[SimulationTests]" ) {

    SECTION("0 -> initSim") {
        // initialise simulation from loaded parameters
        Parameters params;
        Simulation sim(params);
        sim.initialise();
        CHECK(sim.popSize() == 50);

        // new parameters with updated parameters
        params.N = 10;
        sim.setParams(params);
        sim.initialise();
        REQUIRE(sim.popSize() == 10);
    }


    SECTION("1 -> createPopulation") {

    }

    SECTION("2 -> createNeighbourList") {

    }

    SECTION("3 -> moveParticle") {
        Parameters params;
        Simulation sim = Simulation(params);
        sim.initialise();

        std::array<double,2> x_t0 = sim.getParticle(0)->getPosition();

        for (int t = 0; t<100;t++){
            sim.move();
        }

        std::array<double,2> x_t100 = sim.getParticle(0)->getPosition();
        REQUIRE(x_t0 != x_t100);
    }

    SECTION("4 -> particleInteraction") {

    }

    SECTION("5 -> particleDivision") {

    }
    SECTION("6 -> particleDeath") {

    }
}


TEST_CASE( "Test Boundaries", "[BoundaryTests]" ) {

    SECTION("0 -> bounded") {
    }

    SECTION("1 -> periodic") {

    }
}


TEST_CASE( "Input/output", "[InputOutputTests]" ) {

    SECTION("0 -> saveData") {
    }

    SECTION("1 -> loadData") {

    }

    SECTION("2 -> saveVTK") {

    }
}
