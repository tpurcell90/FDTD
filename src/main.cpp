// #include <iostream>
// #include "FDTDField.hpp"
// #include "FDTDFieldTE.hpp"
// #include "FDTDFieldTM.hpp"

#include <FDTD_MANAGER/parallelFDTDField.hpp>
// #include <iomanip>

namespace mpi = boost::mpi;

int main(int argc, char const *argv[])
{
    // Initialize the boost mpi environment and communicator
    mpi::environment env;
    mpiInterface world;
    std::clock_t start;
    std::string filename;
    double duration= 0.0;
    start = std::clock();
    if (argc < 2)
    {
        std::cout << "Provide an input json file" << std::endl;
        exit(1);
    }
    else
    {
        // Take in the file name and strip out all comments for the parser
        filename = argv[1];
        if(world.rank() == 0)
            stripComments(filename);
        else
            filename = "stripped_" + filename;
    }
    world.barrier();
    if(world.rank() == 0)
        std::cout << "Reading input file " << argv[1] << "..." << std::endl;
    //construct the parser and pass it to the inputs
    boost::property_tree::ptree propTree;
    boost::property_tree::json_parser::read_json(filename,propTree);
    parallelProgramInputs IP(propTree, filename);
    world.barrier();
    if(world.rank() == 0)
         boost::filesystem::remove(filename) ;
    if(world.rank() == 0)
        std::cout << "I TOOK ALL THE INPUT PARAMETERS" << std::endl;

    double maxT = IP.tMax();

    parallelFDTDFieldReal FF(IP, world);
    if(world.rank() == 0)
        std::cout << "made" << std::endl;

    while(FF.getTime() < maxT)
        FF.step();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    for(int ii = 0; ii < world.size(); ii ++)
    {
        world.barrier();
        if(world.rank() == ii)
            std::cout << world.rank() << "\t" << duration<<std::endl;
    }
    // Output the final flux for all flux objects, Polarization terms are because of units for continuous boxes defined by the z component (TE off grid is E fields, TM H fields)
    for(auto & flux : FF.fluxArr() )
    {
        // flux->getFlux(FF.E_incd(), FF.H_incd(), FF.H_mn_incd());
        if(FF.Ez_)
            flux->getFlux(FF.E_incd(), FF.H_incd(), FF.H_mn_incd(), true);
        else
            flux->getFlux(FF.E_incd(), FF.H_incd(), FF.E_pl_incd(), false);
    }
    // Power outputs need incident fields for normalization
    for(auto & dtc : FF.dtcFreqArr())
    {
        if(dtc->pow())
            dtc->toFile(FF.E_incd(), FF.dt());
        else
            dtc->toFile();
    }
    // output scaling information
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    for(int ii = 0; ii < world.size(); ii ++)
    {
        world.barrier();
        if(world.rank() == ii)
            std::cout << world.rank() << "\t" << duration<<std::endl;
    }

    return 0;
}