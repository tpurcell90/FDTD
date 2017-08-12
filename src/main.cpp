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
    std::shared_ptr<mpiInterface> gridComm = std::make_shared<mpiInterface>();
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
        if(gridComm->rank() == 0)
            stripComments(filename);
        else
            filename = "stripped_" + filename;
    }
    gridComm->barrier();
    if(gridComm->rank() == 0)
        std::cout << "Reading input file " << argv[1] << "..." << std::endl;
    //construct the parser and pass it to the inputs
    boost::property_tree::ptree propTree;
    boost::property_tree::json_parser::read_json(filename,propTree);
    parallelProgramInputs IP(propTree, filename);
    gridComm->barrier();
    if(gridComm->rank() == 0)
         boost::filesystem::remove(filename) ;
    if(gridComm->rank() == 0)
        std::cout << "I TOOK ALL THE INPUT PARAMETERS" << std::endl;

    double maxT = IP.tMax();

    parallelFDTDFieldReal FF(IP, gridComm);
    if(gridComm->rank() == 0)
        std::cout << "made" << std::endl;

    int nSteps = int(std::ceil( IP.tMax_ / (IP.courant_/IP.res_) ) );
    for(int tt = 0; tt < nSteps; ++tt )
        FF.step();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    for(int ii = 0; ii < gridComm->size(); ii ++)
    {
        gridComm->barrier();
        if(gridComm->rank() == ii)
            std::cout << gridComm->rank() << "\t" << duration<<std::endl;
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
    for(int ii = 0; ii < gridComm->size(); ii ++)
    {
        gridComm->barrier();
        if(gridComm->rank() == ii)
            std::cout << gridComm->rank() << "\t" << duration<<std::endl;
    }

    return 0;
}