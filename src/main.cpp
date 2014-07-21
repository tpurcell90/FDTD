#include <iostream>
#include "Inputs.hpp"
#include "FDTDField.hpp"

using namespace std;
int main(int argc, char const *argv[])
{
    programInputs IP("/home/tap620/git/FDTD/inputs.json");
    FDTDField FF(IP);

    /*cout << "hello \n";
    cout << FF.n_x() <<"\n";
    cout << FF.n_y() <<"\n";
    cout << FF.d_x() <<"\n";
    cout << FF.d_y() <<"\n";
    cout << FF.d_t() <<"\n";
    */
    FF.initializeGrid(IP);
    for(int ii = 0; ii < 3; ii++)
        FF.step();
    cout << "I am always in error\n";
    return 0; 
}
 
