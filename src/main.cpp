#include <iostream>
#include "Inputs.hpp"
#include "FDTDField.hpp"
#include "toBitMap.hpp"

using namespace std;
int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        cout << "Provide an input json file" << endl;
        exit(1);
    }
    programInputs IP(argv[1]);
    FDTDField FF(IP);
    cout<<"I MADE THE FDTDField" << endl;

    /*cout << "hello \n";
    cout << FF.n_x() <<"\n";
    cout << FF.n_y() <<"\n";
    cout << FF.d_x() <<"\n";
    cout << FF.d_y() <<"\n";
    cout << FF.d_t() <<"\n";
    */
    FF.initializeGrid();
    cout<<"AND INITIALIZED THE GRID!" << endl;
    int count = 0;
    while(FF.getTime() <= IP.tMax())
    {
        FF.step();
        //string filename = "fout/Ez/testImageT" + to_string(FF.getTime()) + ".bmp";
        //cout << filename << endl;
        //if (count++ % 20 == 0) fieldToBitMap(FF, filename);
    }
    cout << "I am always in error\n";
    return 0;
}

