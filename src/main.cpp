#include <iostream>
#include "Inputs.hpp"
#include "FDTDField.hpp"
#include "toBitMap.hpp"
#include <cstdio>
#include <ctime>


using namespace std;
int main(int argc, char const *argv[])
{
    clock_t start;
    double duration= 0.0;
    start = clock();

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
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout <<duration<<endl;
    cout << "I am always in error\n";
    return 0;
}

