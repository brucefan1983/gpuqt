/*
    Copyright 2017 Zheyong Fan, Ville Vierimaa, and Ari Harju

    This file is part of GPUQT.

    GPUQT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUQT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUQT.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "lsqt.h"
#include <iostream>
#include <fstream>
using namespace std;


//----------------------------------------------------------------------------80
int main(int argc, char *argv[])
{
    cout << endl;
    cout << "***************************************************************\n";
    cout << "*                  Welcome to use LSQT                        *\n";
    cout << "*          (Linear Scaling Quantum Transport)                 *\n";
    cout << "*        (Author:  Zheyong Fan <brucenju@gmail.com>)          *\n";
    cout << "***************************************************************\n";
    cout << endl;
	
    if (argc != 2)
    {
        cout << "Usage: src/gpuqt input.txt" << std::endl;
        exit(1);
    }
	
    ifstream input(argv[1]); // input = the driver input file
    if (!input.is_open())
    {
        cout << "Failed to open " << argv[1] << endl;
        exit(1);
    }		

    string directory;
    while (getline(input, directory))
    {
        if (directory == "")
             continue;
        cout << endl;
        cout << "===========================================================\n";
        cout << "Run LSQT simulation for " << directory << std::endl; 
        cout << "===========================================================\n";

        clock_t time_begin = clock();
         
        // call the driver function
        lsqt(directory);

        clock_t time_finish = clock();
        double time_used = double(time_finish - time_begin) / CLOCKS_PER_SEC;

        cout << endl;
        cout << "===========================================================\n";
        cout << "Total time used for " << directory << " = " 
             << time_used <<" s" << endl; 
        cout << "===========================================================\n";
    }

    return 0;
}


