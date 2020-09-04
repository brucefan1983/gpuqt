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

/*----------------------------------------------------------------------------80
    The main function of the LSQT code
------------------------------------------------------------------------------*/

#include "lsqt.h"
#include <fstream>
#include <iostream>
using namespace std;

static void print_welcome();
static void check_argc(int);
static void print_start(std::string);
static void print_finish(std::string, real);

int main(int argc, char* argv[])
{
  print_welcome();
  check_argc(argc);
  ifstream input(argv[1]); // input = the driver input file
  if (!input.is_open()) {
    cout << "Failed to open " << argv[1] << endl;
    exit(1);
  }
  string directory;
  while (getline(input, directory)) {
    if (directory == "") {
      continue;
    }
    print_start(directory);
    clock_t time_begin = clock();
    lsqt(directory);
    clock_t time_finish = clock();
    real time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
    print_finish(directory, time_used);
  }
  return 0;
}

static void print_welcome()
{
  cout << endl;
  cout << "***************************************************************\n";
  cout << "*                  Welcome to use LSQT                        *\n";
  cout << "*          (Linear Scaling Quantum Transport)                 *\n";
  cout << "*        (Author:  Zheyong Fan <brucenju@gmail.com>)          *\n";
  cout << "***************************************************************\n";
  cout << endl;
}

static void check_argc(int argc)
{
  if (argc != 2) {
    cout << "Usage: src/gpuqt input.txt" << std::endl;
    exit(1);
  }
}

static void print_start(std::string directory)
{
  cout << endl;
  cout << "===============================================================\n";
  cout << "Run LSQT simulation for " << directory << std::endl;
  cout << "===============================================================\n";
}

static void print_finish(std::string directory, real time)
{
  cout << endl;
  cout << "===============================================================\n";
  cout << "Total time used for " << directory << " = " << time << " s" << endl;
  cout << "===============================================================\n";
}
