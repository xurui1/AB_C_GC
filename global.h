#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include <time.h>      //Call system time libraries to define the integer seed for random numbers
#include "smemory.h"  //Use my custom memory handling class
//#include "mpi.h"     //Use this for MPI parallel implimentation later

using namespace std;

#define Nr 50
#define Nz 50
#define ChainType 3
#define Pi 3.14159


double z1, z2, kappa, kappa2;
double bond_energy, N_bond;
double r_0;
int initial;
