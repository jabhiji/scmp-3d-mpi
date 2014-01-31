#ifndef DOMAINDECOMP_H
#define DOMAINDECOMP_H

// include files

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

// define a struct to store the beginning and ending node numbers inside a process
struct node_range
{
    int beg;
    int end;
};

#endif
