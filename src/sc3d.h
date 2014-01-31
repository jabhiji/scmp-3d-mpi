#ifndef SC3D_H
#define SC3D_H

//    C++ headers

      #include <iostream>     // cout()
      #include <cmath>        // pow()
      #include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC
      #include <mpi.h>        // MPI 

//    data structures

//    define a struct to store the beginning and ending node numbers inside a process
      struct node_range
      {
        int beg;
        int end;
      };

//    set up MPI

      extern void mpiSetup (int argc,               // argc is short for argument count (set depending on how many arguments the user enters at the command line) 
                            char *argv[],           // pointer to a char array  ... argv is short for argument values ... array size set based on the value of argc
                            int* numprocs,          // pointer to an integer - number of distinct MPI processes on which this code will be executed 
                            int* myid,              // pointer to an integer - process ID
                            const int ndims,        // number of spatial dimensions for domain partitioning
                            int* dims,              // pointer to --> dims[0] - number of partitions of the domain along X, Y and Z...dims[0], dims[1], dims[2]
                            int* coords,            // pointer to --> coords[0] - coordinates of this process within the Cartesian topology ...coords[0], coords[1], coords[2]
                            MPI_Comm & CART_COMM,   // name of the Cartesian communicator
                            int* nbr_WEST,          // pointer to --> ID of neighboring process to my west   (i-1,j,k)
                            int* nbr_EAST,          // pointer to --> ID of neighboring process to my east   (i+1,j,k)
                            int* nbr_SOUTH,         // pointer to --> ID of neighboring process to my south  (i,j-1,k)
                            int* nbr_NORTH,         // pointer to --> ID of neighboring process to my north  (i,j+1,k)
                            int* nbr_BOTTOM,        // pointer to --> ID of neighboring process to my bottom (i,j,k-1)
                            int* nbr_TOP);          // pointer to --> ID of neighboring process to my top    (i,j,k+1)

void domainDecomp3D(// inputs
                    const int      & myid,           // MPI rank
                    const MPI_Comm & CART_COMM,      // MPI communicator name
                    const int      * dims,           // number of partitions of the domain along X, Y and Z, dims[0], dims[1] and dims[2]
                    const int      * coords,         // (X,Y,Z) "coordinates" of this partition
                    const int      & nodes_x, 
                    const int      & nodes_y, 
                    const int      & nodes_z,
                    const double   & delta,
                    const double   & x_min,
                    const double   & y_min,
                    const double   & z_min,
                    // outputs
                    node_range & x_range, 
                    node_range & y_range, 
                    node_range & z_range, 
                    double & local_origin_x, 
                    double & local_origin_y, 
                    double & local_origin_z,
                    int & LX,
                    int & LY,
                    int & LZ);

//    initialize all buffers

      extern void initialize(const int NX, const int NY, const int NZ,
                             const double rhoAvg,
                             double* ex, double* ey, double* ez, double* wt,
                             double* rho, double* u, double* v, double* w,
                             double* f, double* f_new, double* f_eq);

//    function to stream PDFs to neighboring lattice points

      extern void streaming(const int NX, const int NY, const int NZ,
                            double* ex, double* ey, double* ez, double tau,
                            double* f, double* f_new, double* f_eq);

//    calculate the change in momentum because of inter-particle forces

      extern void calc_dPdt(const int NX, const int NY, const double NZ,
                            double* ex, double* ey, double* ez, double* G11,
                            double* rho, double* dPdt_x, double* dPdt_y, double* dPdt_z);

//    calculate the density and velocity at all nodes

      extern void updateMacro(const int NX, const int NY, const int NZ,
                              double* ex, double* ey, double* ez, double* wt,
                              double tau,
                              double* rho, double* u, double* v, double* w,
                              double* dPdt_x, double* dPdt_y, double* dPdt_z,
                              double* f);

//    update equilibrium PDFs based on the latest {rho,u,v,w}

      extern void updateEquilibrium(const int NX, const int NY, const int NZ,
                                    double* ex, double* ey, double* ez, double* wt,
                                    const double* rho, 
                                    const double* u, const double* v, const double* w,
                                    double* f_eq);

//    writes data to output files using XDMF + HDF5 format

      extern void writeMesh(const int      NX, 
                            const int      NY, 
                            const int      NZ, 
                            const int      time,
                            const double*  rho);

//    MPI 

      int numprocs;          // total number of processors
      int myid;              // processor id
      const int ndims = 3;   // number of dimensions of the Cartesian space
      int dims[ndims];       // number of domain partitions along X, Y and Z
      MPI_Comm CART_COMM;    // new Cartesian communicator after automatic 3D domain decomposition
      int coords[ndims];     // 3D coordinates of process "myid" after domain decomposition
      int nbr_WEST;          // id of neighbor in location (i-1)
      int nbr_EAST;          // id of neighbor in location (i+1)
      int nbr_SOUTH;         // id of neighbor in location (j-1)
      int nbr_NORTH;         // id of neighbor in location (j+1)
      int nbr_BOTTOM;        // id of neighbor in location (k-1)
      int nbr_TOP;           // id of neighbor in location (k+1)

//    global size

      const int NX = 64;         // number of lattice points along X
      const int NY = 64;         // number of lattice points along Y
      const int NZ = 64;         // number of lattice points along Z

      const double delta = 1.0;  // grid spacing is unity along X and Y

      const double x_min = 0;    // global minimum X coordinate
      const double x_max = NX-1;
      const double y_min = 0;    // global minimum Y coordinate
      const double y_max = NY-1;
      const double z_min = 0;    // global minimum Z coordinate
      const double z_max = NZ-1;

//    local parameters (obtained after MPI domain decomposition)

      int LX;   // number of lattice points along X
      int LY;   // number of lattice points along Y
      int LZ;   // number of lattice points along Z

      double local_origin_x;  // X coordinate of the local domain origin
      double local_origin_y;  // Y coordinate of the local domain origin
      double local_origin_z;  // Z coordinate of the local domain origin

      node_range x_range;
      node_range y_range;
      node_range z_range;

//    LBM parameters

      const double GEE11 = -0.27;   // interaction strength
      const double tau = 1.0;       // relaxation time
      const double rhoAvg = 0.693;  // reference density value

//    D2Q9 directions

//                    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
//                    |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
      double ex[] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0};
      double ey[] = { 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
      double ez[] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};

//    weight factors for the various directions

      double wt[] = {1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
                            1./36., 1./36., 1./36., 1./36., 1./36., 1./36.,
                            1./36., 1./36., 1./36., 1./36., 1./36., 1./36., };

//    cohesive force along various lattice directions

      double G11[] = {0, GEE11, GEE11, GEE11, GEE11, GEE11, GEE11,
                         GEE11/2, GEE11/2, GEE11/2, GEE11/2,
                         GEE11/2, GEE11/2, GEE11/2, GEE11/2,
                         GEE11/2, GEE11/2, GEE11/2, GEE11/2};

#endif
