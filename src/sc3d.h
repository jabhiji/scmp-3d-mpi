#ifndef SC3D_H
#define SC3D_H

//    C++ headers

      #include <iostream>     // cout()
      #include <cmath>        // pow()
      #include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC

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

//    lattice size

      const int NX = 64;         // number of lattice points along X
      const int NY = 64;         // number of lattice points along Y
      const int NZ = 64;         // number of lattice points along Z

//    domain size in lattice units
//    grid spacing is unity along X and Y

      const double xmin = 0;
      const double xmax = NX-1;
      const double ymin = 0;
      const double ymax = NY-1;
      const double zmin = 0;
      const double zmax = NZ-1;

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
