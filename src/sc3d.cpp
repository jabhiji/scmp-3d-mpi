//    ------------------------------------------------------------------
//    Simulation of Single Component Multiphase flow in 2D
//
//    Shan & Chen Model
//
//    Periodic boundary conditions
//
//    Written by: Abhijit Joshi
//    ------------------------------------------------------------------

      #include "sc3d.h"

      int main(int argc, char* argv[])
      {
//      set up MPI and implement Cartesian domain decomposition
//      identify coordinates and neighboring MPI ranks

        mpiSetup(argc, argv, &numprocs, &myid, ndims,
                 &dims[0], &coords[0], CART_COMM,
                 &nbr_WEST, &nbr_EAST,
                 &nbr_SOUTH, &nbr_NORTH,
                 &nbr_BOTTOM, &nbr_TOP);

//      calculate size of local 3D sub-domain handled by this rank

        domainDecomp3D(myid, CART_COMM, dims, coords,
                       NX,
                       NY,
                       NZ,
                       delta,
                       x_min,
                       y_min,
                       z_min,
                       x_range,
                       y_range,
                       z_range,
                       local_origin_x,
                       local_origin_y,
                       local_origin_z,
                       LX,
                       LY,
                       LZ);

//      define buffers

        double *rho    = new double[NX*NY*NZ]; // density
        double *u      = new double[NX*NY*NZ]; // velocity x-component
        double *v      = new double[NX*NY*NZ]; // velocity y-component
        double *w      = new double[NX*NY*NZ]; // velocity z-component
        double *dPdt_x = new double[NX*NY*NZ]; // momentum change along x
        double *dPdt_y = new double[NX*NY*NZ]; // momentum change along y
        double *dPdt_z = new double[NX*NY*NZ]; // momentum change along z
        double *f      = new double[NX*NY*NZ*19]; // PDF
        double *f_eq   = new double[NX*NY*NZ*19]; // PDF
        double *f_new  = new double[NX*NY*NZ*19]; // PDF

//      initialize fields

        initialize(NX, NY, NZ, rhoAvg, &ex[0], &ey[0], &ez[0], &wt[0], 
                   rho, u, v, w, f, f_new, f_eq);

//      time integration

        int time = 0;
        clock_t t0, tN;
        t0 = clock();

//      write initial condition to output files

        writeMesh(NX, NY, NZ, time, rho);

//      time integration loop

        while(time < 1000)
        {
          time++; // increment lattice time

          streaming(NX, NY, NZ, ex, ey, ez, tau, f, f_new, f_eq);

          calc_dPdt(NX, NY, NZ, ex, ey, ez, G11, rho, dPdt_x, dPdt_y, dPdt_z);

          updateMacro(NX, NY, NZ, ex, ey, ez, wt, tau, 
                      rho, u, v, w, dPdt_x, dPdt_y, dPdt_z, f);

          updateEquilibrium(NX, NY, NZ, ex, ey, ez, wt, rho, u, v, w, f_eq);

//        transfer fnew back to f

          for(int f_index = 0; f_index < NX*NY*NZ*19; f_index++)
          {
            f[f_index] = f_new[f_index];
          }

//        write output data using (XDMF+HDF5)

          if(time%100 == 0) 
          {
            writeMesh(NX, NY, NZ, time, rho);
          }

//        calculate the number of lattice time-steps per second

          tN = clock() - t0;

//        std::cout << " lattice time steps per second = " 
//                  << (float) CLOCKS_PER_SEC * time / (float) tN 
//                  << std::endl;
        }

//      clean up

        delete[] rho;
        delete[] u;
        delete[] v;
        delete[] w;
        delete[] dPdt_x;
        delete[] dPdt_y;
        delete[] dPdt_z;
        delete[] f;
        delete[] f_eq;
        delete[] f_new;

//      main program ends

        return 0;
      }
