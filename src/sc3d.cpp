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
                       LX,                // local nodes along X
                       LY,                // local nodes along Y
                       LZ);               // local nodes along Z

//      ghost layer thickness

        const int nn = 1;

//      define local buffers for this MPI rank

        const int size1 = (nn+LX+nn) * (nn+LY+nn) * (nn+LZ+nn);
        const int size2 = size1 * 19;

        double *rho    = new double[size1]; // density
        double *u      = new double[size1]; // velocity x-component
        double *v      = new double[size1]; // velocity y-component
        double *w      = new double[size1]; // velocity z-component
        double *dPdt_x = new double[size1]; // momentum change along x
        double *dPdt_y = new double[size1]; // momentum change along y
        double *dPdt_z = new double[size1]; // momentum change along z

        double *f      = new double[size2]; // PDF
        double *f_eq   = new double[size2]; // PDF
        double *f_new  = new double[size2]; // PDF

//      initialize fields

        initialize(nn, LX, LY, LZ, myid,
                   local_origin_x, local_origin_y, local_origin_z,
                   rhoAvg, &ex[0], &ey[0], &ez[0], &wt[0], 
                   rho, u, v, w, f, f_new, f_eq);

        // fill ghost layers in the macroscopic variable buffers ( rho, u, v, w )

        fillGhostLayersMacVar(nn,              // ghost layer thickness
                              LX,              // number of nodes along X (local for this MPI process)
                              LY,              // number of nodes along Y (local for this MPI process)
                              LZ,              // number of nodes along Z (local for this MPI process)
                              myid,            // MPI process id or rank
                              CART_COMM,       // Cartesian communicator
                              nbr_WEST,        // neighboring MPI process to my west
                              nbr_EAST,        // neighboring MPI process to my east
                              nbr_SOUTH,       // neighboring MPI process to my south
                              nbr_NORTH,       // neighboring MPI process to my north
                              nbr_BOTTOM,      // neighboring MPI process to my bottom
                              nbr_TOP,         // neighboring MPI process to my top
                              rho,            // density
                              u,              // velocity (x-component)
                              v,              // velocity (y-component)
                              w);             // velocity (z-component)

        exchangePDF (nn,                // number of ghost cell layers
                     Q,                 // number of LBM streaming directions
                     LX,                // number of voxels along X in this process
                     LY,                // number of voxels along Y in this process
                     LZ,                // number of voxels along Z in this process
                     myid,              // my process id
                     CART_COMM,         // Cartesian topology communicator
                     nbr_WEST,          // process id of my western neighbor
                     nbr_EAST,          // process id of my eastern neighbor
                     nbr_SOUTH,         // process id of my southern neighbor
                     nbr_NORTH,         // process id of my northern neighbor
                     nbr_BOTTOM,        // process id of my bottom neighbor
                     nbr_TOP,           // process id of my top neighbor
                     f);                // pointer to the 4D array being exchanged (of type double)

        exchangePDF (nn,                // number of ghost cell layers
                     Q,                 // number of LBM streaming directions
                     LX,                // number of voxels along X in this process
                     LY,                // number of voxels along Y in this process
                     LZ,                // number of voxels along Z in this process
                     myid,              // my process id
                     CART_COMM,         // Cartesian topology communicator
                     nbr_WEST,          // process id of my western neighbor
                     nbr_EAST,          // process id of my eastern neighbor
                     nbr_SOUTH,         // process id of my southern neighbor
                     nbr_NORTH,         // process id of my northern neighbor
                     nbr_BOTTOM,        // process id of my bottom neighbor
                     nbr_TOP,           // process id of my top neighbor
                     f_new);            // pointer to the 4D array being exchanged (of type double)

        exchangePDF (nn,                // number of ghost cell layers
                     Q,                 // number of LBM streaming directions
                     LX,                // number of voxels along X in this process
                     LY,                // number of voxels along Y in this process
                     LZ,                // number of voxels along Z in this process
                     myid,              // my process id
                     CART_COMM,         // Cartesian topology communicator
                     nbr_WEST,          // process id of my western neighbor
                     nbr_EAST,          // process id of my eastern neighbor
                     nbr_SOUTH,         // process id of my southern neighbor
                     nbr_NORTH,         // process id of my northern neighbor
                     nbr_BOTTOM,        // process id of my bottom neighbor
                     nbr_TOP,           // process id of my top neighbor
                     f_eq);             // pointer to the 4D array being exchanged (of type double)

//      time integration

        int time = 0;
        clock_t t0, tN;
        t0 = clock();

//      write initial condition to output files

        writeMesh(nn, CART_COMM, myid, 
                  local_origin_x, local_origin_y, local_origin_z, delta, 
                  LX, LY, LZ, time, rho);

//      time integration loop

        while(time < MAXIMUM_TIME)
        {
          time++; // increment lattice time

          streaming(nn, LX, LY, LZ, ex, ey, ez, tau, f, f_new, f_eq);

          calc_dPdt(nn, LX, LY, LZ, ex, ey, ez, G11, rho, dPdt_x, dPdt_y, dPdt_z);

          updateMacro(nn, LX, LY, LZ, ex, ey, ez, wt, tau, 
                      rho, u, v, w, dPdt_x, dPdt_y, dPdt_z, f);

          // fill ghost layers in the macroscopic variable buffers ( rho, u, v, w )

          fillGhostLayersMacVar(nn,              // ghost layer thickness
                                LX,              // number of nodes along X (local for this MPI process)
                                LY,              // number of nodes along Y (local for this MPI process)
                                LZ,              // number of nodes along Z (local for this MPI process)
                                myid,            // MPI process id or rank
                                CART_COMM,       // Cartesian communicator
                                nbr_WEST,        // neighboring MPI process to my west
                                nbr_EAST,        // neighboring MPI process to my east
                                nbr_SOUTH,       // neighboring MPI process to my south
                                nbr_NORTH,       // neighboring MPI process to my north
                                nbr_BOTTOM,      // neighboring MPI process to my bottom
                                nbr_TOP,         // neighboring MPI process to my top
                                rho,            // density
                                u,              // velocity (x-component)
                                v,              // velocity (y-component)
                                w);             // velocity (z-component)

          updateEquilibrium(nn, LX, LY, LZ, ex, ey, ez, wt, rho, u, v, w, f_eq);

          exchangePDF (nn,                // number of ghost cell layers
                       Q,                 // number of LBM streaming directions
                       LX,                // number of voxels along X in this process
                       LY,                // number of voxels along Y in this process
                       LZ,                // number of voxels along Z in this process
                       myid,              // my process id
                       CART_COMM,         // Cartesian topology communicator
                       nbr_WEST,          // process id of my western neighbor
                       nbr_EAST,          // process id of my eastern neighbor
                       nbr_SOUTH,         // process id of my southern neighbor
                       nbr_NORTH,         // process id of my northern neighbor
                       nbr_BOTTOM,        // process id of my bottom neighbor
                       nbr_TOP,           // process id of my top neighbor
                       f_eq);             // pointer to the 4D array being exchanged (of type double)

//        transfer fnew back to f

          const int GX = nn + LX + nn;  // size along X including ghost nodes
          const int GY = nn + LY + nn;  // size along Y including ghost nodes
          const int GZ = nn + LZ + nn;  // size along Z including ghost nodes
          for(int f_index = 0; f_index < GX*GY*GZ*19; f_index++)
          {
            f[f_index] = f_new[f_index];
          }

//        write output data using (XDMF+HDF5)

          if(time%frame_rate == 0) 
          {
             writeMesh(nn, CART_COMM, myid, 
                       local_origin_x, local_origin_y, local_origin_z, delta, 
                       LX, LY, LZ, time, rho);
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

//      MPI clean up

        MPI_Finalize();

//      main program ends

        return 0;
      }
