//    ------------------------------------------------------------------
//    Simulation of Single Component Multiphase flow in 2D
//
//    Shan & Chen Model
//
//    Periodic boundary conditions
//
//    Written by: Abhijit Joshi
//    ------------------------------------------------------------------

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

//    writes data to output files using XDMF + HDF5 format

      extern void writeMesh(const int      NX, 
                            const int      NY, 
                            const int      NZ, 
                            const int      time,
                            const double*  rho);

//    calculate the density and velocity at all nodes

      void updateDensityAndVelocity(const int NX, const int NY, const int NZ,
                                    double* ex, double* ey, double* ez, double* wt,
                                    double tau,
                                    double* rho, double* u, double* v, double* w,
                                    double* dPdt_x, double* dPdt_y, double* dPdt_z,
                                    double* f)
      {
        // update density and velocity
        for(int k = 0; k < NZ-1; k++)
        {  
          for(int j = 0; j < NY-1; j++)
          {  
            for(int i = 0; i < NX-1; i++)
            {
              int N = i + NX*j + NX*NY*k;
              double f_sum = 0;
              double fex_sum = 0;
              double fey_sum = 0;
              double fez_sum = 0;
              for(int id = 0; id < 19; id++)
              {  
                f_sum += f[19*N + id];
                fex_sum += f[19*N + id]*ex[id];
                fey_sum += f[19*N + id]*ey[id];
                fez_sum += f[19*N + id]*ez[id];
              }
              rho[N] = f_sum;
              u[N] = fex_sum / rho[N] + tau * dPdt_x[N] / rho[N];
              v[N] = fey_sum / rho[N] + tau * dPdt_y[N] / rho[N];
              w[N] = fez_sum / rho[N] + tau * dPdt_z[N] / rho[N];
            }
          }
        }

        // periodic B.C. for rho on top face
        for(int i = 0; i < NX-1; i++)
        {
          for(int j = 0; j < NY-1; j++)
          {
            int k;
            k = NZ-1; int N_end = i + NX*j + NX*NY*k;
            k = 0   ; int N_beg = i + NX*j + NX*NY*k;
            rho[N_end] = rho[N_beg];
          }
        }

        // periodic B.C. for rho on north face
        for(int i = 0; i < NX-1; i++)
        {
          for(int k = 0; k < NZ-1; k++)
          {
            int j;
            j = NY-1; int N_end = i + NX*j + NX*NY*k;
            j = 0   ; int N_beg = i + NX*j + NX*NY*k;
            rho[N_end] = rho[N_beg];
          }
        }

        // periodic B.C. for rho on east face
        for(int j = 0; j < NY-1; j++)
        {
          for(int k = 0; k < NZ-1; k++)
          {
            int i;
            i = NX-1; int N_end = i + NX*j + NX*NY*k;
            i = 0   ; int N_beg = i + NX*j + NX*NY*k;
            rho[N_end] = rho[N_beg];
          }
        }

        // periodic B.C. for edge y = NY-1, z = NZ-1
        for(int i = 0; i < NX-1; i++)
        {
          int j = NY-1, k = NZ-1;
          int N_end = i + NX*j + NX*NY*k;
          int N_beg = i + NX*0 + NX*NY*0;
          rho[N_end] = rho[N_beg];
        }

        // periodic B.C. for edge z = NZ-1, x = NX-1
        for(int j = 0; j < NY-1; j++)
        {
          int k = NZ-1, i = NX-1;
          int N_end = i + NX*j + NX*NY*k;
          int N_beg = 0 + NX*j + NX*NY*0;
          rho[N_end] = rho[N_beg];
        }

        // periodic B.C. for edge x = NX-1, y = NY-1
        for(int k = 0; k < NZ-1; k++)
        {
          int i = NX-1, j = NY-1;
          int N_end = i + NX*j + NX*NY*k;
          int N_beg = 0 + NX*0 + NX*NY*k;
          rho[N_end] = rho[N_beg];
        }

        // periodic B.C. for corner
        rho[NX*NY*NZ-1] = rho[0];
      }

//    update equilibrium PDFs based on the latest {rho,u,v,w}

      void updateEquilibrium(const int NX, const int NY, const int NZ,
                             double* ex, double* ey, double* ez, double* wt,
                             const double* rho, 
                             const double* u, const double* v, const double* w,
                             double* f_eq)
      {
        for(int k = 0; k < NZ-1; k++)
        {  
          for(int j = 0; j < NY-1; j++)
          {  
            for(int i = 0; i < NX-1; i++)
            {
              int N = i + NX*j + NX*NY*k;
              double udotu = u[N]*u[N] + v[N]*v[N] + w[N]*w[N];
              for(int id = 0; id < 19; id++)
              {
                int index_f = 19*N + id;
                double edotu = ex[id]*u[N] + ey[id]*v[N] + ez[id]*w[N];
                f_eq[index_f] = wt[id] * rho[N] 
                              * (1 + 3*edotu
                                   + 4.5*edotu*edotu - 1.5*udotu);
              }
            }
          }
        }
      }

//    main program

      int main(void)
      {
//      lattice size

        const int NX = 64;         // number of lattice points along X
        const int NY = 64;         // number of lattice points along Y
        const int NZ = 64;         // number of lattice points along Z

        // domain size in lattice units
        // grid spacing is unity along X and Y

        const double xmin = 0;
        const double xmax = NX-1;
        const double ymin = 0;
        const double ymax = NY-1;
        const double zmin = 0;
        const double zmax = NZ-1;

//      LBM parameters

        const double GEE11 = -0.27;   // interaction strength
        const double tau = 1.0;       // relaxation time
        const double rhoAvg = 0.693;  // reference density value

//      D2Q9 directions

//                      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
//                      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
        double ex[] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0};
        double ey[] = { 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
        double ez[] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};

//      weight factors for the various directions

        double wt[] = {1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
                              1./36., 1./36., 1./36., 1./36., 1./36., 1./36.,
                              1./36., 1./36., 1./36., 1./36., 1./36., 1./36., };

//      cohesive force along various lattice directions

        double G11[] = {0, GEE11, GEE11, GEE11, GEE11, GEE11, GEE11,
                           GEE11/2, GEE11/2, GEE11/2, GEE11/2,
                           GEE11/2, GEE11/2, GEE11/2, GEE11/2,
                           GEE11/2, GEE11/2, GEE11/2, GEE11/2};

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

        initialize(NX, NY, NZ, rhoAvg, 
                   &ex[0], &ey[0], &ez[0], &wt[0], 
                   rho, u, v, w, f, f_new, f_eq);

//      time integration

        int time = 0;
        clock_t t0, tN;
        t0 = clock();

//      write initial condition to output files

        writeMesh(NX, NY, NZ, time, rho);

//      time integration loop

        while(time < 100)
        {
          time++; // increment lattice time

          streaming(NX, NY, NZ, ex, ey, ez, tau, f, f_new, f_eq);

          calc_dPdt(NX, NY, NZ, ex, ey, ez, G11, rho, dPdt_x, dPdt_y, dPdt_z);

          updateDensityAndVelocity(NX, NY, NZ, ex, ey, ez, wt, tau, 
                                   rho, u, v, w, dPdt_x, dPdt_y, dPdt_z, f);

          updateEquilibrium(NX, NY, NZ, ex, ey, ez, wt, rho, u, v, w, f_eq);

//        transfer fnew back to f

          for(int f_index = 0; f_index < NX*NY*NZ*19; f_index++)
          {
            f[f_index] = f_new[f_index];
          }

//        write output data using (XDMF+HDF5)

          if(time%10 == 0) 
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
