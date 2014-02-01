//    function to initialize density, velocity and PDFs

      #include "initialize.h"

      void initialize(const int nn, const int NX, const int NY, const int NZ, const int myid,
                      const double rhoAvg,
                      double* ex, double* ey, double* ez, double* wt,
                      double* rho, double* u, double* v, double* w,
                      double* f, double* f_new, double* f_eq)
      {
        std::cout << "Initializing buffers.....";

//      initialize random seed

        srand (time(NULL) + myid);

//      initialize density and velocity

        const int GX = nn + NX + nn;
        const int GY = nn + NY + nn;
        const int GZ = nn + NZ + nn;

        double rhoVar = 0.01 * rhoAvg;
        for(int k = 0; k < NZ; k++)
        {
          int K = nn+k;
          for(int j = 0; j < NY; j++)
          {
            int J = nn+j;
            for(int i = 0; i < NX; i++)
            {
              int I = nn+i;
              int N = I + GX*J + GX*GY*K;
              rho[N] = rhoAvg - 0.5*rhoVar + rhoVar * rand()/RAND_MAX;
              u[N] = 0.0;
              v[N] = 0.0;
              w[N] = 0.0;
            }
          }
        }

//      initialize distribution functions to their equilibrium value

        for(int k = 0; k < NZ; k++)
        {
          int K = nn+k;
          for(int j = 0; j < NY; j++)
          {
            int J = nn+j;
            for(int i = 0; i < NX; i++)
            {
              int I = nn+i;
              int N = I + GX*J + GX*GY*K;
              double udotu = u[N]*u[N] + v[N]*v[N] + w[N]*w[N];

              for(int id = 0; id < 19; id++)
              {
                int index_f = 19*N + id;
                double edotu = ex[id]*u[N] + ey[id]*v[N] + ez[id]*w[N];
                f_eq[index_f] = wt[id] * rho[N]
                              * (1 + 3*edotu
                                   + 4.5*edotu*edotu - 1.5*udotu);
                f[index_f] = f_eq[index_f];
                f_new[index_f] = f_eq[index_f];
              }
            }
          }
        }

        std::cout << "Done\n";
      }
