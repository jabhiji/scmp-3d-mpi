//    update equilibrium PDFs based on the latest {rho,u,v,w}

      #include "updateEquilibrium.h"

      void updateEquilibrium(const int nn, const int NX, const int NY, const int NZ,
                             double* ex, double* ey, double* ez, double* wt,
                             const double* rho, 
                             const double* u, const double* v, const double* w,
                             double* f_eq)
      {
        const int GX = nn + NX + nn;
        const int GY = nn + NY + nn;

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
              }
            }
          }
        }
      }
