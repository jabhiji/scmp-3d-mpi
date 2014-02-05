//    calculate the density and velocity at all nodes

      #include "updateMacro.h"
              
      void updateMacro(const int nn, const int NX, const int NY, const int NZ,
                       double* ex, double* ey, double* ez, double* wt,
                       double tau,
                       double* rho, double* u, double* v, double* w,
                       double* dPdt_x, double* dPdt_y, double* dPdt_z,
                       double* f)
      { 
        const int GX = nn + NX + nn;
        const int GY = nn + NY + nn;

        // update density and velocity
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
              double f_sum = 0;
              double fex_sum = 0;
              double fey_sum = 0;
              double fez_sum = 0;
              for(int id = 0; id < 19; id++)
              {
                int f_index = 19*N + id;
                f_sum   += f[f_index];
                fex_sum += f[f_index]*ex[id];
                fey_sum += f[f_index]*ey[id];
                fez_sum += f[f_index]*ez[id];
              }
              rho[N] = f_sum;
              u[N] = fex_sum / rho[N] + tau * dPdt_x[N] / rho[N];
              v[N] = fey_sum / rho[N] + tau * dPdt_y[N] / rho[N];
              w[N] = fez_sum / rho[N] + tau * dPdt_z[N] / rho[N];
            }
          }
        }
      }
