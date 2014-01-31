//    calculate the density and velocity at all nodes

      #include "updateMacro.h"
              
      void updateMacro(const int NX, const int NY, const int NZ,
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

