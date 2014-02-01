//    calculate the change in momentum because of inter-particle forces

      #include "calc_dPdt.h"

//    funtion to calculate effective density in the Shan & Chen model

      double psi(double x)
      {
        const double E = 2.71828;
        const double rho0 = 1.0;
        return rho0 * (1 - pow(E, -x/rho0));
      }

      void calc_dPdt(const int nn, const int NX, const int NY, const double NZ,
                     double* ex, double* ey, double* ez, double* G11,
                     double* rho, double* dPdt_x, double* dPdt_y, double* dPdt_z)
      { 
        const int GX = nn + NX + nn;
        const int GY = nn + NY + nn;
        const int GZ = nn + NZ + nn;
        // interparticle forces
        for(int k = 0; k < NZ; k++)
        {  
          int K = nn + k;
          for(int j = 0; j < NY; j++)
          {  
            int J = nn + j;
            for(int i = 0; i < NX; i++)
            {  
              int I = nn + i;
              int N = I + GX*J + GX*GY*K;
              double Gsumx = 0.;
              double Gsumy = 0.;
              double Gsumz = 0.;
              for(int id = 0; id < 19; id++)
              {
                int iflow = I + ex[id];
                int jflow = J + ey[id];
                int kflow = K + ez[id];

                int Nflow = iflow + GX*jflow + GX*GY*kflow;

                double strength = psi(rho[N]) * psi(rho[Nflow]) * G11[id];

                Gsumx += strength * ex[id];
                Gsumy += strength * ey[id];
                Gsumz += strength * ez[id];
              }
              dPdt_x[N] = -Gsumx;
              dPdt_y[N] = -Gsumy;
              dPdt_z[N] = -Gsumz;

          //  printf("calc_dPdt --- %f %f %f\n", Gsumx, Gsumy, Gsumz);

            }
          }
        }
      }
