//    calculate the change in momentum because of inter-particle forces

      #include "calc_dPdt.h"

//    funtion to calculate effective density in the Shan & Chen model

      double psi(double x)
      {
        const double E = 2.71828;
        const double rho0 = 1.0;
        return rho0 * (1 - pow(E, -x/rho0));
      }

      void calc_dPdt(const int NX, const int NY, const double NZ,
                     double* ex, double* ey, double* ez, double* G11,
                     double* rho, double* dPdt_x, double* dPdt_y, double* dPdt_z)
      { 
        // interparticle forces
        for(int k = 0; k < NZ-1; k++)
        {  
          for(int j = 0; j < NY-1; j++)
          {  
            for(int i = 0; i < NX-1; i++)
            {  
              int N = i + NX*j + NX*NY*k;
              double Gsumx = 0.;
              double Gsumy = 0.;
              double Gsumz = 0.;
              for(int id = 0; id < 19; id++)
              {
                int iflow = i + ex[id];
                int jflow = j + ey[id];
                int kflow = k + ez[id];

                // periodic B.C.
                if(iflow == -1) iflow = NX-2;
                if(jflow == -1) jflow = NY-2;
                if(kflow == -1) kflow = NZ-2;
                if(iflow == NX-1) iflow = 0;
                if(jflow == NY-1) jflow = 0;
                if(kflow == NZ-1) kflow = 0;

                int Nflow = iflow + NX*jflow + NX*NY*kflow;

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
