//    function to stream PDFs to neighboring lattice points

      #include "streaming.h"

      void streaming(const int NX, const int NY, const int NZ,
                     double* ex, double* ey, double* ez, double tau,
                     double* f, double* f_new, double* f_eq)
      {
        for(int k = 0; k < NZ-1; k++)
        {
          for(int j = 0; j < NY-1; j++)
          {
            for(int i = 0; i < NX-1; i++)
            {
              int N = i + NX*j + NX*NY*k;
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
                int f_index_beg = 19*N + id; 
                int f_index_end = 19*Nflow + id;
        
                f_new[f_index_end] = f[f_index_beg]
                                   - (f[f_index_beg] - f_eq[f_index_beg])
                                   / tau;

            //  printf("N = %d Nflow = %d f = %f f_eq = %f \n",N, Nflow, f[f_index_beg], f_eq[f_index_beg]);

              }
            }
          }
        }
      }
