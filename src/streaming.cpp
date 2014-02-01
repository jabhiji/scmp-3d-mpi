//    function to stream PDFs to neighboring lattice points

      #include "streaming.h"

      void streaming(const int nn, const int NX, const int NY, const int NZ,
                     double* ex, double* ey, double* ez, double tau,
                     double* f, double* f_new, double* f_eq)
      {

        const int GX = nn + NX + nn;  // size along X including ghost nodes
        const int GY = nn + NY + nn;  // size along Y including ghost nodes

        // stream TO all interior nodes

        for(int k = 0; k < NZ; k++)
        {
          int K = nn + k;

          for(int j = 0; j < NY; j++)
          {
            int J = nn + j;

            for(int i = 0; i < NX; i++)
            {
              int I = nn + i;

              int N = I + GX*J + GX*GY*K;  // streaming destination

              for(int id = 0; id < 19; id++)
              {
                int ifrom = I - ex[id];
                int jfrom = J - ey[id];
                int kfrom = K - ez[id];
       
                int Nfrom = ifrom + GX*jfrom + GX*GY*kfrom;
                int f_index_end = 19*N + id; 
                int f_index_beg = 19*Nfrom + id;
        
                f_new[f_index_end] = f[f_index_beg]
                                   - (f[f_index_beg] - f_eq[f_index_beg])
                                   / tau;
              }
            }
          }
        }
      }
