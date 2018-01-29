/******************************************************************************
* The order is P_0(k), P_2(k), P_4(k), B_0(k1,k2,k3), etc. 
******************************************************************************/
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "gsl/gsl_matrix.h"

int main() {

    double V = 1000;
    double kmax = 0.2;
    int Nk = 20;

    double dk = kmax/Nk;

    double klow[Nk];
    double khigh[Nk];
    double kmid[Nk];

    int multipole[3] = [0, 2, 4];

    int size = Nk*3 + Nk*(Nk + 1)*(Nk + 2)*3/2;

    gsl_matrix * PBFisher = gsl_matrix_alloc(size,size);

    for (int i = 0; i < Nk; i++) {
        klow[i] = dk*i;
        khigh[i] = dk*(i + 1);
        kmid[i] = dk*(i + 0.5);
    }

    double cov;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j <= i; j++) {
            // Both power spectra
            if (i < 3*Nk && j < 3*Nk) {
                int il1 = i/Nk;
                int ik1 = i%Nk;
                int il2 = j/Nk;
                int ik2 = j%Nk;
                if (ik1 != ik2) {
                    cov = 0;
                    break;
                }
                else {
                    cov = Pcov(ik1,il1,il2);
                }
            }
            // Power spectrum vs Bispectrum
            else if (i < 3*Nk) {
                int il1 = i/Nk;
                int ik1 = i%Nk;
                int il2 = (i - 3*Nk)/;
                int ik2 = ;
                int ik3 = ;
                int ik4 = ;
                if (ik1 != ik2 && ik1 != ik3 && ik1 != ik4) {
                    cov = 0;
                    break;
                }
                else {
                    cov = PBcross(ik1,il1,ik2,ik3,ik4,il2);
                }
            }
            // Both Bispectra
            else {

            }
            gsl_matrix_set(PBcov, i, j, cov);
            gsl_matrix_set(PBcov, j, i, cov);
        }
    }

}
