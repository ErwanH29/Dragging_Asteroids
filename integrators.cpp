#include <cmath>
#include "integrators.h"


extern int kick(double time_step);
extern int drift(double time_step);
extern int get_acc_pot_coll();



int yoshida_fourth_order(double time_step){
    // Coefficients from 1990PhLA..150..262Y
    static const double w0 = -cbrt(2.0) / (2.0 - cbrt(2.0));
    static const double w1 = 1.0 / (2.0 - cbrt(2.0));
    static const double c1 = w1 / 2.0;
    static const double c2 = (w0 + w1) / (2.0);
    static const double c3 = (w0 + w1) / (2.0);
    static const double c4 = w1 / 2.0;
    static const double d1 = w1;
    static const double d2 = w0;
    static const double d3 = w1;

    int nparticles = radius.size();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c1 * time_step * vel[i][k];
        }
    }

    get_acc_pot_coll();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += d1 * time_step * acc[i][k];
        }
    }
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c2 * time_step * vel[i][k];
        }
    }

    get_acc_pot_coll();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += d2 * time_step * acc[i][k];
        }
    }
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c3 * time_step * vel[i][k];
        }
    }

    get_acc_pot_coll();
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            vel[i][k] += d3 * time_step * acc[i][k];
        }
    }
    for (int i=0; i<nparticles; i++){
        for (int k=0; k<NDIM; k++){
            pos[i][k] += c4 * time_step * vel[i][k];
        }
    }

    return 0;
}



int leapfrog_kdk(double time_step){
    kick(0.5 * time_step);
    drift(time_step);
    get_acc_pot_coll();
    kick(0.5 * time_step);
    return 0;
}
