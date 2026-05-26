#ifndef INTEGRATORS
#define INTEGRATORS


int leapfrog_kdk(double time_step);
int yoshida_fourth_order(double time_step);
int yoshida_eighth_order(double time_step);

#endif