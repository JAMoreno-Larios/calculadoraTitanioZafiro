#ifdef __cplusplus
extern "C" {
#endif

#ifndef NO_LINEAL_H_INCLUDED
#define NO_LINEAL_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <complex.h>
#include <tgmath.h>
#include "matrices.h"
#include "lineal.h"

long double _Complex propagacionKerr(int pasos, long double L, long double n0,long double n2, long double w_pump,long double _Complex q_in, long double chi, long double kth, long double Cp, long double rho, long double dn_dv,long double P_pump,long double P_laser, long double lambda0);
long double _Complex * propagacionKerrMatriz(int pasos, long double L, long double n0,long double n2, long double w_pump,long double _Complex q_in, long double chi, long double kth, long double Cp, long double rho, long double dn_dv,long double P_pump,long double P_laser, long double lambda0);
long double _Complex * propagacionLibre(long double L,int pasos, long double n0,long double _Complex q_in);
long double _Complex prop_q(matriz * ABCD, long double _Complex q_in, long double n1, long double n2); // Propagación de un haz gaussiano a través de un sistema óptico ABCD - Siegman
long double spot_q(long double _Complex q, long double n, long double lambda);
long double radio_q(long double _Complex q);
void confNoAstigmatica(matriz * A, long double _Complex *casiUno, int *fila, int *columna);
matriz * spot_q_matriz(matriz * q, long double n, long double lambda);


#endif // NO_LINEAL_H_INCLUDED

#ifdef __cplusplus
}
#endif
