#ifdef __cplusplus
extern "C" {
#endif

#ifndef LINEAL_H_INCLUDED
#define LINEAL_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrices.h"
#include <stdbool.h>
//#include <complex.h>
#include <tgmath.h>

// Rutinas para cálculo lineal

void anguloLineal(char *conjugado_corto, char *conjugado_largo, long double L, long double n, long double L1, long double L2, long double f1, long double f2,long double * Angulos);
long double distanciaCristal(char *tipo_conjugado, long double f, long double Longitud, long double L, long double n, long double theta);
long double spot_lineal(matriz * ABCD, long double lambda);
void calculoLineal(char *conjugado_corto, char *conjugado_largo, long double n, long double lambda, long double L1,\
                   long double L2,long double L,long double f1,long double f2, \
                   _Bool hacerEM1, _Bool hacerEM2, char *rutaSpotEM1, char *rutaSpotEM2,\
                   matriz *wtEM1, matriz *wsEM1, matriz *wtEM2, matriz *wsEM2, matriz *epsilon1, matriz *epsilon2,matriz *qtEM1, matriz *qsEM1, matriz *qtEM2, matriz *qsEM2);

#endif // LINEAL_H_INCLUDED

#ifdef __cplusplus
}
#endif
