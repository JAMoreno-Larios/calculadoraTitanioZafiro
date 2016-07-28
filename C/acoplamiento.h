#ifdef __cplusplus
extern "C" {
#endif

#ifndef ACOPLAMIENTO_H_INCLUDED
#define ACOPLAMIENTO_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "matrices.h"
#include "lineal.h"
#include "no_lineal.h"
#include "ctocppqprogressbar.h" // Agrega barra de control

typedef struct
{
    long double _Complex *qOutTan;
    long double _Complex *qOutSag;
    int numSpots;
}spotsPump;

spotsPump * nuevoSpotsPump(int numSpots);
spotsPump * borraSpotsPump( spotsPump * sP);
spotsPump * acopleOptico(long double wFuenteLaserTan, long double wFuenteLaserSag, long double divergencia, long double lambda0,\
                         long double La, long double tL, long double R11, long double R12,\
                         long double nL, long double thetaL, long double Lb, long double tE,\
                         long double R21, long double R22, long double nE, long double thetaE,\
                          long double nC, long double delta1, matriz * epsilon1);
spotsPump * acopleOpticoBarra(long double wFuenteLaserTan, long double wFuenteLaserSag, long double divergencia, long double lambda0,\
                         long double La, long double tL, long double R11, long double R12,\
                         long double nL, long double thetaL, long double Lb, long double tE,\
                         long double R21, long double R22, long double nE, long double thetaE,\
                          long double nC, long double delta1, matriz * epsilon1, CtoCppQProgressBar *barra);
#endif

#ifdef __cplusplus
}
#endif
