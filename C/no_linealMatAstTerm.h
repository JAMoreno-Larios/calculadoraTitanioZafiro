#ifdef __cplusplus
extern "C" {
#endif

#ifndef NO_LINEALMATASTTERM_H_INCLUDED
#define NO_LINEALMATASTTERM_H_INCLUDED
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "matrices.h"
#include "lineal.h"
#include "termico.h"
#include "no_lineal.h"
#include "ctocppqprogressbar.h"

typedef struct
{
    long double _Complex qTan;
    long double _Complex qSag;
}propAstigmatica;

// Original
propAstigmatica * propagacionKerrAcoplada(int pasos, long double L, long double n0,long double n2, long double _Complex qInTan, \
                                            long double _Complex qInSag, long double chi, long double kth, long double Cp, \
                                            long double rho, long double dn_dv,long double P_laser, long double lambda0, \
                                            ajusteTemperaturaCristal *vectorPlano, _Bool ladoBombeo);


matriz ** propNoLinealEM1Astigmatico(matriz *qInTan, matriz *qInSag, ajusteTemperaturaCristal **vectorDeVectores, \
                                    char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, \
                                    long double L, long double f1, long double f2, long double n0, long double n2, \
                                    long double chi, long double kth, long double Cp, \
                                    long double rho,long double dn_dv,long double P_laser,\
                                    long double lambda0, matriz *epsilon1, matriz *epsilon2, int iteraciones, int pasosKerr, long double umbral, \
                                    _Bool guardaSpotIteracion, _Bool guardaVariacionIteracion, char* rutaSpotIterMax,\
                                    char* rutaVarIterTan, char* rutaVarIterSag, CtoCppQProgressBar *bar);

matriz ** propNoLinealEM2Astigmatico(matriz *qInTan, matriz *qInSag, ajusteTemperaturaCristal **vectorDeVectores, \
                                    char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, \
                                    long double L, long double f1, long double f2, long double n0, long double n2, \
                                    long double chi, long double kth, long double Cp, \
                                    long double rho,long double dn_dv,long double P_laser,\
                                    long double lambda0, matriz *epsilon1, matriz *epsilon2, int iteraciones, int pasosKerr, long double umbral, \
                                    _Bool guardaSpotIteracion, _Bool guardaVariacionIteracion, char* rutaSpotIterMax,\
                                    char* rutaVarIterTan, char* rutaVarIterSag, CtoCppQProgressBar *bar);

#endif // NO_LINEALMATASTTERM_H_INCLUDED

#ifdef __cplusplus
}
#endif
