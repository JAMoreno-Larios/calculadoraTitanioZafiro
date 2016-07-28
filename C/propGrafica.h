#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "matrices.h"
#include "lineal.h"
#include "no_lineal.h"
#include "termico.h"


/* Estructura de vector */

typedef struct
{
    int numeroElementos;
    long double posKerrIn; // Posición donde el haz entra al cristal
    long double posKerrOut; // Posición donde el haz sale del cristal
    long double posKerrMid; // Sección media del cristal.
    long double * spotCW; // Radio del haz
    long double * spotML; // Radio del haz
    long double * pos; // Posición en cavidad.
}vectorDato;

void escribeData(vectorDato *datos, char* filenameCW, char* filenameML);

/* Escribe a gnuplot - Todos los datos. */
void gnuplotEscribe(vectorDato *EM1tan, vectorDato *EM1sag, vectorDato *EM2tan, vectorDato *EM2sag);


/* Rutina para gráficar propagación lineal - Conlleva un error por la discreti-
zación de la propagación de rayos (en lugar de una matriz para propagar en espacio,
se utilizan muchas matrices */

void propagacionKerrGrafica(int pasos, long double deltaZ, long double n0,long double n2,
                             long double complex qInTan, long double complex qInSag, long double complex *qTan, \
                                            long double complex *qSag, long double *wTan, long double *wSag, long double chi, long double kth, long double Cp, \
                                            long double rho, long double dn_dv,long double P_laser, long double lambda0, \
                                            ajusteTemperaturaCristal *vectorPlano, _Bool ladoBombeo);
void graficaPropagacion(char *conjugado_corto, char *conjugado_largo, int iteraciones, long double umbral, int N, long double Epsilon1, long double Epsilon2,
                     long double lambdaPump, long double n0, long double n2, long double P_laser, long double P_pump, long double chi, long double L,
                     long double kth, long double Cp, long double rho, long double dn_dv, long double ancho, long double alto, long double w_pump_t,
                     long double w_pump_s, long double nPump, long double lambda0, long double L1, long double L2,
                     long double f1, long double f2);
