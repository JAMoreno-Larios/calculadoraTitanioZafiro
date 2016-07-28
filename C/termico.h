
#ifdef __cplusplus
extern "C" {
#endif

#ifndef TERMICO_H_INCLUDED
#define TERMICO_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "matrices.h"
#include "lineal.h"
#include "no_lineal.h"
#include "acoplamiento.h"
#include "ctocppqprogressbar.h"

// Rutina para resolver un sistema de n ecuaciones algebráicas de n incógnitas
// mediante Gauss-Seidel. La matriz A representa los coeficientes constantes, X
// son las incógnitas y B son los resultados. La ecuación en forma matricial
// es [A]{X}={B} El algoritmo emplea relajación para mejorar la convergencia, esto es, cada
// valor nuevo de x calculado y modificado por un promedio ponderado de los resultados de
// la iteración anterior y actual: xi_nuevo=lambda*xi_nuevo'+(1-lambda)*xi_viejo.
// lambda se determina de forma empírica. Si lambda=0, no hay modificación al sistema.

// Función Gauss-Seidel

long double * gaussSeidel(long double **A, long double *B, int n, int iteraciones, long double umbral, long double lambda);
// Estructura para guardar la información relevante a la distribución de temperatura en un plano
// Almacena ejes cartesianos, valores del plano, y profundiad en el material usada.

typedef struct
{
    int N; // Órden de la matriz a usar.
    int paso; // Número de paso. dz=L/(pasos-1); xi(i)=dz*paso(i). exp(-alpha*xi)
    long double deltaX; // Diferencial en el plano tangencial
    long double deltaY; // Dirección en el plano sagital
    long double * X; // Vector de posición tangencial
    long double * Y; // Vector de posición sagital
    long double * T;   // Mapa de datos
}distTempPlano;

// Estructura para guardar los coeficientes del ajuste parabólico de la distribución de temperatura
// Almacena los coeficientes de T=a0+a1*X^2+a2*Y^2 y la posición correspondiente (Pos=L/xi)
// Los vectores de los coeficientes tendrán un total de pasosTotales pasos.

typedef struct
{
    long double a0;
    long double a1;
    long double a2;
    int paso;
}ajusteTemperaturaPlano;

// Estructura que guarda coeficientes en forma de arreglos.
typedef struct
{
    long double *a0;
    long double *a1;
    long double *a2;
    int *paso;
}ajusteTemperaturaCristal;

distTempPlano * planoNuevo(int N, int paso, long double ancho, long double alto);
distTempPlano * borraPlano(distTempPlano * plano);
distTempPlano * distTempPlanoCristal(int N, long double ancho, long double alto, long double deltaZ, long double P_pump, long double chi, long double L, long double kth, long double Cp, long double rho, long double w_pump_t,long double w_pump_s, long double alpha,\
                                      int paso, int iteraciones, long double umbral, long double lambda_relax);
ajusteTemperaturaPlano * ajusteCuadraticoPlano(distTempPlano * plano, int iteraciones, long double umbral, long double lambda_relax);
ajusteTemperaturaPlano * ajusteCuadraticoPlanoPonderado(distTempPlano * plano, int iteraciones, long double umbral, long double lambda_relax);
ajusteTemperaturaCristal * vectorAjusteCristal(long double _Complex qInTan, long double _Complex qInSag, long double alpha, long double lambdaIn, long double nLin, \
    long double Power, long double chi, long double L, long double kth, long double Cp, long double rho, long double dn_dT, \
    long double ancho, long double alto, int iteraciones, int pasos, long double umbral, long double lambda_relax, int N, _Bool termico);
void escribePlanoArchivoCompleto(distTempPlano * plano, char *nombre);
void escribePlanoArchivo(distTempPlano * plano, char *nombre);
void escribePlanoArchivoOrigin(distTempPlano * plano, char *nombre);
void escribeCoeficientesAjuste(ajusteTemperaturaCristal *ajuste, int pasos, char *nombre);
ajusteTemperaturaCristal * borraVectorAjuste(ajusteTemperaturaCristal * vector);
ajusteTemperaturaCristal ** matrizAjusteCristal(spotsPump * spots, matriz * epsilon1, long double alpha, long double lambdaIn, long double nLin, \
    long double Power, long double chi, long double L, long double kth, long double Cp, long double rho, long double dn_dT, \
    long double ancho, long double alto, int iteraciones, int pasos, long double umbral, long double lambda_relax, int N, _Bool termico,
                                                CtoCppQProgressBar *bar);
/*ajusteTemperaturaCristal ** matrizAjusteCristal(spotsPump * spots, matriz * epsilon1, long double alpha, long double lambdaIn, long double nLin, \
    long double Power, long double chi, long double L, long double kth, long double Cp, long double rho, long double dn_dT, \
    long double ancho, long double alto, int iteraciones, int pasos, long double umbral, long double lambda_relax, int N, _Bool termico);*/
ajusteTemperaturaCristal ** borraMatrizAjusteCristal(ajusteTemperaturaCristal ** matriz, int numVectores);


#endif // TERMICO_H_INCLUDED

#ifdef __cplusplus
}
#endif
