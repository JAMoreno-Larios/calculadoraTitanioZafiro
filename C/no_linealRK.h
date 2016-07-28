#ifdef __cplusplus
extern "C" {
#endif

#ifndef NO_LINEALRK_H_INCLUDED
#define NO_LINEALRK_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <_Complex.h>
#include <string.h>
#include <assert.h>

// Resolver el sistema de ecuaciones diferenciales acopladas mostradas en el artículo de D. Huang.
// Emplea Runge-Kutta de 4to orden. Se define a p=1/q, por lo que pReal=Re(p) y pImag=Im(p).
// xi queda definido por sqrt(1-P_laser/P_cr). Registrar la posición al momento de evaluar

// Ecuaciones diferenciales

long double dpReal(long double pReal, long double pImag, long double n0, long double n2, long double lambda0, long double P_laser);
long double dpImag(long double pReal, long double pImag, long double n0);
long double * pasoKerrRK(long double h, long double pReal_ini, long double pImag_ini, long double n0, long double n2, long double lambda0, long double P_laser);
long double _Complex propagacionKerrRK(long double _Complex q_in, long double lambda0, long double L, long double n0, long double n2, long double P_laser, int pasos);
matriz * propNoLinealRK_tanEM1(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_sagEM1(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_tanEM2(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_sagEM2(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_AnilloRuta1Tan(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_AnilloRuta1Sag(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_AnilloRuta2Tan(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);
matriz * propNoLinealRK_AnilloRuta2Sag(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral);



#endif // NO_LINEALRK_H_INCLUDED

#ifdef __cplusplus
}
#endif
