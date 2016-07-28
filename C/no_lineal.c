#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <complex.h>
#include <tgmath.h>
#include "matrices.h"
#include "lineal.h"
#include "error_iteraciones.h"

// Constantes para propagación en medio Kerr

    int pasosKerr=1000;


// Funciones para propagacion en rutina no lineal

// Nota: Las rutinas correran indefinidamente hasta que converjan.


long double _Complex prop_q(matriz * ABCD, long double _Complex q_in, long double n1, long double n2) // Propagación de un haz gaussiano a través de un sistema óptico ABCD - Siegman
{
    assert(ABCD);
    assert(q_in);
    assert(n1!=0&&n2!=0);
    long double _Complex q_out;
    q_out=n2*(obtieneElemento(ABCD,1,1)*q_in+n1*obtieneElemento(ABCD,1,2))/(obtieneElemento(ABCD,2,1)*q_in+n1*obtieneElemento(ABCD,2,2));
    return q_out;
}


long double spot_q(long double _Complex q, long double n, long double lambda)
{
    long double w=sqrtl(-lambda/(n*M_PI*cimagl(1/q)));
    return w;
}

long double radio_q(long double _Complex q)
{
    long double R=creal(1/creall(1/q));
    return R;
}

matriz * spot_q_matriz(matriz * q, long double n, long double lambda)
{
    matriz *w=nuevaMatriz(q->filas,q->columnas);
    long double _Complex q_ij;
    long double w_ij;
    for(int i=q->filas;i>0;i--)
        for(int j=q->columnas;j>0;j--)
            {
                q_ij=obtieneElemento(q,i,j);
                w_ij=creal(sqrtl(-lambda/(n*M_PI*cimag(1/q_ij))));
                if(isfinite(w_ij)==0)
                    w_ij=INFINITY;
                else
                fijaElemento(w,i,j,w_ij);
            }
    return w;
}


// Propagación de haz gaussiano en medio Kerr - Ecuaciones desacopladas, método matricial.
// Original
long double _Complex propagacionKerr(int pasos, long double L, long double n0,long double n2, long double w_pump,long double _Complex q_in, long double chi, long double kth, long double Cp, long double rho, long double dn_dv,long double P_pump,long double P_laser, long double lambda0)
{
        assert(q_in);
        // Variables locales
        long double dz=L/pasos; // Grosor de lámina
        long double lambda1=lambda0/n0, w_in, parabola;
        long double _Complex *qProp = (long double _Complex*)calloc(pasos,sizeof(long double _Complex));
        long double _Complex A, B, C, D;
        assert(qProp);
        matriz * M_i;
        M_i=nuevaMatriz(2,2);
        int i=pasos-1;
        qProp[i]=q_in;

        for (i=pasos-1;i>0;i--)
        {
            w_in=sqrtl(-lambda1/(M_PI*cimagl(1/qProp[i])));
            // Factor de parábola: parábola=1/h^2
            //parabola=1/n0*dn_dv*(chi*P_pump/L)/(M_PI*powl(w_pump,2)*kth*Cp*rho)+1/n0*(8*n2*P_laser)/(M_PI*powl(w_in,4)); // Kerr y térmico
            //parabola=1/n0*(8*n2*P_laser)/(M_PI*powl(w_in,4)); // Sólo efecto Kerr - truncado por Taylor
            parabola=1/n0*(n2*P_laser)/(2.0*M_PI*powl(w_in,4)); // Sólo efecto Kerr - Ajuste cuadrático
            // Calculando parámetros de la matriz
            A=1.0;
            B=dz/n0;
            C=-n0*parabola*dz;
            D=1.0;
            llenaMatrizSinCrear(M_i,A,B,C,D);
            qProp[i-1]=prop_q(M_i,qProp[i],n0,n0);
        }

        borraMatriz(M_i);
        long double _Complex regresar=qProp[0];
        assert(qProp!=NULL);
        free(qProp);
        qProp=NULL;
        return regresar;
}


long double _Complex * propagacionKerrMatriz(int pasos, long double L, long double n0,long double n2, long double w_pump,long double _Complex q_in, long double chi, long double kth, long double Cp, long double rho, long double dn_dv,long double P_pump,long double P_laser, long double lambda0)
{
        assert(q_in);
        // Variables locales
        long double dz=L/pasos; // Grosor de lámina
        long double lambda1=lambda0/n0, w_in, parabola;
        long double _Complex *qProp = (long double _Complex*)calloc(pasos,sizeof(long double _Complex));
        long double _Complex A, B, C, D;
        assert(qProp);
        matriz * M_i;
        M_i=nuevaMatriz(2,2);
        qProp[0]=q_in;
        int i=0;
        for (i=0;i<pasos-1;i++)
        {
            w_in=sqrtl(-lambda1/(M_PI*cimag(1/ qProp[i])));
            // Factor de parábola: parábola=1/h^2
            //parabola=1/n0*dn_dv*(chi*P_pump/L)/(M_PI*powl(w_pump,2)*kth*Cp*rho)+1/n0*(8*n2*P_laser)/(M_PI*powl(w_in,4)); // Kerr y térmico
            parabola=1/n0*(8.0*n2*P_laser)/(M_PI*powl(w_in,4.0)); // Sólo efecto Kerr
            // Calculando parámetros de la matriz
            A=1.0;
            B=dz/n0;
            C=-n0*parabola*dz;
            D=1.0;
            llenaMatrizSinCrear(M_i,A,B,C,D);
            // Propagando el haz
            qProp[i+1]=prop_q(M_i,qProp[i],n0,n0);
        }
        borraMatriz(M_i);
        return qProp;
}

// Busca el valor menos astigmático en la matriz.
void confNoAstigmatica(matriz * A, long double _Complex *casiUno, int *fila, int *columna)
{
    *casiUno = obtieneElemento(A,1,1);
    *fila=1;
    *columna=1;
    for(int i=1;i<A->filas;++i)
    {
        for(int j=1;j<A->columnas;++j)
            {
                if(cabsl(obtieneElemento(A,i,j)-1)<cabsl(*casiUno-1))
                {
                    *casiUno=obtieneElemento(A,i,j);
                    *fila=i;
                    *columna=j;
                }
            }
    }
    assert(obtieneElemento(A,*fila,*columna)==*casiUno);
}
