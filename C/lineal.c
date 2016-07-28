#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrices.h"
//#include <complex.h>
#include <tgmath.h>
#include <stdbool.h>

// Rutinas para cálculo lineal

void anguloLineal(char *conjugado_corto, char *conjugado_largo, long double L, long double n, long double L1, long double L2, long double f1, long double f2,long double * Angulos)
{
    long double A=L/2*(1/powl(n,3)-1/n);
    long double x_1=0.0, x_2=0.0;
    char inf[]="inf";
    char fin[]="fin";
    if(strcmp(conjugado_corto,inf)==0)
        x_1=(A+sqrtl(powl(A,2)+4*powl(f1,2)))/(2*f1);
        else if(strcmp(conjugado_corto,fin)==0)
        x_1=(A*(powl(L1,2)+powl(f1,2))+sqrtl(4*powl(f1,2)*powl(L1,4)+powl(A,2)*powl((powl(f1,2)-powl(L1,2)),2)))/(2*(f1*L1*(L1+A)));
        //else return -2;
    if(strcmp(conjugado_largo,inf)==0)
        x_2=(A+sqrtl(powl(A,2)+4*powl(f2,2)))/(2*f2);
        else if(strcmp(conjugado_largo,fin)==0)
        x_2=(A*(powl(L2,2)+powl(f2,2))+sqrtl(4*powl(f2,2)*powl(L2,4)+powl(A,2)*powl((powl(f2,2)-powl(L2,2)),2)))/(2*(f2*L2*(L2+A)));
        //else return -2;
    Angulos[0]=acosl(x_1);
    Angulos[1]=acosl(x_2);
}

long double distanciaCristal(char *tipo_conjugado, long double f, long double Longitud, long double L, long double n, long double theta)
{
    char inf[]="inf";
    char fin[]="fin";
    long double delta=0.0;
    if(strcmp(tipo_conjugado,inf)==0)
        delta=f/cosl(theta)-L/(2*n);
    else if(strcmp(tipo_conjugado,fin)==0)
        delta=Longitud*f/cosl(theta)/(Longitud-f/cosl(theta))-L/(2*n);
    else
        printf("\nError\n");
    return delta;

}

long double spot_lineal(matriz * ABCD, long double lambda)
{
    assert(ABCD);
    long double estabilidad, A, B, D, w;
    A=obtieneElemento(ABCD,1,1);
    B=obtieneElemento(ABCD,1,2);
    D=obtieneElemento(ABCD,2,2);
    estabilidad=0.5*(A+D);
    if(estabilidad<=-1||estabilidad>=1) w=0;
    else
        w=sqrtl(lambda*fabs(B)/(M_PI*sqrtl(1-powl(((A+D)/2),2))));
    if(isfinite(w)==0)
            w=0;
    return w;
}
void calculoLineal(char *conjugado_corto, char *conjugado_largo, long double n, long double lambda, long double L1,\
                   long double L2,long double L,long double f1,long double f2, \
                   _Bool hacerEM1, _Bool hacerEM2, char *rutaSpotEM1, char *rutaSpotEM2,\
                   matriz *wtEM1, matriz *wsEM1, matriz *wtEM2, matriz *wsEM2, matriz *epsilon1, matriz *epsilon2,matriz *qtEM1, matriz *qsEM1, matriz *qtEM2, matriz *qsEM2)
{
    long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
    assert(strlen(conjugado_corto)==3);
    assert(strlen(conjugado_largo)==3);
    assert(wtEM1);
    assert(wsEM1);
    assert(wtEM2);
    assert(wsEM2);
    assert(qtEM1);
    assert(qsEM1);
    assert(qtEM2);
    assert(qsEM2);

    anguloLineal(conjugado_corto,conjugado_largo,L,n,L1,L2,f1,f2,angulos);
    delta1=distanciaCristal(conjugado_corto,f1,L1,L,n,angulos[0]);
    delta2=distanciaCristal(conjugado_largo,f2,L2,L,n,angulos[1]);
    f1t=f1*cosl(angulos[0]);
    f1s=f1/cosl(angulos[0]);
    f2t=f2*cosl(angulos[1]);
    f2s=f2/cosl(angulos[1]);

    // Bloque EM1
    if(hacerEM1)
    {
        matriz *EM1tan1,*EM1tan2,*EM1tan3,*EM1tan4,*EM1tan5,*EM1tan6,*EM1tan7,*EM1tan8,*EM1tan9,*EM1tan10,*EM1tan11,*EM1tan12,*EM1tan13,*EM1tan14,*EM1tan15;
        matriz *EM1sag1,*EM1sag2,*EM1sag3,*EM1sag4,*EM1sag5,*EM1sag6,*EM1sag7,*EM1sag8,*EM1sag9,*EM1sag10,*EM1sag11,*EM1sag12,*EM1sag13,*EM1sag14,*EM1sag15;

        /* Llenando matrices "estáticas"*/
        // EM1, Caso tangencial
        EM1tan1=llenaMatriz(1,L1,0,1);
        EM1tan2=llenaMatriz(1,0,-1/f1t,1);
        EM1tan3=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
        EM1tan4=llenaMatriz(1,L/powl(n,3),0,1);
        EM1tan5=nuevaMatriz(2,2);//EM1tan5 pendiente
        EM1tan6=llenaMatriz(1,0,-1/f2t,1);
        EM1tan7=llenaMatriz(1,L2,0,1);
        EM1tan8=llenaMatriz(1,0,0,1);
        EM1tan9=llenaMatriz(1,L2,0,1);
        EM1tan10=llenaMatriz(1,0,-1/f2t,1);
        EM1tan11=nuevaMatriz(2,2);    //EM1tan11 pendiente
        EM1tan12=llenaMatriz(1,L/powl(n,3),0,1);
        EM1tan13=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
        EM1tan14=llenaMatriz(1,0,-1/f1t,1);
        EM1tan15=llenaMatriz(1,L1,0,1);

        // EM1, Caso sagital
        EM1sag1=llenaMatriz(1,L1,0,1);
        EM1sag2=llenaMatriz(1,0,-1/f1s,1);
        EM1sag3=nuevaMatriz(2,2);
        EM1sag4=llenaMatriz(1,L/n,0,1);
        EM1sag5=nuevaMatriz(2,2);   //EM1sag5 pendiente
        EM1sag6=llenaMatriz(1,0,-1/f2s,1);
        EM1sag7=llenaMatriz(1,L2,0,1);
        EM1sag8=llenaMatriz(1,0,0,1);
        EM1sag9=llenaMatriz(1,L2,0,1);
        EM1sag10=llenaMatriz(1,0,-1/f2s,1);
        EM1sag11=nuevaMatriz(2,2);    //EM1sag11 pendiente
        EM1sag12=llenaMatriz(1,L/n,0,1);
        EM1sag13=nuevaMatriz(2,2);
        EM1sag14=llenaMatriz(1,0,-1/f1s,1);
        EM1sag15=llenaMatriz(1,L1,0,1);

        // Ciclo
        for(int index1=epsilon1->columnas;index1>0;index1--)
        {
            for(int index2=epsilon2->columnas;index2>0;index2--)
            {

                // Llena las matrices y realiza las multiplicaciones

                // Para epsilon1 EM1: 3 y 13; EM2: 5 y 11
                llenaMatrizSinCrear(EM1tan3,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM1tan13,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM1sag3,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM1sag13,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                // Para epsilon2  EM1: 5 y 11; EM2: 3 y 13
                llenaMatrizSinCrear(EM1tan5,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                llenaMatrizSinCrear(EM1tan11,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                llenaMatrizSinCrear(EM1sag5,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                llenaMatrizSinCrear(EM1sag11,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                matriz *EM1tan_piezas[]={EM1tan15,EM1tan14,EM1tan13,EM1tan12,EM1tan11,EM1tan10,EM1tan9,EM1tan8,EM1tan7,EM1tan6,EM1tan5,EM1tan4,EM1tan3,EM1tan2,EM1tan1};
                matriz *EM1sag_piezas[]={EM1sag15,EM1sag14,EM1sag13,EM1sag12,EM1sag11,EM1sag10,EM1sag9,EM1sag8,EM1sag7,EM1sag6,EM1sag5,EM1sag4,EM1sag3,EM1sag2,EM1sag1};
                matriz *EM1tan, *EM1sag;
                EM1tan=variasMultMatriciales(EM1tan_piezas,15);
                EM1sag=variasMultMatriciales(EM1sag_piezas,15);
                long double spot_EM1tan=spot_lineal(EM1tan,lambda);
                long double spot_EM1sag=spot_lineal(EM1sag,lambda);
                long double _Complex q_EM1tan=1/((-I*lambda)/(M_PI*powl(spot_EM1tan,2)));
                long double _Complex q_EM1sag=1/((-I*lambda)/(M_PI*powl(spot_EM1sag,2)));

                if(spot_EM1tan>0.1){
                    fijaElemento(wtEM1,index1,index2,0);
                    fijaElemento(qtEM1,index1,index2,0);}
                else{
                    fijaElemento(wtEM1,index1,index2,spot_EM1tan);
                    fijaElemento(qtEM1,index1,index2,q_EM1tan);
                }
                if(spot_EM1sag>0.1){
                    fijaElemento(wsEM1,index1,index2,0);
                    fijaElemento(qsEM1,index1,index2,0);}
                else{
                    fijaElemento(wsEM1,index1,index2,spot_EM1sag);
                    fijaElemento(qsEM1,index1,index2,q_EM1sag);
                }
                // Limpia matrices
                EM1tan=borraMatriz(EM1tan);
                EM1sag=borraMatriz(EM1sag);
            }
        }

        // Escritura

        imprimeListaReal(epsilon1,epsilon2,wtEM1,wsEM1,rutaSpotEM1);
        //imprimeListaReal(epsilon1,epsilon2,qtEM1,qsEM1,rutaQEM1);

        // Limpieza

        EM1tan1=borraMatriz(EM1tan1);
        EM1tan2=borraMatriz(EM1tan2);
        EM1tan3=borraMatriz(EM1tan3);
        EM1tan4=borraMatriz(EM1tan4);
        EM1tan5=borraMatriz(EM1tan5);
        EM1tan6=borraMatriz(EM1tan6);
        EM1tan7=borraMatriz(EM1tan7);
        EM1tan8=borraMatriz(EM1tan8);
        EM1tan9=borraMatriz(EM1tan9);
        EM1tan10=borraMatriz(EM1tan10);
        EM1tan11=borraMatriz(EM1tan11);
        EM1tan12=borraMatriz(EM1tan12);
        EM1tan13=borraMatriz(EM1tan13);
        EM1tan14=borraMatriz(EM1tan14);
        EM1tan15=borraMatriz(EM1tan15);

        EM1sag1=borraMatriz(EM1sag1);
        EM1sag2=borraMatriz(EM1sag2);
        EM1sag3=borraMatriz(EM1sag3);
        EM1sag4=borraMatriz(EM1sag4);
        EM1sag5=borraMatriz(EM1sag5);
        EM1sag6=borraMatriz(EM1sag6);
        EM1sag7=borraMatriz(EM1sag7);
        EM1sag8=borraMatriz(EM1sag8);
        EM1sag9=borraMatriz(EM1sag9);
        EM1sag10=borraMatriz(EM1sag10);
        EM1sag11=borraMatriz(EM1sag11);
        EM1sag12=borraMatriz(EM1sag12);
        EM1sag13=borraMatriz(EM1sag13);
        EM1sag14=borraMatriz(EM1sag14);
        EM1sag15=borraMatriz(EM1sag15);

    }

    if(hacerEM2)
    {
        matriz *EM2tan1,*EM2tan2,*EM2tan3,*EM2tan4,*EM2tan5,*EM2tan6,*EM2tan7,*EM2tan8,*EM2tan9,*EM2tan10,*EM2tan11,*EM2tan12,*EM2tan13,*EM2tan14,*EM2tan15;
        matriz *EM2sag1,*EM2sag2,*EM2sag3,*EM2sag4,*EM2sag5,*EM2sag6,*EM2sag7,*EM2sag8,*EM2sag9,*EM2sag10,*EM2sag11,*EM2sag12,*EM2sag13,*EM2sag14,*EM2sag15;

        // EM2, Caso tangencial
        EM2tan1=llenaMatriz(1,L2,0,1);
        EM2tan2=llenaMatriz(1,0,-1/f2t,1);
        EM2tan3=nuevaMatriz(2,2);//EM2tan3=llenaMatriz(1,delta2+eps,0,1);
        EM2tan4=llenaMatriz(1,L/powl(n,3),0,1);
        EM2tan5=nuevaMatriz(2,2);
        EM2tan6=llenaMatriz(1,0,-1/f1t,1);
        EM2tan7=llenaMatriz(1,L1,0,1);
        EM2tan8=llenaMatriz(1,0,0,1);
        EM2tan9=llenaMatriz(1,L1,0,1);
        EM2tan10=llenaMatriz(1,0,-1/f1t,1);
        EM2tan11=nuevaMatriz(2,2);
        EM2tan12=llenaMatriz(1,L/powl(n,3),0,1);
        EM2tan13=nuevaMatriz(2,2);//EM2tan13=llenaMatriz(1,delta1,0,1);
        EM2tan14=llenaMatriz(1,0,-1/f2t,1);
        EM2tan15=llenaMatriz(1,L2,0,1);

        // EM2, Caso sagital
        EM2sag1=llenaMatriz(1,L2,0,1);
        EM2sag2=llenaMatriz(1,0,-1/f2s,1);
        EM2sag3=nuevaMatriz(2,2);   //EM2sag3=llenaMatriz(1,delta2+e,0,1);
        EM2sag4=llenaMatriz(1,L/n,0,1);
        EM2sag5=nuevaMatriz(2,2);
        EM2sag6=llenaMatriz(1,0,-1/f1s,1);
        EM2sag7=llenaMatriz(1,L1,0,1);
        EM2sag8=llenaMatriz(1,0,0,1);
        EM2sag9=llenaMatriz(1,L1,0,1);
        EM2sag10=llenaMatriz(1,0,-1/f1s,1);
        EM2sag11=nuevaMatriz(2,2);
        EM2sag12=llenaMatriz(1,L/n,0,1);
        EM2sag13=nuevaMatriz(2,2);  //EM2sag13=llenaMatriz(1,delta2+e,0,1);
        EM2sag14=llenaMatriz(1,0,-1/f2s,1);
        EM2sag15=llenaMatriz(1,L2,0,1);

        for(int index1=epsilon1->columnas;index1>0;index1--)
        {
            for(int index2=epsilon2->columnas;index2>0;index2--)
            {
                llenaMatrizSinCrear(EM2tan5,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM2tan11,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM2sag5,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM2sag11,1,delta1+obtieneElemento(epsilon1,1,index1),0,1);
                llenaMatrizSinCrear(EM2tan3,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                llenaMatrizSinCrear(EM2tan13,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                llenaMatrizSinCrear(EM2sag3,1,delta2+obtieneElemento(epsilon2,1,index2),0,1);
                llenaMatrizSinCrear(EM2sag13,1,delta2+obtieneElemento(epsilon2,1,index2),0,1); // Optimizar asignaciones.
                matriz *EM2tan_piezas[]={EM2tan15,EM2tan14,EM2tan13,EM2tan12,EM2tan11,EM2tan10,EM2tan9,EM2tan8,EM2tan7,EM2tan6,EM2tan5,EM2tan4,EM2tan3,EM2tan2,EM2tan1};
                matriz *EM2sag_piezas[]={EM2sag15,EM2sag14,EM2sag13,EM2sag12,EM2sag11,EM2sag10,EM2sag9,EM2sag8,EM2sag7,EM2sag6,EM2sag5,EM2sag4,EM2sag3,EM2sag2,EM2sag1};
                matriz *EM2tan, *EM2sag;
                EM2tan=variasMultMatriciales(EM2tan_piezas,15);
                EM2sag=variasMultMatriciales(EM2sag_piezas,15);
                long double spot_EM2tan=spot_lineal(EM2tan,lambda);
                long double spot_EM2sag=spot_lineal(EM2sag,lambda);
                long double _Complex q_EM2tan=1/((-I*lambda)/(M_PI*powl(spot_EM2tan,2)));
                long double _Complex q_EM2sag=1/((-I*lambda)/(M_PI*powl(spot_EM2sag,2)));
                if(spot_EM2tan>0.1){
                    fijaElemento(wtEM2,index1,index2,0);
                    fijaElemento(qtEM2,index1,index2,0);}
                else{
                    fijaElemento(wtEM2,index1,index2,spot_EM2tan);
                    fijaElemento(qtEM2,index1,index2,q_EM2tan);
                }
                if(spot_EM2sag>0.1){
                    fijaElemento(wsEM2,index1,index2,0);
                    fijaElemento(qsEM2,index1,index2,0);}
                else{
                    fijaElemento(wsEM2,index1,index2,spot_EM2sag);
                    fijaElemento(qsEM2,index1,index2,q_EM2sag);
                }
                EM2tan=borraMatriz(EM2tan);
                EM2sag=borraMatriz(EM2sag);

            }
        }

        // Escritura
        imprimeListaReal(epsilon1,epsilon2,wtEM2,wsEM2,rutaSpotEM2);
        //imprimeListaReal(epsilon1,epsilon2,qtEM2,qsEM2,rutaQEM2);

        // Limpieza

        EM2tan1=borraMatriz(EM2tan1);
        EM2tan2=borraMatriz(EM2tan2);
        EM2tan3=borraMatriz(EM2tan3);
        EM2tan4=borraMatriz(EM2tan4);
        EM2tan5=borraMatriz(EM2tan5);
        EM2tan6=borraMatriz(EM2tan6);
        EM2tan7=borraMatriz(EM2tan7);
        EM2tan8=borraMatriz(EM2tan8);
        EM2tan9=borraMatriz(EM2tan9);
        EM2tan10=borraMatriz(EM2tan10);
        EM2tan11=borraMatriz(EM2tan11);
        EM2tan12=borraMatriz(EM2tan12);
        EM2tan13=borraMatriz(EM2tan13);
        EM2tan14=borraMatriz(EM2tan14);
        EM2tan15=borraMatriz(EM2tan15);

        EM2sag1=borraMatriz(EM2sag1);
        EM2sag2=borraMatriz(EM2sag2);
        EM2sag3=borraMatriz(EM2sag3);
        EM2sag4=borraMatriz(EM2sag4);
        EM2sag5=borraMatriz(EM2sag5);
        EM2sag6=borraMatriz(EM2sag6);
        EM2sag7=borraMatriz(EM2sag7);
        EM2sag8=borraMatriz(EM2sag8);
        EM2sag9=borraMatriz(EM2sag9);
        EM2sag10=borraMatriz(EM2sag10);
        EM2sag11=borraMatriz(EM2sag11);
        EM2sag12=borraMatriz(EM2sag12);
        EM2sag13=borraMatriz(EM2sag13);
        EM2sag14=borraMatriz(EM2sag14);
        EM2sag15=borraMatriz(EM2sag15);

    }

}

/* Las rutas de propagación del haz se denominan como (1) y (2). La ruta (1) contempla la propagación del espejo plano 1 (FM1) al espejo curvo 1 (CM1),
de CM1 al espejo curvo 2 (CM2) a través del medio Kerr, de CM2 al espejo plano 2 (FM2) y de FM2 de vuelta a FM1.
La ruta (2) describe la propagación de FM1 a FM2, de FM2 a CM2, de CM2 a CM1 atravezando el medio Kerr y de CM1 a FM1.

Se realizarán cálculos para los planos tangencial y sagital de ambas rutas, compensando el astigmatismo lineal.
El conjugado corto será equivalente al tramo entre FM1 a CM1 (tramo a), y el conjugado largo será la suma de tramos de CM2 a FM2 y de FM2 a FM1.
Para satisfacer la condición de autoconsistencia, sólo se admitirán conjugados del mismo tipo, esto es, infinito-infinito o finito-finito. */
