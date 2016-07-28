#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "matrices.h"
#include "lineal.h"
#include "no_lineal.h"
#include "ctocppqprogressbar.h"


typedef struct
{
    long double _Complex *qOutTan;
    long double _Complex *qOutSag;
    int numSpots;
}spotsPump;

spotsPump * nuevoSpotsPump(int numSpots)
{
    if(numSpots<=0) return NULL;
    spotsPump *sP = (spotsPump*)calloc(1,sizeof(spotsPump));
    sP->qOutTan=(long double _Complex *)calloc(numSpots,sizeof(long double _Complex));
    sP->qOutSag=(long double _Complex *)calloc(numSpots,sizeof(long double _Complex));
    sP->numSpots=numSpots;
}

spotsPump * borraSpotsPump( spotsPump * sP)
{
    assert(sP);
    // Libera datos de la matriz
    assert(sP->qOutTan);
    free(sP->qOutTan);
    sP->qOutTan=NULL;
    assert(sP->qOutSag);
    free(sP->qOutSag);
    sP->qOutSag=NULL;
    free(sP);
    sP=NULL;
    return sP;
}

spotsPump * acopleOptico(long double wFuenteLaserTan, long double wFuenteLaserSag, long double divergencia, long double lambda0,\
                         long double La, long double tL, long double R11, long double R12,\
                         long double nL, long double thetaL, long double Lb, long double tE,\
                         long double R21, long double R22, long double nE, long double thetaE,\
                          long double nC, long double delta1, matriz * epsilon1)
                        {
                            long double spotWaistTan=2*lambda0/M_PI/divergencia;
                            long double spotWaistSag=2*lambda0/M_PI/divergencia;
                            long double LtoInTan=sqrtl(powl(wFuenteLaserTan/spotWaistTan,2)-1)*(M_PI*powl(spotWaistTan,2))/lambda0;
                            long double LtoInSag=sqrtl(powl(wFuenteLaserSag/spotWaistSag,2)-1)*(M_PI*powl(spotWaistSag,2))/lambda0;
                            long double _Complex qInTan=La+LtoInTan+I*M_PI*1*powl(spotWaistTan,2)/lambda0;
                            long double _Complex qInSag=La+LtoInSag+I*M_PI*1*powl(spotWaistSag,2)/lambda0;

                            assert(epsilon1->filas==1);
                            long double Lc;

                            // Prepara estructura de salida
                            spotsPump *Salida=nuevoSpotsPump(epsilon1->columnas);

                            // Variables de trabajo

                            long double tanA11, tanB11, tanC11, tanD11, tanA12, tanB12, tanC12, tanD12,\
                             tanA21, tanB21, tanC21, tanD21, tanA22, tanB22, tanC22, tanD22;
                            // Entrada y salida de lentes gruesas - tangencial
                            // Lente 1
                            tanA11=sqrtl(1-powl(sinl(thetaL)/nL,2));
                            tanB11=0;
                            tanC11=(sqrtl(1-powl(sinl(thetaL)/nL,2))-nL*cosl(thetaL))/(R11*cosl(thetaL)*sqrtl(1-powl(sinl(thetaL)/nL,2)));
                            tanD11=1/tanA11;
                            tanA12=1/tanA11;
                            tanB12=0;
                            tanC12=-tanC11;
                            tanD12=tanA11;

                            //Lente 2 (Espejo cóncavo) -tangencial

                            tanA21=sqrtl(1-powl(sinl(thetaL)/nL,2));
                            tanB21=0;
                            tanC21=(sqrtl(1-powl(sinl(thetaL)/nL,2))-nL*cosl(thetaL))/(R12*cosl(thetaL)*sqrtl(1-powl(sinl(thetaL)/nL,2)));
                            tanD21=1/tanA11;
                            tanA22=1/tanA11;
                            tanB22=0;
                            tanC22=-tanC11;
                            tanD22=tanA11;

                            long double sagA11, sagB11, sagC11, sagD11, sagA12, sagB12, sagC12, sagD12,\
                             sagA21, sagB21, sagC21, sagD21, sagA22, sagB22, sagC22, sagD22;
                            // Entrada y salida de lentes gruesas - sagital
                            // Lente 1
                            sagA11=sqrtl(1-powl(sinl(thetaE)/nE,2));
                            sagB11=0;
                            sagC11=(sqrtl(1-powl(sinl(thetaE)/nE,2))-nE*cosl(thetaE))/(R21*cosl(thetaE)*sqrtl(1-powl(sinl(thetaE)/nE,2)));
                            sagD11=1/sagA11;
                            sagA12=1/sagA11;
                            sagB12=0;
                            sagC12=-sagC11;
                            sagD12=sagA11;

                            //Lente 2 (Espejo cóncavo) - sagital

                            sagA21=sqrtl(1-powl(sinl(thetaE)/nE,2));
                            sagB21=0;
                            sagC21=(sqrtl(1-powl(sinl(thetaE)/nE,2))-nE*cosl(thetaE))/(R22*cosl(thetaE)*sqrtl(1-powl(sinl(thetaE)/nE,2)));
                            sagD21=1/sagA11;
                            sagA22=1/sagA11;
                            sagB22=0;
                            sagC22=-sagC11;
                            sagD22=sagA11;

                            // Declarando matrices
                            matriz *Tan1, *Tan2, *Tan3, *Tan4, *Tan5, *Tan6, *Tan7,\
                             *Tan8, *Tan9, *Tan10;
                            matriz *Sag1, *Sag2, *Sag3, *Sag4, *Sag5, *Sag6, *Sag7,\
                             *Sag8, *Sag9, *Sag10;

                            // Plano tangencial

                            Tan1=llenaMatriz(1, La, 0, 1);
                            Tan2=llenaMatriz(tanA11,tanB11,tanC11,tanD11);
                            Tan3=llenaMatriz(1,tL/nL,0,1);
                            Tan4=llenaMatriz(tanA12,tanB12,tanC12,tanD12);
                            Tan5=llenaMatriz(1,Lb,0,1);
                            Tan6=llenaMatriz(tanA21,tanB21,tanC21,tanD21);
                            Tan7=llenaMatriz(1,tE/nE,0,1);
                            Tan8=llenaMatriz(tanA22,tanB22,tanC22,tanD22);
                            Tan9=nuevaMatriz(2,2); // Variable
                            Tan10=llenaMatriz(nC,0,0,1/nC);

                            // Plano tangencial

                            Sag1=llenaMatriz(1, La, 0, 1);
                            Sag2=llenaMatriz(sagA11,sagB11,sagC11,sagD11);
                            Sag3=llenaMatriz(1,tL/nL,0,1);
                            Sag4=llenaMatriz(sagA12,sagB12,sagC12,sagD12);
                            Sag5=llenaMatriz(1,Lb,0,1);
                            Sag6=llenaMatriz(sagA21,sagB21,sagC21,sagD21);
                            Sag7=llenaMatriz(1,tE/nE,0,1);
                            Sag8=llenaMatriz(sagA22,sagB22,sagC22,sagD22);
                            Sag9=nuevaMatriz(2,2); // Variable
                            Sag10=llenaMatriz(1,0,0,1);

                            int i;
                            for(i=0;i<epsilon1->columnas;i++)
                            {
                                Lc=delta1+creall(obtieneElemento(epsilon1,1,i+1));
                                llenaMatrizSinCrear(Tan9,1,Lc,0,1);
                                llenaMatrizSinCrear(Sag9,1,Lc,0,1);
                                matriz *propTanPiezas[10]={Tan10,Tan9,Tan8,Tan7,Tan6,Tan5,Tan4,Tan3,Tan2,Tan1};
                                matriz *propSagPiezas[10]={Sag10,Sag9,Sag8,Sag7,Sag6,Sag5,Sag4,Sag3,Sag2,Sag1};
                                matriz *propTan, *propSag;
                                propTan=variasMultMatriciales(propTanPiezas,10);
                                propSag=variasMultMatriciales(propSagPiezas,10);

                                Salida->qOutTan[i]=prop_q(propTan,qInTan,1,nC);
                                Salida->qOutSag[i]=prop_q(propSag,qInSag,1,nC);

                                // Limpieza local
                                propTan=borraMatriz(propTan);
                                propSag=borraMatriz(propSag);
                            }

                            // Limpieza final
                            Tan1=borraMatriz(Tan1);
                            Tan2=borraMatriz(Tan2);
                            Tan3=borraMatriz(Tan3);
                            Tan4=borraMatriz(Tan4);
                            Tan5=borraMatriz(Tan5);
                            Tan6=borraMatriz(Tan6);
                            Tan7=borraMatriz(Tan7);
                            Tan8=borraMatriz(Tan8);
                            Tan9=borraMatriz(Tan9);
                            Tan10=borraMatriz(Tan10);
                            Sag1=borraMatriz(Sag1);
                            Sag2=borraMatriz(Sag2);
                            Sag3=borraMatriz(Sag3);
                            Sag4=borraMatriz(Sag4);
                            Sag5=borraMatriz(Sag5);
                            Sag6=borraMatriz(Sag6);
                            Sag7=borraMatriz(Sag7);
                            Sag8=borraMatriz(Sag8);
                            Sag9=borraMatriz(Sag9);
                            Sag10=borraMatriz(Sag10);

                            return Salida;
                        }


spotsPump * acopleOpticoBarra(long double wFuenteLaserTan, long double wFuenteLaserSag, long double divergencia, long double lambda0,\
                         long double La, long double tL, long double R11, long double R12,\
                         long double nL, long double thetaL, long double Lb, long double tE,\
                         long double R21, long double R22, long double nE, long double thetaE,\
                          long double nC, long double delta1, matriz * epsilon1, CtoCppQProgressBar * barra)
                        {
                            long double spotWaistTan=2*lambda0/M_PI/divergencia;
                            long double spotWaistSag=2*lambda0/M_PI/divergencia;
                            long double LtoInTan=sqrtl(powl(wFuenteLaserTan/spotWaistTan,2)-1)*(M_PI*powl(spotWaistTan,2))/lambda0;
                            long double LtoInSag=sqrtl(powl(wFuenteLaserSag/spotWaistSag,2)-1)*(M_PI*powl(spotWaistSag,2))/lambda0;
                            long double _Complex qInTan=La+LtoInTan+I*M_PI*1*powl(spotWaistTan,2)/lambda0;
                            long double _Complex qInSag=La+LtoInSag+I*M_PI*1*powl(spotWaistSag,2)/lambda0;

                            assert(epsilon1->filas==1);
                            long double Lc;

                            // Prepara estructura de salida
                            spotsPump *Salida=nuevoSpotsPump(epsilon1->columnas);

                            // Variables de trabajo

                            long double tanA11, tanB11, tanC11, tanD11, tanA12, tanB12, tanC12, tanD12,\
                             tanA21, tanB21, tanC21, tanD21, tanA22, tanB22, tanC22, tanD22;
                            // Entrada y salida de lentes gruesas - tangencial
                            // Lente 1
                            tanA11=sqrtl(1-powl(sinl(thetaL)/nL,2));
                            tanB11=0;
                            tanC11=(sqrtl(1-powl(sinl(thetaL)/nL,2))-nL*cosl(thetaL))/(R11*cosl(thetaL)*sqrtl(1-powl(sinl(thetaL)/nL,2)));
                            tanD11=1/tanA11;
                            tanA12=1/tanA11;
                            tanB12=0;
                            tanC12=-tanC11;
                            tanD12=tanA11;

                            //Lente 2 (Espejo cóncavo) -tangencial

                            tanA21=sqrtl(1-powl(sinl(thetaL)/nL,2));
                            tanB21=0;
                            tanC21=(sqrtl(1-powl(sinl(thetaL)/nL,2))-nL*cosl(thetaL))/(R12*cosl(thetaL)*sqrtl(1-powl(sinl(thetaL)/nL,2)));
                            tanD21=1/tanA11;
                            tanA22=1/tanA11;
                            tanB22=0;
                            tanC22=-tanC11;
                            tanD22=tanA11;

                            long double sagA11, sagB11, sagC11, sagD11, sagA12, sagB12, sagC12, sagD12,\
                             sagA21, sagB21, sagC21, sagD21, sagA22, sagB22, sagC22, sagD22;
                            // Entrada y salida de lentes gruesas - sagital
                            // Lente 1
                            sagA11=sqrtl(1-powl(sinl(thetaE)/nE,2));
                            sagB11=0;
                            sagC11=(sqrtl(1-powl(sinl(thetaE)/nE,2))-nE*cosl(thetaE))/(R21*cosl(thetaE)*sqrtl(1-powl(sinl(thetaE)/nE,2)));
                            sagD11=1/sagA11;
                            sagA12=1/sagA11;
                            sagB12=0;
                            sagC12=-sagC11;
                            sagD12=sagA11;

                            //Lente 2 (Espejo cóncavo) - sagital

                            sagA21=sqrtl(1-powl(sinl(thetaE)/nE,2));
                            sagB21=0;
                            sagC21=(sqrtl(1-powl(sinl(thetaE)/nE,2))-nE*cosl(thetaE))/(R22*cosl(thetaE)*sqrtl(1-powl(sinl(thetaE)/nE,2)));
                            sagD21=1/sagA11;
                            sagA22=1/sagA11;
                            sagB22=0;
                            sagC22=-sagC11;
                            sagD22=sagA11;

                            // Declarando matrices
                            matriz *Tan1, *Tan2, *Tan3, *Tan4, *Tan5, *Tan6, *Tan7,\
                             *Tan8, *Tan9, *Tan10;
                            matriz *Sag1, *Sag2, *Sag3, *Sag4, *Sag5, *Sag6, *Sag7,\
                             *Sag8, *Sag9, *Sag10;

                            // Plano tangencial

                            Tan1=llenaMatriz(1, La, 0, 1);
                            Tan2=llenaMatriz(tanA11,tanB11,tanC11,tanD11);
                            Tan3=llenaMatriz(1,tL/nL,0,1);
                            Tan4=llenaMatriz(tanA12,tanB12,tanC12,tanD12);
                            Tan5=llenaMatriz(1,Lb,0,1);
                            Tan6=llenaMatriz(tanA21,tanB21,tanC21,tanD21);
                            Tan7=llenaMatriz(1,tE/nE,0,1);
                            Tan8=llenaMatriz(tanA22,tanB22,tanC22,tanD22);
                            Tan9=nuevaMatriz(2,2); // Variable
                            Tan10=llenaMatriz(nC,0,0,1/nC);

                            // Plano tangencial

                            Sag1=llenaMatriz(1, La, 0, 1);
                            Sag2=llenaMatriz(sagA11,sagB11,sagC11,sagD11);
                            Sag3=llenaMatriz(1,tL/nL,0,1);
                            Sag4=llenaMatriz(sagA12,sagB12,sagC12,sagD12);
                            Sag5=llenaMatriz(1,Lb,0,1);
                            Sag6=llenaMatriz(sagA21,sagB21,sagC21,sagD21);
                            Sag7=llenaMatriz(1,tE/nE,0,1);
                            Sag8=llenaMatriz(sagA22,sagB22,sagC22,sagD22);
                            Sag9=nuevaMatriz(2,2); // Variable
                            Sag10=llenaMatriz(1,0,0,1);

                            int i;
                            CtoCppQProgressBar_setMinimum(barra,0);
                            CtoCppQProgressBar_setMaximum(barra,epsilon1->columnas);
                            for(i=0;i<epsilon1->columnas;i++)
                            {
                                Lc=delta1+creall(obtieneElemento(epsilon1,1,i+1));
                                llenaMatrizSinCrear(Tan9,1,Lc,0,1);
                                llenaMatrizSinCrear(Sag9,1,Lc,0,1);
                                matriz *propTanPiezas[10]={Tan10,Tan9,Tan8,Tan7,Tan6,Tan5,Tan4,Tan3,Tan2,Tan1};
                                matriz *propSagPiezas[10]={Sag10,Sag9,Sag8,Sag7,Sag6,Sag5,Sag4,Sag3,Sag2,Sag1};
                                matriz *propTan, *propSag;
                                propTan=variasMultMatriciales(propTanPiezas,10);
                                propSag=variasMultMatriciales(propSagPiezas,10);

                                Salida->qOutTan[i]=prop_q(propTan,qInTan,1,nC);
                                Salida->qOutSag[i]=prop_q(propSag,qInSag,1,nC);

                                // Limpieza local
                                propTan=borraMatriz(propTan);
                                propSag=borraMatriz(propSag);

                                // Actualiza barra
                                CtoCppQProgressBar_setValue(barra,i+1);
                                CtoCppQProgressBar_repaint(barra);
                            }

                            // Limpieza final
                            Tan1=borraMatriz(Tan1);
                            Tan2=borraMatriz(Tan2);
                            Tan3=borraMatriz(Tan3);
                            Tan4=borraMatriz(Tan4);
                            Tan5=borraMatriz(Tan5);
                            Tan6=borraMatriz(Tan6);
                            Tan7=borraMatriz(Tan7);
                            Tan8=borraMatriz(Tan8);
                            Tan9=borraMatriz(Tan9);
                            Tan10=borraMatriz(Tan10);
                            Sag1=borraMatriz(Sag1);
                            Sag2=borraMatriz(Sag2);
                            Sag3=borraMatriz(Sag3);
                            Sag4=borraMatriz(Sag4);
                            Sag5=borraMatriz(Sag5);
                            Sag6=borraMatriz(Sag6);
                            Sag7=borraMatriz(Sag7);
                            Sag8=borraMatriz(Sag8);
                            Sag9=borraMatriz(Sag9);
                            Sag10=borraMatriz(Sag10);

                            return Salida;
                        }
