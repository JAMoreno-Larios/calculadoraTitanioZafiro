#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
//#include <tgmath.h>
#include "matrices.h"
#include "lineal.h"
#include "termico.h"
#include "no_lineal.h"
#include "error_iteraciones.h"
#include "ctocppqprogressbar.h"

// Programa para la propagación no lineal considerando haces astigmáticos para efecto Kerr y térmico

// Constantes para propagación en medio Kerr

    //int pasosKerr=1000;

// Estructura para sacar dos números complejos de la propagación Kerr acoplada

typedef struct
{
    long double _Complex qTan;
    long double _Complex qSag;
}propAstigmatica;

// Original
propAstigmatica * propagacionKerrAcoplada(int pasos, long double L, long double n0,long double n2, long double _Complex qInTan, \
                                            long double _Complex qInSag, long double chi, long double kth, long double Cp, \
                                            long double rho, long double dn_dv,long double P_laser, long double lambda0, \
                                            ajusteTemperaturaCristal *vectorPlano, _Bool ladoBombeo)
{
        propAstigmatica *salida=(propAstigmatica*)calloc(1,sizeof(propAstigmatica)); // Reserva de espacio
        long double deltaZ=L/(pasos-1); // Diferencia de distancia en la dirección de propagación
        long double wTan=spot_q(qInTan,n0,lambda0);
        long double wSag=spot_q(qInSag,n0,lambda0);

        long double _Complex *qPropTan = (long double _Complex*)calloc(pasos,sizeof(long double _Complex));
        long double _Complex *qPropSag = (long double _Complex*)calloc(pasos,sizeof(long double _Complex));
        long double _Complex At, Bt, Ct, Dt;
        long double _Complex As, Bs, Cs, Ds;
        assert(qPropTan);
        assert(qPropSag);
        matriz *MTan, *MSag;
        MTan=nuevaMatriz(2,2);
        MSag=nuevaMatriz(2,2);
        // Cálculo inicial
        qPropTan[0]=qInTan;
        qPropSag[0]=qInSag;
        long double n0Prop=n0;
        // Identificar si se propaga desde el lado bombeado o no, cambiar coeficientes
        // de ajuste acorde a esto.
        if(ladoBombeo==1)
        {
            for (int i=0;i<pasos;i++)
            {
                // Cálculo de spots
                long double wTan=spot_q(qPropTan[i],n0Prop,lambda0);
                long double wSag=spot_q(qPropSag[i],n0Prop,lambda0);
                // Cálculo de plano
                //n0Prop=n0+vectorPlano->a0[i]*dn_dv+n2*P_laser/(M_PI*wTan*wSag)*3.0/4.0; // Factores lineal, térmico y Kerr
                n0Prop=n0; // Descartados términos lineales por ser pequeños.
                // Formación de matrices
                // Factor de parábola=1/h^2
                long double parabola1=-2.0*vectorPlano->a1[i]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wTan,3)*wSag);
                long double parabola2=-2.0*vectorPlano->a2[i]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wSag,3)*wTan);
                // Coeficientes
                At=1;
                Bt=deltaZ/n0Prop;
                Ct=-n0Prop*parabola1*deltaZ;
                Dt=1;
                As=1;
                Bs=deltaZ/n0Prop;
                Cs=-n0Prop*parabola2*deltaZ;
                Ds=1;
                llenaMatrizSinCrear(MTan,At,Bt,Ct,Dt);
                llenaMatrizSinCrear(MSag,As,Bs,Cs,Ds);
                if(i+1<pasos)
                {
                    qPropTan[i+1]=prop_q(MTan,qPropTan[i],n0Prop,n0Prop);
                    qPropSag[i+1]=prop_q(MSag,qPropSag[i],n0Prop,n0Prop);
                }
            }
        }
        else
        {
            int indiceTermico=pasos-1;
            for (int i=0;i<pasos;i++)
            {
                // Cálculo de spots
                long double w_pump_t=spot_q(qPropTan[i],n0Prop,lambda0);
                long double w_pump_s=spot_q(qPropSag[i],n0Prop,lambda0);
                // Cálculo de plano
                n0Prop=n0;
                // Formación de matrices
                // Factor de parábola=1/h^2
                long double parabola1=-2.0*vectorPlano->a1[indiceTermico]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wTan,3)*wSag);
                long double parabola2=-2.0*vectorPlano->a2[indiceTermico]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wSag,3)*wTan);
                // Coeficientes
                At=1;
                Bt=deltaZ/n0Prop;
                Ct=-n0Prop*parabola1*deltaZ;
                Dt=1;
                As=1;
                Bs=deltaZ/n0Prop;
                Cs=-n0Prop*parabola2*deltaZ;
                Ds=1;
                llenaMatrizSinCrear(MTan,At,Bt,Ct,Dt);
                llenaMatrizSinCrear(MSag,As,Bs,Cs,Ds);
                if(i+1<pasos)
                {
                    qPropTan[i+1]=prop_q(MTan,qPropTan[i],n0Prop,n0Prop);
                    qPropSag[i+1]=prop_q(MSag,qPropSag[i],n0Prop,n0Prop);
                }
                indiceTermico--;
            }
        }
        salida->qTan=qPropTan[pasos-1];
        salida->qSag=qPropSag[pasos-1];
        //Limpieza
        borraMatriz(MTan);
        borraMatriz(MSag);
        free(qPropTan);
        free(qPropSag);
        return salida;
    }


// Propagación no lineal (matrices)
// Termina el lazo si alcanza el umbral solicitado
// Para epsilon1 y epsilon2

// EM1 - Busca que sean estables en ambos planos.

matriz ** propNoLinealEM1Astigmatico(matriz *qInTan, matriz *qInSag, ajusteTemperaturaCristal **vectorDeVectores, \
                                    char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, \
                                    long double L, long double f1, long double f2, long double n0, long double n2, \
                                    long double chi, long double kth, long double Cp, \
                                    long double rho,long double dn_dv,long double P_laser,\
                                    long double lambda0, matriz *epsilon1, matriz *epsilon2, int iteraciones, int pasosKerr, long double umbral, \
                                    _Bool guardaSpotIteracion, _Bool guardaVariacionIteracion, char* rutaSpotIterMax,\
                                    char* rutaVarIterTan, char* rutaVarIterSag, CtoCppQProgressBar *bar)
{
    // Variables
    //imprimeMatriz(q_in);
    CtoCppQProgressBar_setMaximum(bar,epsilon1->columnas*epsilon2->columnas);
    CtoCppQProgressBar_setMinimum(bar,0);
    int progreso=0;
    long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
    _Bool termino=0;
    _Bool escribeCSV=0;
    _Bool reallocFallido=0;
    long double *spotsTan, *spotsSag;
    int * iteracionActualTan, *iteracionActualSag;
    int cuenta=0;
    long double wIterTanViejo=0, wIterTanNuevo=0, wIterSagViejo=0, wIterSagNuevo=0, diferenciaTan=0, diferenciaSag=0; // Umbral no porcentual
    char* prefijoTan="EM1_tan";
    char* prefijoSag="EM1_sag";
    assert(strlen(conjugado_corto)==3);
    assert(strlen(conjugado_largo)==3);
    anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
    delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
    delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

    // Cálculo de distancias focales
    f1t=f1*cos(angulos[0]);
    f2t=f2*cos(angulos[1]);
    f1s=f1/cos(angulos[0]);
    f2s=f2/cos(angulos[1]);

    matriz ** q_new=(matriz**)calloc(2,sizeof(matriz*)); // Espacio para dos matrices
    q_new[0]=nuevaMatriz(qInTan->filas,qInTan->columnas); // Creando matriz tangencial
    q_new[1]=nuevaMatriz(qInSag->filas,qInSag->columnas); // Creando matriz tangencial
    long double _Complex *qPropTan=(long double _Complex*)calloc(6,sizeof(long double _Complex)); // Parámetros para propagación
    long double _Complex *qPropSag=(long double _Complex*)calloc(6,sizeof(long double _Complex));

    // Crea matrices para la propagación

    matriz *EM1tan1,*EM1tan2,*EM1tan3,*EM1tan4,*EM1tan5,*EM1tan6,*EM1tan7,*EM1tan8,*EM1tan9,*EM1tan10,*EM1tan11,*EM1tan12,*EM1tan13,*EM1tan14,*EM1tan15,*EM1tan16,*EM1tan17;
    matriz *EM1sag1,*EM1sag2,*EM1sag3,*EM1sag4,*EM1sag5,*EM1sag6,*EM1sag7,*EM1sag8,*EM1sag9,*EM1sag10,*EM1sag11,*EM1sag12,*EM1sag13,*EM1sag14,*EM1sag15,*EM1sag16,*EM1sag17;
    // EM1, Caso tangencial
    EM1tan1=llenaMatriz(1,L1,0,1);
    EM1tan2=llenaMatriz(1,0,-1/f1t,1);
    EM1tan3=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM1tan4=llenaMatriz(n0,0,0,1/n0);
    EM1tan5=llenaMatriz(1/n0,0,0,n0);
    EM1tan6=nuevaMatriz(2,2);//EM1tan5 pendiente
    EM1tan7=llenaMatriz(1,0,-1/f2t,1);
    EM1tan8=llenaMatriz(1,L2,0,1);
    EM1tan9=llenaMatriz(1,0,0,1);
    EM1tan10=llenaMatriz(1,L2,0,1);
    EM1tan11=llenaMatriz(1,0,-1/f2t,1);
    EM1tan12=nuevaMatriz(2,2);    //EM1tan11 pendiente
    EM1tan13=llenaMatriz(n0,0,0,1/n0);
    EM1tan14=llenaMatriz(1/n0,0,0,n0);
    EM1tan15=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM1tan16=llenaMatriz(1,0,-1/f1t,1);
    EM1tan17=llenaMatriz(1,L1,0,1);
    // EM1, Caso saggencial
    EM1sag1=llenaMatriz(1,L1,0,1);
    EM1sag2=llenaMatriz(1,0,-1/f1s,1);
    EM1sag3=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM1sag4=llenaMatriz(1,0,0,1);
    EM1sag5=llenaMatriz(1,0,0,1);
    EM1sag6=nuevaMatriz(2,2);//EM1sag5 pendiente
    EM1sag7=llenaMatriz(1,0,-1/f2s,1);
    EM1sag8=llenaMatriz(1,L2,0,1);
    EM1sag9=llenaMatriz(1,0,0,1);
    EM1sag10=llenaMatriz(1,L2,0,1);
    EM1sag11=llenaMatriz(1,0,-1/f2s,1);
    EM1sag12=nuevaMatriz(2,2);    //EM1sag11 pendiente
    EM1sag13=llenaMatriz(1,0,0,1);
    EM1sag14=llenaMatriz(1,0,0,1);
    EM1sag15=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM1sag16=llenaMatriz(1,0,-1/f1s,1);
    EM1sag17=llenaMatriz(1,L1,0,1);


    // Ciclo
    matriz * q_iter = nuevaMatriz(2,iteraciones); // 1 para tangencial, 2 para sagital
    for(int index1=1;index1<=epsilon1->columnas;index1++)
    {
        for(int index2=1;index2<=epsilon2->columnas;index2++)
        {
            spotsTan=(long double*)calloc(1,sizeof(long double));
            spotsSag=(long double*)calloc(1,sizeof(long double));
            iteracionActualTan=(int *)calloc(1,sizeof(int));
            iteracionActualSag=(int *)calloc(1,sizeof(int));
            cuenta=1;
            iteracionActualTan[0]=0;
            iteracionActualSag[0]=0;
            spotsTan[0]=spot_q(obtieneElemento(qInTan,index1,index2),1,lambda0);
            spotsSag[0]=spot_q(obtieneElemento(qInSag,index1,index2),1,lambda0);
            // Llenando matrices variables:
            llenaMatrizSinCrear(EM1tan3,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM1tan6,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM1tan12,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM1tan15,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM1sag3,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM1sag6,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM1sag12,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM1sag15,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);

            if((elem(qInTan,index1,index2)!=0)&&(elem(qInSag,index1,index2)!=0)) // Si no son cero, realiza rutina
            {
                fijaElemento(q_iter,1,1,obtieneElemento(qInTan,index1,index2));
                fijaElemento(q_iter,2,1,obtieneElemento(qInSag,index1,index2));
                escribeCSV=1;
                for(;;)
                {
                    int ahora=1;
                    int siguiente=2;
                    matriz * Prop1Tan_piezas[]={EM1tan4,EM1tan3,EM1tan2,EM1tan1};
                    matriz * Prop1Tan=variasMultMatriciales(Prop1Tan_piezas,4);
                    matriz * Prop1Sag_piezas[]={EM1sag4,EM1sag3,EM1sag2,EM1sag1};
                    matriz * Prop1Sag=variasMultMatriciales(Prop1Sag_piezas,4);
                    // Almacenando primer q
                    qPropTan[0]=obtieneElemento(q_iter,1,ahora);
                    qPropSag[0]=obtieneElemento(q_iter,2,ahora);
                    // Cálculo de q2
                    qPropTan[1]=prop_q(Prop1Tan,qPropTan[0],1,n0);
                    qPropSag[1]=prop_q(Prop1Sag,qPropSag[0],1,n0);
                    // Cálculo de q3
                    propAstigmatica *q3;
                    q3=propagacionKerrAcoplada(pasosKerr,L,n0,n2,qPropTan[1],qPropSag[1],chi,kth,Cp,rho,dn_dv,P_laser,lambda0,vectorDeVectores[index1-1],1);
                    qPropTan[2]=q3->qTan;
                    qPropSag[2]=q3->qSag;
                    free(q3);
                    // Más matrices
                    matriz * Prop2Tan_piezas[]={EM1tan13,EM1tan12,EM1tan11,EM1tan10,EM1tan9,EM1tan8,EM1tan7,EM1tan6,EM1tan5}; // 5 a 13
                    matriz * Prop2Tan=variasMultMatriciales(Prop2Tan_piezas,9);
                    matriz * Prop2Sag_piezas[]={EM1sag13,EM1sag12,EM1sag11,EM1sag10,EM1sag9,EM1sag8,EM1sag7,EM1sag6,EM1sag5}; // 5 a 13
                    matriz * Prop2Sag=variasMultMatriciales(Prop2Sag_piezas,9);
                    // Cálculo de q4
                    qPropTan[3]=prop_q(Prop2Tan,qPropTan[2],n0,n0);
                    qPropSag[3]=prop_q(Prop2Sag,qPropSag[2],n0,n0);
                    // Propagación de regreso Kerr, cálculo q5
                    propAstigmatica *q5=propagacionKerrAcoplada(pasosKerr,L,n0,n2,qPropTan[3],qPropSag[3],chi,kth,Cp,rho,dn_dv,P_laser,lambda0,vectorDeVectores[index1-1],0);
                    qPropTan[4]=q5->qTan;
                    qPropSag[4]=q5->qSag;
                    free(q5);
                    // Más matrices
                    matriz * Prop3Tan_piezas[]={EM1tan17,EM1tan16,EM1tan15,EM1tan14}; // 14 a 17
                    matriz * Prop3Tan=variasMultMatriciales(Prop3Tan_piezas,4);
                    matriz * Prop3Sag_piezas[]={EM1sag17,EM1sag16,EM1sag15,EM1sag14}; // 14 a 17
                    matriz * Prop3Sag=variasMultMatriciales(Prop3Sag_piezas,4);
                    // Cálculo de q6
                    qPropTan[5]=prop_q(Prop3Tan,qPropTan[4],n0,1);
                    qPropSag[5]=prop_q(Prop3Sag,qPropSag[4],n0,1);
                    // Guarda elementos
                    fijaElemento(q_iter,1,ahora,qPropTan[5]);
                    fijaElemento(q_iter,2,ahora,qPropSag[5]);

                    // Limpieza local
                    Prop1Tan=borraMatriz(Prop1Tan);
                    Prop2Tan=borraMatriz(Prop2Tan);
                    Prop3Tan=borraMatriz(Prop3Tan);
                    Prop1Sag=borraMatriz(Prop1Sag);
                    Prop2Sag=borraMatriz(Prop2Sag);
                    Prop3Sag=borraMatriz(Prop3Sag);

                    // Calculando spots y comparación
                    wIterTanViejo=spot_q(qPropTan[0],1,lambda0);
                    wIterTanNuevo=spot_q(qPropTan[5],1,lambda0);
                    wIterSagViejo=spot_q(qPropSag[0],1,lambda0);
                    wIterSagNuevo=spot_q(qPropSag[5],1,lambda0);

                    if(cuenta>iteraciones)  break; // Por si acaso
                    ++cuenta;
                    //printf("%i\n",cuenta);
                    void *tmp_sptTan, *tmp_sptSag, *tmp_iterTan, *tmp_iterSag;
                    tmp_sptTan=(long double*)realloc(spotsTan,cuenta*sizeof(long double));
                    tmp_sptSag=(long double*)realloc(spotsSag,cuenta*sizeof(long double));
                    tmp_iterTan=(int *)realloc(iteracionActualTan,cuenta*sizeof(int));
                    tmp_iterSag=(int *)realloc(iteracionActualSag,cuenta*sizeof(int));
                    if(tmp_sptTan!=NULL && tmp_sptSag!=NULL && tmp_iterTan!=NULL && tmp_iterSag!=NULL)
                    {
                        spotsTan=tmp_sptTan;
                        spotsSag=tmp_sptSag;
                        iteracionActualTan=tmp_iterTan;
                        iteracionActualSag=tmp_iterSag;
                        spotsTan[cuenta-1]=wIterTanNuevo;
                        spotsSag[cuenta-1]=wIterSagNuevo;
                        iteracionActualTan[cuenta-1]=cuenta-1;
                        iteracionActualSag[cuenta-1]=cuenta-1;
                    }
                    else
                    {
                        reallocFallido=1;
                        break;
                    }
                    if(wIterTanViejo>0.1||wIterTanNuevo>0.1||wIterSagViejo>0.1||wIterSagNuevo>0.1 ||\
                       isnan(wIterTanViejo)||isnan(wIterTanNuevo)||isnan(wIterSagViejo)||isnan(wIterSagNuevo))
                       {
                            termino=0;
                            escribeCSV=1;
                            fijaElemento(q_new[0],index1,index2,0);
                            fijaElemento(q_new[1],index1,index2,0);
                            break;
                       }
                    diferenciaTan=fabsl((wIterTanNuevo-wIterTanViejo)/wIterTanNuevo)*100.00;
                    diferenciaSag=fabsl((wIterSagNuevo-wIterSagViejo)/wIterSagNuevo)*100.00;
                    if(diferenciaTan<=umbral&&diferenciaSag<=umbral) // Termina la iteración
                    {
                        fijaElemento(q_new[0],index1,index2,obtieneElemento(q_iter,1,ahora));
                        fijaElemento(q_new[1],index1,index2,obtieneElemento(q_iter,2,ahora));
                        termino=1;
                        escribeCSV=1;
                        break;
                    }
                    else
                    {
                        if(diferenciaTan<=umbral)
                            fijaElemento(q_new[0],index1,index2,obtieneElemento(q_iter,1,ahora));
                        else
                            fijaElemento(q_new[0],index1,index2,0);
                        if(diferenciaSag<=umbral)
                            fijaElemento(q_new[1],index1,index2,obtieneElemento(q_iter,2,ahora));
                        else
                            fijaElemento(q_new[1],index1,index2,0);
                        termino=0;
                    }
                }
            }
            else
            {
                if(elem(qInTan,index1,index2)==0)
                    fijaElemento(q_new[0],index1,index2,0.0);
                if(elem(qInSag,index1,index2)==0)
                    fijaElemento(q_new[1],index1,index2,0.0);
                escribeCSV=0;
                reallocFallido=0;
            }
            if(reallocFallido==1)
            {
                termino=0;
                free(spotsTan);
                free(spotsSag);
                free(iteracionActualTan);
                free(iteracionActualSag);
                printf("Se acabó la memoria\n");
                reallocFallido=0;
            }
            if(guardaSpotIteracion)
                {
                    guardaSpotIteraciones(spotsTan[cuenta-1],spotsSag[cuenta-1],cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon1,1,index2)),rutaSpotIterMax,termino);
                }
            if(escribeCSV==1&&guardaVariacionIteracion)
            {

                guardaVariaciones(spotsTan,iteracionActualTan,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),prefijoTan,rutaVarIterTan,termino);
                guardaVariaciones(spotsSag,iteracionActualSag,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),prefijoSag,rutaVarIterSag,termino);
            }
            progreso++;
            CtoCppQProgressBar_setValue(bar,progreso);
            if(termino==0)
                printf("Cavidad en X, EM1: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo tangencial: %Le\nError relativo sagital: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferenciaTan,diferenciaSag);
            else
                printf("Cavidad en X, EM1: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo tangencial: %Le\nError relativo sagital: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferenciaTan,diferenciaSag);
            free(spotsTan);
            free(spotsSag);
            free(iteracionActualTan);
            free(iteracionActualSag);
        }
    }
    // Limpieza
    borraMatriz(q_iter);
    borraMatriz(EM1tan1);
    borraMatriz(EM1tan2);
    borraMatriz(EM1tan3);
    borraMatriz(EM1tan4);
    borraMatriz(EM1tan5);
    borraMatriz(EM1tan6);
    borraMatriz(EM1tan7);
    borraMatriz(EM1tan8);
    borraMatriz(EM1tan9);
    borraMatriz(EM1tan10);
    borraMatriz(EM1tan11);
    borraMatriz(EM1tan12);
    borraMatriz(EM1tan13);
    borraMatriz(EM1tan14);
    borraMatriz(EM1tan15);
    borraMatriz(EM1tan16);
    borraMatriz(EM1tan17);
    borraMatriz(EM1sag1);
    borraMatriz(EM1sag2);
    borraMatriz(EM1sag3);
    borraMatriz(EM1sag4);
    borraMatriz(EM1sag5);
    borraMatriz(EM1sag6);
    borraMatriz(EM1sag7);
    borraMatriz(EM1sag8);
    borraMatriz(EM1sag9);
    borraMatriz(EM1sag10);
    borraMatriz(EM1sag11);
    borraMatriz(EM1sag12);
    borraMatriz(EM1sag13);
    borraMatriz(EM1sag14);
    borraMatriz(EM1sag15);
    borraMatriz(EM1sag16);
    borraMatriz(EM1sag17);
    return q_new;
}

matriz ** propNoLinealEM2Astigmatico(matriz *qInTan, matriz *qInSag, ajusteTemperaturaCristal **vectorDeVectores, \
                                    char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, \
                                    long double L, long double f1, long double f2, long double n0, long double n2, \
                                    long double chi, long double kth, long double Cp, \
                                    long double rho,long double dn_dv,long double P_laser,\
                                    long double lambda0, matriz *epsilon1, matriz *epsilon2, int iteraciones, int pasosKerr, long double umbral, \
                                    _Bool guardaSpotIteracion, _Bool guardaVariacionIteracion, char* rutaSpotIterMax,\
                                    char* rutaVarIterTan, char* rutaVarIterSag, CtoCppQProgressBar *bar)
{
    // Variables
    //imprimeMatriz(q_in);
    CtoCppQProgressBar_setMaximum(bar,epsilon1->columnas*epsilon2->columnas);
    CtoCppQProgressBar_setMinimum(bar,0);
    int progreso=0;
    long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
    _Bool termino=0;
    _Bool escribeCSV=0;
    _Bool reallocFallido=0;
    long double *spotsTan, *spotsSag;
    int * iteracionActualTan, *iteracionActualSag;
    int cuenta=0;
    long double wIterTanViejo=0, wIterTanNuevo=0, wIterSagViejo=0, wIterSagNuevo=0, diferenciaTan=0, diferenciaSag=0; // Umbral no porcentual
    char* prefijoTan="EM2_tan";
    char* prefijoSag="EM2_sag";
    assert(strlen(conjugado_corto)==3);
    assert(strlen(conjugado_largo)==3);
    anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
    delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
    delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

    // Cálculo de distancias focales
    f1t=f1*cos(angulos[0]);
    f2t=f2*cos(angulos[1]);
    f1s=f1/cos(angulos[0]);
    f2s=f2/cos(angulos[1]);

    matriz ** q_new=(matriz**)calloc(2,sizeof(matriz*)); // Espacio para dos matrices
    q_new[0]=nuevaMatriz(qInTan->filas,qInTan->columnas); // Creando matriz tangencial
    q_new[1]=nuevaMatriz(qInSag->filas,qInSag->columnas); // Creando matriz tangencial
    long double _Complex *qPropTan=(long double _Complex*)calloc(6,sizeof(long double _Complex)); // Parámetros para propagación
    long double _Complex *qPropSag=(long double _Complex*)calloc(6,sizeof(long double _Complex));

    // Crea matrices para la propagación

    matriz *EM2tan1,*EM2tan2,*EM2tan3,*EM2tan4,*EM2tan5,*EM2tan6,*EM2tan7,*EM2tan8,*EM2tan9,*EM2tan10,*EM2tan11,*EM2tan12,*EM2tan13,*EM2tan14,*EM2tan15,*EM2tan16,*EM2tan17;
    matriz *EM2sag1,*EM2sag2,*EM2sag3,*EM2sag4,*EM2sag5,*EM2sag6,*EM2sag7,*EM2sag8,*EM2sag9,*EM2sag10,*EM2sag11,*EM2sag12,*EM2sag13,*EM2sag14,*EM2sag15,*EM2sag16,*EM2sag17;

    // EM2, Caso tangencial
    EM2tan1=llenaMatriz(1,L2,0,1);
    EM2tan2=llenaMatriz(1,0,-1/f2t,1);
    EM2tan3=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM2tan4=llenaMatriz(n0,0,0,1/n0);
    EM2tan5=llenaMatriz(1/n0,0,0,n0);
    EM2tan6=nuevaMatriz(2,2);//EM2tan5 pendiente
    EM2tan7=llenaMatriz(1,0,-1/f1t,1);
    EM2tan8=llenaMatriz(1,L1,0,1);
    EM2tan9=llenaMatriz(1,0,0,1);
    EM2tan10=llenaMatriz(1,L1,0,1);
    EM2tan11=llenaMatriz(1,0,-1/f1t,1);
    EM2tan12=nuevaMatriz(2,2);    //EM2tan11 pendiente
    EM2tan13=llenaMatriz(n0,0,0,1/n0);
    EM2tan14=llenaMatriz(1/n0,0,0,n0);
    EM2tan15=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM2tan16=llenaMatriz(1,0,-1/f2t,1);
    EM2tan17=llenaMatriz(1,L2,0,1);
    // EM1, Caso saggencial
    EM2sag1=llenaMatriz(1,L2,0,1);
    EM2sag2=llenaMatriz(1,0,-1/f2s,1);
    EM2sag3=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM2sag4=llenaMatriz(1,0,0,1);
    EM2sag5=llenaMatriz(1,0,0,1);
    EM2sag6=nuevaMatriz(2,2);//EM2sag5 pendiente
    EM2sag7=llenaMatriz(1,0,-1/f1s,1);
    EM2sag8=llenaMatriz(1,L1,0,1);
    EM2sag9=llenaMatriz(1,0,0,1);
    EM2sag10=llenaMatriz(1,L1,0,1);
    EM2sag11=llenaMatriz(1,0,-1/f1s,1);
    EM2sag12=nuevaMatriz(2,2);    //EM2sag11 pendiente
    EM2sag13=llenaMatriz(1,0,0,1);
    EM2sag14=llenaMatriz(1,0,0,1);
    EM2sag15=nuevaMatriz(2,2); //llenaMatriz(1,delta1,0,1);
    EM2sag16=llenaMatriz(1,0,-1/f2s,1);
    EM2sag17=llenaMatriz(1,L2,0,1);


    // Ciclo
    matriz * q_iter = nuevaMatriz(2,iteraciones); // 1 para tangencial, 2 para sagital
    for(int index1=1;index1<=epsilon1->columnas;index1++)
    {
        for(int index2=1;index2<=epsilon2->columnas;index2++)
        {
            spotsTan=(long double*)calloc(1,sizeof(long double));
            spotsSag=(long double*)calloc(1,sizeof(long double));
            iteracionActualTan=(int *)calloc(1,sizeof(int));
            iteracionActualSag=(int *)calloc(1,sizeof(int));
            cuenta=1;
            iteracionActualTan[0]=0;
            iteracionActualSag[0]=0;
            spotsTan[0]=spot_q(obtieneElemento(qInTan,index1,index2),1,lambda0);
            spotsSag[0]=spot_q(obtieneElemento(qInSag,index1,index2),1,lambda0);
            // Llenando matrices variables:
            llenaMatrizSinCrear(EM2tan3,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM2tan6,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM2tan12,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM2tan15,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM2sag3,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
            llenaMatrizSinCrear(EM2sag6,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM2sag12,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
            llenaMatrizSinCrear(EM2sag15,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);

            if((elem(qInTan,index1,index2)!=0)&&(elem(qInSag,index1,index2)!=0)) // Si no son cero, realiza rutina
            {
                fijaElemento(q_iter,1,1,obtieneElemento(qInTan,index1,index2));
                fijaElemento(q_iter,2,1,obtieneElemento(qInSag,index1,index2));
                spotsTan=(long double*)calloc(1,sizeof(long double));
                spotsSag=(long double*)calloc(1,sizeof(long double));
                iteracionActualTan=(int *)calloc(1,sizeof(int));
                iteracionActualSag=(int *)calloc(1,sizeof(int));
                cuenta=1;
                iteracionActualTan[0]=0;
                iteracionActualSag[0]=0;
                spotsTan[0]=spot_q(obtieneElemento(qInTan,index1,index2),1,lambda0);
                spotsSag[0]=spot_q(obtieneElemento(qInSag,index1,index2),1,lambda0);
                escribeCSV=1;
                for(;;)
                {
                    int ahora=1;
                    int siguiente=2;
                    matriz * Prop1Tan_piezas[]={EM2tan4,EM2tan3,EM2tan2,EM2tan1};
                    matriz * Prop1Tan=variasMultMatriciales(Prop1Tan_piezas,4);
                    matriz * Prop1Sag_piezas[]={EM2sag4,EM2sag3,EM2sag2,EM2sag1};
                    matriz * Prop1Sag=variasMultMatriciales(Prop1Sag_piezas,4);
                    // Almacenando primer q
                    qPropTan[0]=obtieneElemento(q_iter,1,ahora);
                    qPropSag[0]=obtieneElemento(q_iter,2,ahora);
                    // Cálculo de q2
                    qPropTan[1]=prop_q(Prop1Tan,qPropTan[0],1,n0);
                    qPropSag[1]=prop_q(Prop1Sag,qPropSag[0],1,n0);
                    // Cálculo de q3
                    propAstigmatica *q3=propagacionKerrAcoplada(pasosKerr,L,n0,n2,qPropTan[1],qPropSag[1],chi,kth,Cp,rho,dn_dv,P_laser,lambda0,vectorDeVectores[index1-1],0);
                    qPropTan[2]=q3->qTan;
                    qPropSag[2]=q3->qSag;
                    free(q3);
                    // Más matrices
                    matriz * Prop2Tan_piezas[]={EM2tan13,EM2tan12,EM2tan11,EM2tan10,EM2tan9,EM2tan8,EM2tan7,EM2tan6,EM2tan5}; // 5 a 13
                    matriz * Prop2Tan=variasMultMatriciales(Prop2Tan_piezas,9);
                    matriz * Prop2Sag_piezas[]={EM2sag13,EM2sag12,EM2sag11,EM2sag10,EM2sag9,EM2sag8,EM2sag7,EM2sag6,EM2sag5}; // 5 a 13
                    matriz * Prop2Sag=variasMultMatriciales(Prop2Sag_piezas,9);
                    // Cálculo de q4
                    qPropTan[3]=prop_q(Prop2Tan,qPropTan[2],n0,n0);
                    qPropSag[3]=prop_q(Prop2Sag,qPropSag[2],n0,n0);
                    // Propagación de regreso Kerr, cálculo q5
                    propAstigmatica *q5=propagacionKerrAcoplada(pasosKerr,L,n0,n2,qPropTan[3],qPropSag[3],chi,kth,Cp,rho,dn_dv,P_laser,lambda0,vectorDeVectores[index1-1],1);
                    qPropTan[4]=q5->qTan;
                    qPropSag[4]=q5->qSag;
                    free(q5);
                    // Más matrices
                    matriz * Prop3Tan_piezas[]={EM2tan17,EM2tan16,EM2tan15,EM2tan14}; // 14 a 17
                    matriz * Prop3Tan=variasMultMatriciales(Prop3Tan_piezas,4);
                    matriz * Prop3Sag_piezas[]={EM2sag17,EM2sag16,EM2sag15,EM2sag14}; // 14 a 17
                    matriz * Prop3Sag=variasMultMatriciales(Prop3Sag_piezas,4);
                    // Cálculo de q6
                    qPropTan[5]=prop_q(Prop3Tan,qPropTan[4],n0,1);
                    qPropSag[5]=prop_q(Prop3Sag,qPropSag[4],n0,1);
                    // Guarda elementos
                    fijaElemento(q_iter,1,ahora,qPropTan[5]);
                    fijaElemento(q_iter,2,ahora,qPropSag[5]);

                    // Limpieza local
                    Prop1Tan=borraMatriz(Prop1Tan);
                    Prop2Tan=borraMatriz(Prop2Tan);
                    Prop3Tan=borraMatriz(Prop3Tan);
                    Prop1Sag=borraMatriz(Prop1Sag);
                    Prop2Sag=borraMatriz(Prop2Sag);
                    Prop3Sag=borraMatriz(Prop3Sag);

                    // Calculando spots y comparación
                    wIterTanViejo=spot_q(qPropTan[0],1,lambda0);
                    wIterTanNuevo=spot_q(qPropTan[5],1,lambda0);
                    wIterSagViejo=spot_q(qPropSag[0],1,lambda0);
                    wIterSagNuevo=spot_q(qPropSag[5],1,lambda0);

                    if(cuenta>iteraciones)  break; // Por si acaso
                    ++cuenta;
                    //printf("%i\n",cuenta);
                    void *tmp_sptTan, *tmp_sptSag, *tmp_iterTan, *tmp_iterSag;
                    tmp_sptTan=(long double*)realloc(spotsTan,cuenta*sizeof(long double));
                    tmp_sptSag=(long double*)realloc(spotsSag,cuenta*sizeof(long double));
                    tmp_iterTan=(int *)realloc(iteracionActualTan,cuenta*sizeof(int));
                    tmp_iterSag=(int *)realloc(iteracionActualSag,cuenta*sizeof(int));
                    if(tmp_sptTan!=NULL && tmp_sptSag!=NULL && tmp_iterTan!=NULL && tmp_iterSag!=NULL)
                    {
                        spotsTan=tmp_sptTan;
                        spotsSag=tmp_sptSag;
                        iteracionActualTan=tmp_iterTan;
                        iteracionActualSag=tmp_iterSag;
                        spotsTan[cuenta-1]=wIterTanNuevo;
                        spotsSag[cuenta-1]=wIterSagNuevo;
                        iteracionActualTan[cuenta-1]=cuenta-1;
                        iteracionActualSag[cuenta-1]=cuenta-1;
                    }
                    else
                    {
                        reallocFallido=1;
                        break;
                    }
                    if(wIterTanViejo>0.1||wIterTanNuevo>0.1||wIterSagViejo>0.1||wIterSagNuevo>0.1 ||\
                       isnan(wIterTanViejo)||isnan(wIterTanNuevo)||isnan(wIterSagViejo)||isnan(wIterSagNuevo))
                       {
                            termino=0;
                            escribeCSV=1;
                            fijaElemento(q_new[0],index1,index2,0);
                            fijaElemento(q_new[1],index1,index2,0);
                            break;
                       }
                    diferenciaTan=fabsl((wIterTanNuevo-wIterTanViejo)/wIterTanNuevo)*100.00;
                    diferenciaSag=fabsl((wIterSagNuevo-wIterSagViejo)/wIterSagNuevo)*100.00;
                    if(diferenciaTan<=umbral&&diferenciaSag<=umbral) // Termina la iteración
                    {
                        fijaElemento(q_new[0],index1,index2,obtieneElemento(q_iter,1,ahora));
                        fijaElemento(q_new[1],index1,index2,obtieneElemento(q_iter,2,ahora));
                        termino=1;
                        escribeCSV=1;
                        break;
                    }
                    else
                    {
                        if(diferenciaTan<=umbral)
                            fijaElemento(q_new[0],index1,index2,obtieneElemento(q_iter,1,ahora));
                        else
                            fijaElemento(q_new[0],index1,index2,0);
                        if(diferenciaSag<=umbral)
                            fijaElemento(q_new[1],index1,index2,obtieneElemento(q_iter,2,ahora));
                        else
                            fijaElemento(q_new[1],index1,index2,0);
                        termino=0;
                    }
                }
            }
            else
            {
                if(elem(qInTan,index1,index2)==0)
                    fijaElemento(q_new[0],index1,index2,0.0);
                if(elem(qInSag,index1,index2)==0)
                    fijaElemento(q_new[1],index1,index2,0.0);
                escribeCSV=0;
                reallocFallido=0;
            }
            if(reallocFallido==1)
            {
                termino=0;
                free(spotsTan);
                free(spotsSag);
                free(iteracionActualTan);
                free(iteracionActualSag);
                printf("Se acabó la memoria\n");
                reallocFallido=0;
            }
            if(guardaSpotIteracion)
                {
                    guardaSpotIteraciones(spotsTan[cuenta-1],spotsSag[cuenta-1],cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon1,1,index2)),rutaSpotIterMax,termino);
                }
            if(escribeCSV==1&&guardaVariacionIteracion)
            {

                guardaVariaciones(spotsTan,iteracionActualTan,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),prefijoTan,rutaVarIterTan,termino);
                guardaVariaciones(spotsSag,iteracionActualSag,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),prefijoSag,rutaVarIterSag,termino);
            }
            progreso++;
            CtoCppQProgressBar_setValue(bar,progreso);
            if(termino==0)
                printf("Cavidad en X, EM2: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo tangencial: %Le\nError relativo sagital: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferenciaTan,diferenciaSag);
            else
                printf("Cavidad en X, EM2: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo tangencial: %Le\nError relativo sagital: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferenciaTan,diferenciaSag);
            free(spotsTan);
            free(spotsSag);
            free(iteracionActualTan);
            free(iteracionActualSag);
        }
    }
    // Limpieza
    free(qPropTan);
    free(qPropSag);
    borraMatriz(q_iter);
    borraMatriz(EM2tan1);
    borraMatriz(EM2tan2);
    borraMatriz(EM2tan3);
    borraMatriz(EM2tan4);
    borraMatriz(EM2tan5);
    borraMatriz(EM2tan6);
    borraMatriz(EM2tan7);
    borraMatriz(EM2tan8);
    borraMatriz(EM2tan9);
    borraMatriz(EM2tan10);
    borraMatriz(EM2tan11);
    borraMatriz(EM2tan12);
    borraMatriz(EM2tan13);
    borraMatriz(EM2tan14);
    borraMatriz(EM2tan15);
    borraMatriz(EM2tan16);
    borraMatriz(EM2tan17);
    borraMatriz(EM2sag1);
    borraMatriz(EM2sag2);
    borraMatriz(EM2sag3);
    borraMatriz(EM2sag4);
    borraMatriz(EM2sag5);
    borraMatriz(EM2sag6);
    borraMatriz(EM2sag7);
    borraMatriz(EM2sag8);
    borraMatriz(EM2sag9);
    borraMatriz(EM2sag10);
    borraMatriz(EM2sag11);
    borraMatriz(EM2sag12);
    borraMatriz(EM2sag13);
    borraMatriz(EM2sag14);
    borraMatriz(EM2sag15);
    borraMatriz(EM2sag16);
    borraMatriz(EM2sag17);
    return q_new;
}

