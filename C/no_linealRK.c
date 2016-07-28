     #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <complex.h>
    #include <string.h>
    #include <assert.h>
    #include "matrices.h"
    #include "lineal.h"
    #include "no_lineal.h"
    #include "error_iteraciones.h"

    int pasosKerrRK=1000;
    // Resolver el sistema de ecuaciones diferenciales acopladas mostradas en el artículo de D. Huang.
    // Emplea Runge-Kutta de 4to orden. Se define a p=1/q, por lo que pReal=Re(p) y pImag=Im(p).
    // xi queda definido por sqrt(1-P_laser/P_cr). Registrar la posición al momento de evaluar

    // Ecuaciones diferenciales

    long double dpReal(long double pReal, long double pImag, long double n0, long double n2, long double lambda0, long double P_laser) //
    {
        long double evaluacion;
        //evaluacion=-powl(pReal,2)+powl(pImag,2)*(1-(M_PI*n0*8*n2*P_laser)/lambda0);
        evaluacion=-powl(pReal,2)+powl(pImag,2)*(1-(M_PI*n0*n2*P_laser)/(2*powl(lambda0,2))); // Ajuste cuadrático
        return evaluacion;
    }

    long double dpImag(long double pReal, long double pImag, long double n0)
    {
        long double evaluacion;
        evaluacion=-2*pReal*pImag;
        return evaluacion;
    }

    // Runge-Kutta 4to orden aplicado a ec. diferenciales de Huang, h es el tamaño de paso.

    long double * pasoKerrRK(long double h, long double pReal_ini, long double pImag_ini, long double n0, long double n2, long double lambda0, long double P_laser)
    {
        long double k[2][4], pReal, pReal_out, pImag, pImag_out;
        long double *salida=(long double *)calloc(2,sizeof(long double));
        memset(k,0,sizeof k);
        pReal=pReal_ini;
        pImag=pImag_ini;
        k[0][0]=dpReal(pReal,pImag,n0,n2,lambda0,P_laser);
        k[1][0]=dpImag(pReal,pImag,n0);
        pReal=pReal_ini+k[0][0]*h/2;
        pImag=pImag_ini+k[1][0]*h/2;
        k[0][1]=dpReal(pReal,pImag,n0,n2,lambda0,P_laser);
        k[1][1]=dpImag(pReal,pImag,n0);
        pReal=pReal_ini+k[0][1]*h/2;
        pImag=pImag_ini+k[1][1]*h/2;
        k[0][2]=dpReal(pReal,pImag,n0,n2,lambda0,P_laser);
        k[1][2]=dpImag(pReal,pImag,n0);
        pReal=pReal_ini+k[0][2]*h;
        pImag=pImag_ini+k[1][2]*h;
        k[0][3]=dpReal(pReal,pImag,n0,n2,lambda0,P_laser);
        k[1][3]=dpImag(pReal,pImag,n0);
        pReal_out=pReal_ini+(k[0][0]+2*(k[0][1]+k[0][2])+k[0][3])*h/6;
        pImag_out=pImag_ini+(k[1][0]+2*(k[1][1]+k[1][2])+k[1][3])*h/6;
        salida[0]=pReal_out;
        salida[1]=pImag_out;
        return salida;
    }

    long double _Complex propagacionKerrRK(long double _Complex q_in, long double lambda0, long double L, long double n0, long double n2, long double P_laser, int pasos)
    {
        long double h=L/pasos;
        long double _Complex p_out, q_out;
        long double *pReal = (long double *)calloc(pasos,sizeof(long double));
        long double *pImag = (long double *)calloc(pasos,sizeof(long double));
        long double *datosPaso;
        //long double *datosPaso = (long double *)calloc(2,sizeof(long double));
        assert(pReal);
        assert(pImag);
        int i=pasos-1;
        pReal[i]=creall(1/q_in);
        pImag[i]=cimagl(1/q_in);
        for(i=pasos-1;i>0;i--)
        {
            datosPaso=pasoKerrRK(h,pReal[i],pImag[i],n0,n2,lambda0,P_laser);
            pReal[i-1]=datosPaso[0];
            pImag[i-1]=datosPaso[1];
            free(datosPaso);
        }
        p_out=pReal[0]+I*pImag[0];
        q_out=1/p_out;
        free(pReal);
        free(pImag);
        return q_out;
    }

    matriz * propNoLinealRK_tanEM1(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        _Bool escribeCSV=0;
        _Bool reallocFallido=0;
        long double *spots;
        int * iteracionActual;
        int cuenta=0;
        char tipo[]="EM1TanRK";
        char subCarpeta[]="EM1Tan_RK";
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4,q_5,q_6; // Parámetros complejos para propagación

        // Crea matrices para la propagación

        matriz *EM1tan1,*EM1tan2,*EM1tan3,*EM1tan4,*EM1tan5,*EM1tan6,*EM1tan7,*EM1tan8,*EM1tan9,*EM1tan10,*EM1tan11,*EM1tan12,*EM1tan13,*EM1tan14,*EM1tan15,*EM1tan16,*EM1tan17;
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

        // Ciclo
        matriz * q_iter = nuevaMatriz(1,2);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                spots=(long double*)calloc(1,sizeof(long double));
                iteracionActual=(int*)calloc(1,sizeof(int));
                cuenta=1;
                iteracionActual[0]=0;
                spots[0]=spot_q(obtieneElemento(q_in,index1,index2),1,lambda0);
                // Llenando matrices variables:
                llenaMatrizSinCrear(EM1tan3,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(EM1tan6,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(EM1tan12,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(EM1tan15,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);

                if(elem(q_in,index1,index2)==0)
                       {elem(q_new,index1,index2)=0;
                    escribeCSV=0;
                    reallocFallido=0;
                    }
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    escribeCSV=1;                //for(int iter=1;iter<iteraciones;iter++)
                   for(;;)
                    {
                        int ahora=1;
                        int siguiente=2;
                        matriz * Prop1_piezas[]={EM1tan4,EM1tan3,EM1tan2,EM1tan1};
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        //q_2=prop_q(Prop1,obtieneElemento(q_iter,1,ahora),1,n0);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={EM1tan13,EM1tan12,EM1tan11,EM1tan10,EM1tan9,EM1tan8,EM1tan7,EM1tan6,EM1tan5}; // 5 a 13
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,9);
                        q_4=prop_q(Prop2,q_3,n0,n0);
                        //q_5=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_4,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_5=propagacionKerrRK(q_4,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop3_piezas[]={EM1tan17,EM1tan16,EM1tan15,EM1tan14}; // 14 a 17
                        matriz * Prop3=variasMultMatriciales(Prop3_piezas,4);
                        q_6=prop_q(Prop3,q_5,n0,1);
                        fijaElemento(q_iter,1,ahora,q_6);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieaza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        Prop3=borraMatriz(Prop3);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_6,1,lambda0);// Prepara vectores
                        if(cuenta>iteraciones)  break; // Por si acaso
                        ++cuenta;
                        //printf("%i\n",cuenta);
                        void *tmp_spt, *tmp_iter;
                        tmp_spt=(long double*)realloc(spots,cuenta*sizeof(long double));
                        tmp_iter=(int *)realloc(iteracionActual,cuenta*sizeof(int));
                        if(tmp_spt!=NULL && tmp_iter!=NULL)
                        {
                            spots=tmp_spt;
                            iteracionActual=tmp_iter;
                            spots[cuenta-1]=w_iter_nuevo;
                            iteracionActual[cuenta-1]=cuenta-1;
                        }
                        else
                        {
                            reallocFallido=1;
                            break;
                        }
                        if(w_iter_viejo>0.1||w_iter_nuevo>0.1||isnan(w_iter_viejo) ||isnan(w_iter_nuevo))
                        {
                            termino=0;
                            escribeCSV=1;
                            fijaElemento(q_new,index1,index2,0);
                            break;
                        }  // Si el spot mide 20 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo)*100.00;
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            //fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,ahora));
                            long double _Complex valor=obtieneElemento(q_iter,1,1);
                            fijaElemento(q_new,index1,index2,valor);
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(reallocFallido==1)
                {
                    termino=0;
                    free(spots);
                    free(iteracionActual);
                    printf("Se acabó la memoria\n");
                    reallocFallido=0;
                }
                guardaSpotIteraciones(spots[cuenta-1],cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,termino);
                if(escribeCSV==1)
                    guardaVariaciones(spots,iteracionActual,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,subCarpeta,termino);
                if(termino==0)
                    printf("Cavidad en X, EM1tanRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en X, EM1tanRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                free(spots);
                free(iteracionActual);
            }
        }
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
        return q_new;
    }

    // EM1 sagital

    matriz * propNoLinealRK_sagEM1(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        _Bool escribeCSV=0;
        _Bool reallocFallido=0;
        long double *spots;
        int * iteracionActual;
        int cuenta=0;
        char tipo[]="EM1SagRK";
        char subCarpeta[]="EM1Sag_RK";
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4,q_5,q_6; // Parámetros complejos para propagación

        // Crea matrices para la propagación

        matriz *EM1sag1,*EM1sag2,*EM1sag3,*EM1sag4,*EM1sag5,*EM1sag6,*EM1sag7,*EM1sag8,*EM1sag9,*EM1sag10,*EM1sag11,*EM1sag12,*EM1sag13,*EM1sag14,*EM1sag15,*EM1sag16,*EM1sag17;
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
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                spots=(long double*)calloc(1,sizeof(long double));
                iteracionActual=(int*)calloc(1,sizeof(int));
                cuenta=1;
                iteracionActual[0]=0;
                spots[0]=spot_q(obtieneElemento(q_in,index1,index2),1,lambda0);
                // Llenando matrices variables:
                llenaMatrizSinCrear(EM1sag3,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(EM1sag6,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(EM1sag12,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(EM1sag15,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);

                if(elem(q_in,index1,index2)==0)
                    {elem(q_new,index1,index2)=0;
                    escribeCSV=0;
                    reallocFallido=0;
                    }
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    escribeCSV=1;                //for(int iter=1;iter<iteraciones;iter++)
                     for(;;)
                    {
                        int ahora=1;
                        int siguiente=2;
                        matriz * Prop1_piezas[]={EM1sag4,EM1sag3,EM1sag2,EM1sag1};
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        //q_2=prop_q(Prop1,obtieneElemento(q_iter,1,ahora),1,n0);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={EM1sag13,EM1sag12,EM1sag11,EM1sag10,EM1sag9,EM1sag8,EM1sag7,EM1sag6,EM1sag5}; // 5 a 13
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,9);
                        q_4=prop_q(Prop2,q_3,n0,n0);
                        //q_5=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_4,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_5=propagacionKerrRK(q_4,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop3_piezas[]={EM1sag17,EM1sag16,EM1sag15,EM1sag14}; // 14 a 17
                        matriz * Prop3=variasMultMatriciales(Prop3_piezas,4);
                        q_6=prop_q(Prop3,q_5,n0,1);
                        fijaElemento(q_iter,1,ahora,q_6);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieaza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        Prop3=borraMatriz(Prop3);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_6,1,lambda0);// Prepara vectores
                        if(cuenta>iteraciones)  break; // Por si acaso
                        ++cuenta;
                        //printf("%i\n",cuenta);
                        void *tmp_spt, *tmp_iter;
                        tmp_spt=(long double*)realloc(spots,cuenta*sizeof(long double));
                        tmp_iter=(int *)realloc(iteracionActual,cuenta*sizeof(int));
                        if(tmp_spt!=NULL && tmp_iter!=NULL)
                        {
                            spots=tmp_spt;
                            iteracionActual=tmp_iter;
                            spots[cuenta-1]=w_iter_nuevo;
                            iteracionActual[cuenta-1]=cuenta-1;
                        }
                        else
                        {
                            reallocFallido=1;
                            break;
                        }
                        if(w_iter_viejo>0.1||w_iter_nuevo>0.1||isnan(w_iter_viejo) ||isnan(w_iter_nuevo))
                        {
                            termino=0;
                            escribeCSV=1;
                            fijaElemento(q_new,index1,index2,0);
                            break;
                        }  // Si el spot mide 20 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo)*100.00;
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,ahora));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(reallocFallido==1)
                {
                    termino=0;
                    free(spots);
                    free(iteracionActual);
                    printf("Se acabó la memoria\n");
                    reallocFallido=0;
                }
                guardaSpotIteraciones(spots[cuenta-1],cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,termino);
                if(escribeCSV==1)
                    guardaVariaciones(spots,iteracionActual,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,subCarpeta,termino);
                if(termino==0)
                    printf("Cavidad en X, EM1sagRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en X, EM1sagRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                free(spots);
                free(iteracionActual);
            }
        }
        borraMatriz(q_iter);
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

    // EM2 tangencial
    matriz * propNoLinealRK_tanEM2(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        _Bool escribeCSV=0;
        _Bool reallocFallido=0;
        long double *spots;
        int * iteracionActual;
        int cuenta=0;
        char tipo[]="EM2TanRK";
        char subCarpeta[]="EM2Tan_RK";
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4,q_5,q_6; // Parámetros complejos para propagación

        // Crea matrices para la propagación

        matriz *EM2tan1,*EM2tan2,*EM2tan3,*EM2tan4,*EM2tan5,*EM2tan6,*EM2tan7,*EM2tan8,*EM2tan9,*EM2tan10,*EM2tan11,*EM2tan12,*EM2tan13,*EM2tan14,*EM2tan15,*EM2tan16,*EM2tan17;
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

        // Ciclo
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                spots=(long double*)calloc(1,sizeof(long double));
                iteracionActual=(int*)calloc(1,sizeof(int));
                cuenta=1;
                iteracionActual[0]=0;
                spots[0]=spot_q(obtieneElemento(q_in,index1,index2),1,lambda0);
                // Llenando matrices variables:
                llenaMatrizSinCrear(EM2tan3,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(EM2tan6,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(EM2tan12,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(EM2tan15,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);

                if(elem(q_in,index1,index2)==0)
                {elem(q_new,index1,index2)=0;
                    escribeCSV=0;
                    reallocFallido=0;
                    }
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    escribeCSV=1;                //for(int iter=1;iter<iteraciones;iter++)
                    for(;;)
                    {
                        int ahora=1;
                        int siguiente=2;

                        matriz * Prop1_piezas[]={EM2tan4,EM2tan3,EM2tan2,EM2tan1};
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        //q_2=prop_q(Prop1,obtieneElemento(q_iter,1,ahora),1,n0);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={EM2tan13,EM2tan12,EM2tan11,EM2tan10,EM2tan9,EM2tan8,EM2tan7,EM2tan6,EM2tan5}; // 5 a 13
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,9);
                        q_4=prop_q(Prop2,q_3,n0,n0);
                        //q_5=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_4,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_5=propagacionKerrRK(q_4,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop3_piezas[]={EM2tan17,EM2tan16,EM2tan15,EM2tan14}; // 14 a 17
                        matriz * Prop3=variasMultMatriciales(Prop3_piezas,4);
                        q_6=prop_q(Prop3,q_5,n0,1);
                        fijaElemento(q_iter,1,ahora,q_6);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieaza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        Prop3=borraMatriz(Prop3);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_6,1,lambda0);// Prepara vectores
                        if(cuenta>iteraciones)  break; // Por si acaso
                        ++cuenta;
                        //printf("%i\n",cuenta);
                        void *tmp_spt, *tmp_iter;
                        tmp_spt=(long double*)realloc(spots,cuenta*sizeof(long double));
                        tmp_iter=(int *)realloc(iteracionActual,cuenta*sizeof(int));
                        if(tmp_spt!=NULL && tmp_iter!=NULL)
                        {
                            spots=tmp_spt;
                            iteracionActual=tmp_iter;
                            spots[cuenta-1]=w_iter_nuevo;
                            iteracionActual[cuenta-1]=cuenta-1;
                        }
                        else
                        {
                            reallocFallido=1;
                            break;
                        }
                        if(w_iter_viejo>0.1||w_iter_nuevo>0.1||isnan(w_iter_viejo) ||isnan(w_iter_nuevo))
                        {
                            termino=0;
                            escribeCSV=1;
                            fijaElemento(q_new,index1,index2,0);
                            break;
                        }  // Si el spot mide 20 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo)*100.00;
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,ahora));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(reallocFallido==1)
                {
                    termino=0;
                    free(spots);
                    free(iteracionActual);
                    printf("Se acabó la memoria\n");
                    reallocFallido=0;
                }
                guardaSpotIteraciones(spots[cuenta-1],cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,termino);
                if(escribeCSV==1)
                    guardaVariaciones(spots,iteracionActual,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,subCarpeta,termino);
                if(termino==0)
                    printf("Cavidad en X, EM2tanRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en X, EM2tanRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                free(spots);
                free(iteracionActual);
            }
        }
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
        return q_new;
    }

    matriz * propNoLinealRK_sagEM2(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double L1, long double L2, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        _Bool escribeCSV=0;
        _Bool reallocFallido=0;
        long double *spots;
        int * iteracionActual;
        int cuenta=0;
        char tipo[]="EM2SagRK";
        char subCarpeta[]="EM2Sag_RK";
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4,q_5,q_6; // Parámetros complejos para propagación

        // Crea matrices para la propagación

        matriz *EM2sag1,*EM2sag2,*EM2sag3,*EM2sag4,*EM2sag5,*EM2sag6,*EM2sag7,*EM2sag8,*EM2sag9,*EM2sag10,*EM2sag11,*EM2sag12,*EM2sag13,*EM2sag14,*EM2sag15,*EM2sag16,*EM2sag17;
        // EM2, Caso saggencial
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
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                spots=(long double*)calloc(1,sizeof(long double));
                iteracionActual=(int*)calloc(1,sizeof(int));
                cuenta=1;
                iteracionActual[0]=0;
                spots[0]=spot_q(obtieneElemento(q_in,index1,index2),1,lambda0);
                // Llenando matrices variables:
                llenaMatrizSinCrear(EM2sag3,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(EM2sag6,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(EM2sag12,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(EM2sag15,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);


                if(elem(q_in,index1,index2)==0)
                    {elem(q_new,index1,index2)=0;
                    escribeCSV=0;
                    reallocFallido=0;
                    }
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    escribeCSV=1;                //for(int iter=1;iter<iteraciones;iter++)
                     for(;;)
                    {
                        int ahora=1;
                        int siguiente=2;

                        matriz * Prop1_piezas[]={EM2sag4,EM2sag3,EM2sag2,EM2sag1};
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        //q_2=prop_q(Prop1,obtieneElemento(q_iter,1,ahora),1,n0);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={EM2sag13,EM2sag12,EM2sag11,EM2sag10,EM2sag9,EM2sag8,EM2sag7,EM2sag6,EM2sag5}; // 5 a 13
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,9);
                        q_4=prop_q(Prop2,q_3,n0,n0);
                        //q_5=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_4,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_5=propagacionKerrRK(q_4,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop3_piezas[]={EM2sag17,EM2sag16,EM2sag15,EM2sag14}; // 14 a 17
                        matriz * Prop3=variasMultMatriciales(Prop3_piezas,4);
                        q_6=prop_q(Prop3,q_5,n0,1);
                                            fijaElemento(q_iter,1,ahora,q_6);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieaza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        Prop3=borraMatriz(Prop3);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_6,1,lambda0);// Prepara vectores
                        if(cuenta>iteraciones)  break; // Por si acaso
                        ++cuenta;
                        //printf("%i\n",cuenta);
                        void *tmp_spt, *tmp_iter;
                        tmp_spt=(long double*)realloc(spots,cuenta*sizeof(long double));
                        tmp_iter=(int *)realloc(iteracionActual,cuenta*sizeof(int));
                        if(tmp_spt!=NULL && tmp_iter!=NULL)
                        {
                            spots=tmp_spt;
                            iteracionActual=tmp_iter;
                            spots[cuenta-1]=w_iter_nuevo;
                            iteracionActual[cuenta-1]=cuenta-1;
                        }
                        else
                        {
                            reallocFallido=1;
                            break;
                        }
                        if(w_iter_viejo>0.1||w_iter_nuevo>0.1||isnan(w_iter_viejo) ||isnan(w_iter_nuevo))
                        {
                            termino=0;
                            escribeCSV=1;
                            fijaElemento(q_new,index1,index2,0);
                            break;
                        }  // Si el spot mide 20 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo);
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,ahora));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(reallocFallido==1)
                {
                    termino=0;
                    free(spots);
                    free(iteracionActual);
                    printf("Se acabó la memoria\n");
                    reallocFallido=0;
                }
                guardaSpotIteraciones(spots[cuenta-1],cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,termino);
                if(escribeCSV==1)
                    guardaVariaciones(spots,iteracionActual,cuenta,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),tipo,subCarpeta,termino);
                if(termino==0)
                    printf("Cavidad en X, EM2sagRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en X, EM2sagRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                free(spots);
                free(iteracionActual);
            }
        }
        borraMatriz(q_iter);
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

    matriz * propNoLinealRK_AnilloRuta1Tan(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        long double L1=a;
        long double L2=b+c;
        int pasosKerrRK=1000;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        assert(strcmp(conjugado_corto,conjugado_largo)==0);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4; // Parámetros complejos para propagación

        // Crea matrices para la propagación
        matriz *Ruta1tan1,*Ruta1tan2,*Ruta1tan3,*Ruta1tan4,*Ruta1tan5,*Ruta1tan6,*Ruta1tan7,*Ruta1tan8;
       // Ruta 1, caso tangencial
        Ruta1tan1=llenaMatriz(1,L1,0,1); // L1=a
        Ruta1tan2=llenaMatriz(1,0,-1/f1t,1);
        Ruta1tan3=nuevaMatriz(2,2); // Pendiente
        Ruta1tan4=llenaMatriz(n0,0,0,1/n0); // Medio Kerr (entrada)
        Ruta1tan5=llenaMatriz(1/n0,0,0,n0); // Medio Kerr (salida)
        Ruta1tan6=nuevaMatriz(2,2);//Ruta1tan6 pendiente
        Ruta1tan7=llenaMatriz(1,0,-1/f2t,1);
        Ruta1tan8=llenaMatriz(1,L2,0,1); // L2=b+c

        // Ciclo
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                // Llenando matrices variables:
                llenaMatrizSinCrear(Ruta1tan3,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(Ruta1tan6,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);;
                if(elem(q_in,index1,index2)==0)
                    {elem(q_new,index1,index2)=0;}
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    //for(int iter=1;iter<iteraciones;iter++)
                    for(int iter=iteraciones-1;iter>0;iter--)
                    {
                         //int ahora=iter;
                        //int siguiente=iter+1;
                        int ahora=iteraciones-iter;
                        int siguiente=ahora+1;
                        matriz * Prop1_piezas[]={Ruta1tan4,Ruta1tan3,Ruta1tan2,Ruta1tan1}; //FM1 a refracción en medio Kerr
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={Ruta1tan8,Ruta1tan7,Ruta1tan6,Ruta1tan5}; // 5 a 8
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,4);
                        q_4=prop_q(Prop2,q_3,n0,1.0);
                        fijaElemento(q_iter,1,siguiente,q_4);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_4,1,lambda0);
                        if(w_iter_viejo>0.2||w_iter_nuevo>0.2) break;  // Si el spot mide 40 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo)*100.00;
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,siguiente));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(termino==0)
                    printf("Cavidad en anillo, Ruta1TanRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en anillo, Ruta1TanRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
            }
        }
        q_iter=borraMatriz(q_iter);
        Ruta1tan1=borraMatriz(Ruta1tan1);
        Ruta1tan2=borraMatriz(Ruta1tan2);
        Ruta1tan3=borraMatriz(Ruta1tan3);
        Ruta1tan4=borraMatriz(Ruta1tan4);
        Ruta1tan5=borraMatriz(Ruta1tan5);
        Ruta1tan6=borraMatriz(Ruta1tan6);
        Ruta1tan7=borraMatriz(Ruta1tan7);
        Ruta1tan8=borraMatriz(Ruta1tan8);
        return q_new;
    }

    matriz * propNoLinealRK_AnilloRuta1Sag(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        long double L1=a;
        long double L2=b+c;
        int pasosKerrRK=1000;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        assert(strcmp(conjugado_corto,conjugado_largo)==0);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4; // Parámetros complejos para propagación

        // Crea matrices para la propagación
        matriz *Ruta1sag1,*Ruta1sag2,*Ruta1sag3,*Ruta1sag4,*Ruta1sag5,*Ruta1sag6,*Ruta1sag7,*Ruta1sag8;
       // Ruta 1, caso sagital
        Ruta1sag1=llenaMatriz(1,L1,0,1); // L1=a
        Ruta1sag2=llenaMatriz(1,0,-1/f1s,1);
        Ruta1sag3=nuevaMatriz(2,2); // Pendiente
        Ruta1sag4=llenaMatriz(1,0,0,1); // Medio Kerr (entrada)
        Ruta1sag5=llenaMatriz(1,0,0,1); // Medio Kerr (salida)
        Ruta1sag6=nuevaMatriz(2,2);//Ruta1sag6 pendiente
        Ruta1sag7=llenaMatriz(1,0,-1/f2s,1);
        Ruta1sag8=llenaMatriz(1,L2,0,1); // L2=b+c

        // Ciclo
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                // Llenando matrices variables:
                llenaMatrizSinCrear(Ruta1sag3,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);
                llenaMatrizSinCrear(Ruta1sag6,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);;
                if(elem(q_in,index1,index2)==0)
                    {elem(q_new,index1,index2)=0;}
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    //for(int iter=1;iter<iteraciones;iter++)
                    for(int iter=iteraciones-1;iter>0;iter--)
                    {
                         //int ahora=iter;
                        //int siguiente=iter+1;
                        int ahora=iteraciones-iter;
                        int siguiente=ahora+1;
                        matriz * Prop1_piezas[]={Ruta1sag4,Ruta1sag3,Ruta1sag2,Ruta1sag1}; //FM1 a refracción en medio Kerr
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={Ruta1sag8,Ruta1sag7,Ruta1sag6,Ruta1sag5}; // 5 a 8
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,4);
                        q_4=prop_q(Prop2,q_3,n0,1);
                        fijaElemento(q_iter,1,siguiente,q_4);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_4,1,lambda0);
                        if(w_iter_viejo>0.2||w_iter_nuevo>0.2) break;  // Si el spot mide 40 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo)*100.00;
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            //fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,iter+1));
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,siguiente));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(termino==0)
                    printf("Cavidad en anillo, Ruta1SagRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en anillo, Ruta1SagRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
            }
        }
        q_iter=borraMatriz(q_iter);
        Ruta1sag1=borraMatriz(Ruta1sag1);
        Ruta1sag2=borraMatriz(Ruta1sag2);
        Ruta1sag3=borraMatriz(Ruta1sag3);
        Ruta1sag4=borraMatriz(Ruta1sag4);
        Ruta1sag5=borraMatriz(Ruta1sag5);
        Ruta1sag6=borraMatriz(Ruta1sag6);
        Ruta1sag7=borraMatriz(Ruta1sag7);
        Ruta1sag8=borraMatriz(Ruta1sag8);
        return q_new;
    }

    matriz * propNoLinealRK_AnilloRuta2Tan(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        long double L1=a;
        long double L2=b+c;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        assert(strcmp(conjugado_corto,conjugado_largo)==0);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4; // Parámetros complejos para propagación

        // Crea matrices para la propagación
        matriz *Ruta2tan1,*Ruta2tan2,*Ruta2tan3,*Ruta2tan4,*Ruta2tan5,*Ruta2tan6,*Ruta2tan7,*Ruta2tan8;
       // Ruta 2, caso tangencial
        Ruta2tan1=llenaMatriz(1,L2,0,1); // L2=b+c
        Ruta2tan2=llenaMatriz(1,0,-1/f2t,1);
        Ruta2tan3=nuevaMatriz(2,2); // Pendiente
        Ruta2tan4=llenaMatriz(n0,0,0,1/n0); // Medio Kerr (entrada)
        Ruta2tan5=llenaMatriz(1/n0,0,0,n0); // Medio Kerr (salida)
        Ruta2tan6=nuevaMatriz(2,2);//Ruta2tan6 pendiente
        Ruta2tan7=llenaMatriz(1,0,-1/f1t,1);
        Ruta2tan8=llenaMatriz(1,L1,0,1); // L1=a

        // Ciclo
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                // Llenando matrices variables:
                llenaMatrizSinCrear(Ruta2tan3,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(Ruta2tan6,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);;
                if(elem(q_in,index1,index2)==0)
                    {elem(q_new,index1,index2)=0;}
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                    //for(int iter=1;iter<iteraciones;iter++)
                    for(int iter=iteraciones-1;iter>0;iter--)
                    {
                         //int ahora=iter;
                        //int siguiente=iter+1;
                        int ahora=iteraciones-iter;
                        int siguiente=ahora+1;
                        matriz * Prop1_piezas[]={Ruta2tan4,Ruta2tan3,Ruta2tan2,Ruta2tan1}; //FM1 a refracción en medio Kerr
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={Ruta2tan8,Ruta2tan7,Ruta2tan6,Ruta2tan5}; // 5 a 8
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,4);
                        q_4=prop_q(Prop2,q_3,n0,1);
                        fijaElemento(q_iter,1,siguiente,q_4);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_4,1,lambda0);
                        if(w_iter_viejo>0.2||w_iter_nuevo>0.2) break;  // Si el spot mide 40 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo);
                        if(diferencia<=umbral) // Termina la iteración
                        {
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,siguiente));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(termino==0)
                    printf("Cavidad en anillo, Ruta2TanRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en anillo, Ruta2TanRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
            }
        }
        q_iter=borraMatriz(q_iter);
        Ruta2tan1=borraMatriz(Ruta2tan1);
        Ruta2tan2=borraMatriz(Ruta2tan2);
        Ruta2tan3=borraMatriz(Ruta2tan3);
        Ruta2tan4=borraMatriz(Ruta2tan4);
        Ruta2tan5=borraMatriz(Ruta2tan5);
        Ruta2tan6=borraMatriz(Ruta2tan6);
        Ruta2tan7=borraMatriz(Ruta2tan7);
        Ruta2tan8=borraMatriz(Ruta2tan8);
        return q_new;
    }

    matriz * propNoLinealRK_AnilloRuta2Sag(matriz *q_in, char *conjugado_corto,char *conjugado_largo,long double a, long double b, long double c, long double L, long double f1, long double f2, long double n0, long double n2, long double w_pump,long double chi, long double kth, long double Cp, long double rho,long double dn_dv,long double P_pump,long double P_laser,long double lambda0,matriz *epsilon1, matriz *epsilon2, int iteraciones, long double umbral)
    {
        // Variables
        //imprimeMatriz(q_in);
        long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
        _Bool termino=0;
        long double w_iter_viejo=0, w_iter_nuevo=0, diferencia=0;
        long double L1=a;
        long double L2=b+c;
        assert(strlen(conjugado_corto)==3);
        assert(strlen(conjugado_largo)==3);
        assert(strcmp(conjugado_corto,conjugado_largo)==0);
        anguloLineal(conjugado_corto,conjugado_largo,L,n0,L1,L2,f1,f2,angulos);
        delta1=distanciaCristal(conjugado_corto,f1,L1,L,n0,angulos[0]);
        delta2=distanciaCristal(conjugado_largo,f2,L2,L,n0,angulos[1]);

        // Cálculo de distancias focales
        f1t=f1*cosl(angulos[0]);
        f2t=f2*cosl(angulos[1]);
        f1s=f1/cosl(angulos[0]);
        f2s=f2/cosl(angulos[1]);

        matriz * q_new=nuevaMatriz(q_in->filas,q_in->columnas); // Creando matriz de salida
        long double _Complex q_2,q_3,q_4; // Parámetros complejos para propagación

        // Crea matrices para la propagación
        matriz *Ruta2sag1,*Ruta2sag2,*Ruta2sag3,*Ruta2sag4,*Ruta2sag5,*Ruta2sag6,*Ruta2sag7,*Ruta2sag8;
       // Ruta 2, caso sagital
        Ruta2sag1=llenaMatriz(1,L2,0,1); // L2=b+c
        Ruta2sag2=llenaMatriz(1,0,-1/f2s,1);
        Ruta2sag3=nuevaMatriz(2,2); // Pendiente
        Ruta2sag4=llenaMatriz(1,0,0,1); // Medio Kerr (entrada)
        Ruta2sag5=llenaMatriz(1,0,0,1); // Medio Kerr (salida)
        Ruta2sag6=nuevaMatriz(2,2);//Ruta2sag6 pendiente
        Ruta2sag7=llenaMatriz(1,0,-1/f1s,1);
        Ruta2sag8=llenaMatriz(1,L1,0,1); // L1=a

        // Ciclo
        matriz * q_iter = nuevaMatriz(1,iteraciones);
        for(int index1=1;index1<=epsilon1->columnas;index1++)
        {
            for(int index2=1;index2<=epsilon2->columnas;index2++)
            {
                // Llenando matrices variables:
                llenaMatrizSinCrear(Ruta2sag3,1,delta2+creall(obtieneElemento(epsilon2,1,index2)),0,1);
                llenaMatrizSinCrear(Ruta2sag6,1,delta1+creall(obtieneElemento(epsilon1,1,index1)),0,1);;
                if(elem(q_in,index1,index2)==0)
                    {elem(q_new,index1,index2)=0;}
                else
                {
                    fijaElemento(q_iter,1,1,obtieneElemento(q_in,index1,index2));
                //for(int iter=1;iter<iteraciones;iter++)
                    for(int iter=iteraciones-1;iter>0;iter--)
                    {
                         //int ahora=iter;
                        //int siguiente=iter+1;
                        int ahora=iteraciones-iter;
                        int siguiente=ahora+1;
                        matriz * Prop1_piezas[]={Ruta2sag4,Ruta2sag3,Ruta2sag2,Ruta2sag1}; //FM1 a refracción en medio Kerr
                        matriz * Prop1=variasMultMatriciales(Prop1_piezas,4);
                        long double _Complex q1_hold=obtieneElemento(q_iter,1,ahora);
                        q_2=prop_q(Prop1,q1_hold,1,n0);
                        //q_3=propagacionKerr(pasosKerrRK,L,n0,n2,w_pump,q_2,chi,kth,Cp,rho,dn_dv,P_pump,P_laser,lambda0);
                        q_3=propagacionKerrRK(q_2,lambda0,L,n0,n2,P_laser,pasosKerrRK);
                        matriz * Prop2_piezas[]={Ruta2sag8,Ruta2sag7,Ruta2sag6,Ruta2sag5}; // 5 a 8
                        matriz * Prop2=variasMultMatriciales(Prop2_piezas,4);
                        q_4=prop_q(Prop2,q_3,n0,1);
                        fijaElemento(q_iter,1,siguiente,q_4);
                        //printf("Actual: %f+i %f, siguiente: %f+i%f",creal(obtieneElemento(q_iter,1,iter)),cimag(obtieneElemento(q_iter,1,iter)),creal(obtieneElemento(q_iter,1,iter+1)),cimag(obtieneElemento(q_iter,1,iter+1)));
                        // Limpieza
                        Prop1=borraMatriz(Prop1);
                        Prop2=borraMatriz(Prop2);
                        // Calculando spots y comparación
                        w_iter_viejo=spot_q(q1_hold,1,lambda0);
                        w_iter_nuevo=spot_q(q_4,1,lambda0);
                        if(w_iter_viejo>0.2||w_iter_nuevo>0.2) break;  // Si el spot mide 40 cm es demasiado
                        diferencia=fabsl((w_iter_nuevo-w_iter_viejo)/w_iter_nuevo)*100.00;
                        if(diferencia<=umbral) // Termina la iteración
                       {
                            fijaElemento(q_new,index1,index2,obtieneElemento(q_iter,1,siguiente));
                            termino=1;
                            break;
                        }
                        else
                        {
                            fijaElemento(q_new,index1,index2,0);
                            termino=0;
                        }
                    }
                }
                if(termino==0)
                    printf("Cavidad en anillo, Ruta2SagRK: No cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
                else
                    printf("Cavidad en anillo, Ruta2SagRK: Cumple con el umbral de %Le para epsilon1=%Le, epsilon2=%Le.\nError relativo: %Le\n",umbral,creall(obtieneElemento(epsilon1,1,index1)),creall(obtieneElemento(epsilon2,1,index2)),diferencia);
            }
        }
        q_iter=borraMatriz(q_iter);
        Ruta2sag1=borraMatriz(Ruta2sag1);
        Ruta2sag2=borraMatriz(Ruta2sag2);
        Ruta2sag3=borraMatriz(Ruta2sag3);
        Ruta2sag4=borraMatriz(Ruta2sag4);
        Ruta2sag5=borraMatriz(Ruta2sag5);
        Ruta2sag6=borraMatriz(Ruta2sag6);
        Ruta2sag7=borraMatriz(Ruta2sag7);
        Ruta2sag8=borraMatriz(Ruta2sag8);
        return q_new;
    }

