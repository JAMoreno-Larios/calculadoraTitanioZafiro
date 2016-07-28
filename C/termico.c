#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include "matrices.h"
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

long double * gaussSeidel(long double **A, long double *B, int n, int iteraciones, long double umbral, long double lambda)
{
    long double *X=(long double *)calloc(n,sizeof(long double));
    for(int i=0;i<n;i++)
    {
        long double dummy=A[i][i];
        for(int j=0;j<n;j++)
        {
            A[i][j]/=dummy; //A[i][j]=A[i][j]/A[i][i];
        }
        B[i]/=dummy; //B[i]=B[i]/A[i][i];
    }
    for(int i=0;i<n;i++)
    {
        long double suma=B[i];
        for(int j=0;j<n;j++)
        {
            if(i!=j)
                suma-=A[i][j]*X[j];
        }
        X[i]=suma;
    }
    for(int paso=0;paso<iteraciones;paso++)
    {
        long double error=0;
        _Bool centinela=1;
        for(int i=0;i<n;i++)
        {
            long double viejo=X[i];
            long double suma=B[i];
            for(int j=0;j<n;j++)
            {
                if(i!=j) suma-=A[i][j]*X[j];
            }
            X[i]=lambda*suma+(1.0-lambda)*viejo;
            if(centinela==1&&X[i]!=0)
            {
                error=100.0*fabsl((X[i]-viejo)/X[i]);
                if(error>umbral) centinela=0;
            }
        }
        if(centinela==1) break;
    }
    return X;
}

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


distTempPlano * planoNuevo(int N, int paso, long double ancho, long double alto) // Crea estructura y almacena datos
{
    distTempPlano *dT =(distTempPlano*)calloc(1,sizeof(distTempPlano));
    dT->N=N;
    dT->paso=paso;
    dT->deltaX=ancho/(N-1);
    dT->deltaY=alto/(N-1);
    dT->X=(long double*)calloc(N,sizeof(long double));
    dT->Y=(long double*)calloc(N,sizeof(long double));
    dT->X[N-1]=ancho/2; // Se hace esa resta por cuestión de índices
    dT->Y[N-1]=alto/2;
    for(int i=N-1;i>0;i--) // Creación de vectores de posición
    {
        dT->X[i-1]=dT->X[i]-dT->deltaX;
        dT->Y[i-1]=dT->Y[i]-dT->deltaY;
    }
    // Plano
    dT->T=(long double *)calloc(N*N,sizeof(long double)); //Reserva para mapa de temperatura
    return dT;
}

distTempPlano * borraPlano(distTempPlano * plano) // Borra estructura de plano
{
    assert(plano);
    // Libera datos
    assert(plano->T);
    assert(plano->X);
    assert(plano->Y);
    free(plano->X);
    free(plano->Y);
    free(plano->T);
    plano->T=NULL;
    plano->X=NULL;
    plano->Y=NULL;
    // Libera a plano
    free(plano);
    plano=NULL;
    return plano;
}

distTempPlano * distTempPlanoCristal(int N, long double ancho, long double alto, long double deltaZ, long double P_pump, long double chi, long double L, long double kth, long double Cp, long double rho, long double w_pump_t,long double w_pump_s, long double alpha,\
                                      int paso, int iteraciones, long double umbral, long double lambda_relax)
{
    assert((N%2==0)!=1);
    distTempPlano *dT=planoNuevo(N,paso,ancho,alto);
    long double xi=(long double)dT->paso*deltaZ;
        // Llena condiciones de frontera
        dT->T[0+0]=dT->T[N-1]=dT->T[(N-1)*N]=dT->T[(N-1)*N+N-1]=20.0+273.15; // Valores no usados en el algoritmo
        for(int i=1;i<N-1;i++)
        {
                dT->T[0*N+i]=20.0+273.15; // Frontera izquierda
                dT->T[i*N+0]=20.0+273.15; // Frontera inferior
                dT->T[(N-1)*N+i]=20.0+273.15; // Frontera derecha
                dT->T[i*N+N-1]=20.0+273.15; // Frontera superior
        }
        // Iterando
        for(int paso=0;paso<iteraciones;paso++)
        {
            _Bool centinela=1;
            for(int j=1;j<N-1;j++)
            {
                for(int i=1;i<N-1;i++)
                {
                    long double func=(-2.0*chi*P_pump)/(L*Cp*rho*kth*M_PI*w_pump_t*w_pump_s)*expl(-2.0*powl(dT->X[i],2)/powl(w_pump_t,2)-2.0*powl(dT->Y[j],2)/powl(w_pump_s,2))*expl(-alpha*xi);
                    long double viejo=dT->T[i*N+j];
                    dT->T[i*N+j]=(powl(dT->deltaY,2)*(dT->T[(i+1)*N+j]+dT->T[(i-1)*N+j])+powl(dT->deltaX,2)*(dT->T[i*N+(j+1)]+dT->T[i*N+(j-1)])-powl(dT->deltaX,2)*powl(dT->deltaY,2)*func)/(2*(powl(dT->deltaX,2)+powl(dT->deltaY,2)));
                    dT->T[i*N+j]=lambda_relax*dT->T[i*N+j]+(1.0-lambda_relax)*viejo;
                    if(centinela==1&&dT->T[i*N+j]!=0)
                    {
                        long double error=fabsl((dT->T[i*N+j]-viejo)/dT->T[i*N+j])*100.0;
                        if(error>umbral) centinela=0;
                    }
                }
            }
            if(centinela==1)
                {break;}
        }

        // Obteniento deltaT
        int centro=(N+1)/2;
        long double Tmax=dT->T[(centro-1)*N+centro-1]; // Centro ubicado en N-2, N-2. Se resta 1 por cuestiones de índices.
        for(int j=N-1;j>=0;j--)
        {
            for(int i=N-1;i>=0;i--)
            {
                dT->T[i*N+j]=dT->T[i*N+j]-Tmax;
            }
        }
    //}
    // Imprimiendo valores
    /*for(int j=1;j<N-1;j++)
    {
        for(int i=1;i<N-1;i++)
        {
            printf("Delta T(%i,%i)= %Lf [K]\n",i,j,dT->T[i*N+j]);
        }
    }*/
    return dT;
}

ajusteTemperaturaPlano * ajusteCuadraticoPlano(distTempPlano * plano, int iteraciones, long double umbral, long double lambda_relax)
{
    // Creación de estructura
    ajusteTemperaturaPlano *ajuste =(ajusteTemperaturaPlano*)calloc(1,sizeof(ajusteTemperaturaPlano));
    // Cálculo de coeficientes
    int N=plano->N;
    long double AA,BB,CC,DD,EE,FF,GG,HH,II; //
    AA=0.0, BB=0.0, CC=0.0, DD=0.0, EE=0.0, FF=0.0, GG=0.0, HH=0.0, II=0.0;
    int i=0,j=0,conteo=0;
    for(j=1;j<N-1;j++)
    {
        for(i=1;i<N-1;i++)
        {
            AA++;
            BB+=powl(plano->X[i],2);
            CC+=powl(plano->Y[j],2);
            DD+=plano->T[i*N+j];
            EE+=powl(plano->X[i],4);
            FF+=powl(plano->X[i],2)*powl(plano->Y[j],2);
            GG+=plano->T[i*N+j]*powl(plano->X[i],2);
            HH+=powl(plano->Y[j],4);
            II+=plano->T[i*N+j]*powl(plano->Y[j],2);
            conteo++;
        }
    }

    // Construyendo matriz

    //int iteraciones=100000;
    //long double umbral=0.1e-6; // Umbral en porcentaje
    // Creación de matriz A
    int orden=3; // Matriz de ajuste
    long double **A=(long double**)calloc(orden,sizeof(long double *));
    for(i=0;i<orden;i++)
        A[i]=(long double *)calloc(orden,sizeof(long double));
    // Llenado de A
    A[0][0]=AA;
    A[0][1]=BB;
    A[0][2]=CC;
    A[1][0]=BB;
    A[1][1]=EE;
    A[1][2]=FF;
    A[2][0]=CC;
    A[2][1]=FF;
    A[2][2]=HH;

    // Creación de B
    long double *B=(long double *)calloc(orden,sizeof(long double));

    // Llenado de B
    B[0]=DD;
    B[1]=GG;
    B[2]=II;

    long double *X=gaussSeidel(A,B,orden,iteraciones,umbral,lambda_relax);
    /*printf("Soluciones\n");
    printf("a0=%Lf\n",X[0]);
    printf("a1=%Lf\n",X[1]);
    printf("a2=%Lf\n",X[2]);*/

    /*FILE *archivo;
    archivo = fopen("Ajuste cuadratico.csv","w");
    fprintf(archivo, "Coeficiente,Valor,\n");
    for(i=0;i<3;i++)
    {
        fprintf(archivo, "a %i,%.15Le,\n",i,X[i]);
    }
    fclose(archivo);*/
    ajuste->a0=X[0];
    ajuste->a1=X[1];
    ajuste->a2=X[2];
    ajuste->paso=plano->paso;

    //Limpieza
    for(int i=orden-1;i>=0;--i)
        free(A[i]);
    free(A);
    free(B);
    free(X);
    return ajuste;
}

// Ajuste cuadrático ponderado a T original
ajusteTemperaturaPlano * ajusteCuadraticoPlanoPonderado(distTempPlano * plano, int iteraciones, long double umbral, long double lambda_relax)
{
    // Creación de estructura
    ajusteTemperaturaPlano *ajustePonderado =(ajusteTemperaturaPlano*)calloc(1,sizeof(ajusteTemperaturaPlano));
    // Cálculo de coeficientes
    int N=plano->N;
    long double AA,BB,CC,DD,EE,FF,GG,HH,II; //
    AA=0.0, BB=0.0, CC=0.0, DD=0.0, EE=0.0, FF=0.0, GG=0.0, HH=0.0, II=0.0;
    int i=0,j=0,conteo=0;
    for(j=1;j<N-1;j++)
    {
        for(i=1;i<N-1;i++)
        {
            AA+=plano->T[i*N+j];
            BB+=powl(plano->X[i],2)*plano->T[i*N+j];
            CC+=powl(plano->Y[j],2)*plano->T[i*N+j];
            DD+=plano->T[i*N+j]*plano->T[i*N+j];
            EE+=powl(plano->X[i],4)*plano->T[i*N+j];
            FF+=powl(plano->X[i],2)*powl(plano->Y[j],2)*plano->T[i*N+j];
            GG+=plano->T[i*N+j]*powl(plano->X[i],2)*plano->T[i*N+j];
            HH+=powl(plano->Y[j],4)*plano->T[i*N+j];
            II+=plano->T[i*N+j]*powl(plano->Y[j],2)*plano->T[i*N+j];
            conteo++;
        }
    }

    // Construyendo matriz

    //int iteraciones=100000;
    //long double umbral=0.1e-6; // Umbral en porcentaje
    // Creación de matriz A
    int orden=3; // Matriz de ajuste
    long double **A=(long double**)calloc(orden,sizeof(long double *));
    for(i=0;i<orden;i++)
        A[i]=(long double *)calloc(orden,sizeof(long double));
    // Llenado de A
    A[0][0]=AA;
    A[0][1]=BB;
    A[0][2]=CC;
    A[1][0]=BB;
    A[1][1]=EE;
    A[1][2]=FF;
    A[2][0]=CC;
    A[2][1]=FF;
    A[2][2]=HH;

    // Creación de B
    long double *B=(long double *)calloc(orden,sizeof(long double));

    // Llenado de B
    B[0]=DD;
    B[1]=GG;
    B[2]=II;

    long double *X=gaussSeidel(A,B,orden,iteraciones,umbral,lambda_relax);

    ajustePonderado->a0=X[0];
    ajustePonderado->a1=X[1];
    ajustePonderado->a2=X[2];
    ajustePonderado->paso=plano->paso;

    //Limpieza
    for(int i=orden-1;i>=0;--i)
        free(A[i]);
    free(A);
    free(B);
    free(X);
    return ajustePonderado;
}

ajusteTemperaturaCristal * vectorAjusteCristal(long double complex qInTan, long double complex qInSag, long double alpha, long double lambdaIn, long double nLin, \
    long double Power, long double chi, long double L, long double kth, long double Cp, long double rho, long double dn_dT, \
    long double ancho, long double alto, int iteraciones, int pasos, long double umbral, long double lambda_relax, int N, _Bool termico)
    {
        long double deltaZ=L/(pasos-1); // Diferencia de distancia en la dirección de propagación
        long double wTan=spot_q(qInTan,nLin,lambdaIn);
        long double wSag=spot_q(qInSag,nLin,lambdaIn);

        //Creación de estructura de salida
        ajusteTemperaturaCristal *ajuste=(ajusteTemperaturaCristal*)calloc(1,sizeof(ajusteTemperaturaCristal));
        ajuste->a0=(long double*)calloc(pasos,sizeof(long double));
        ajuste->a1=(long double*)calloc(pasos,sizeof(long double));
        ajuste->a2=(long double*)calloc(pasos,sizeof(long double));
        ajuste->paso=(int*)calloc(pasos,sizeof(int));

        // Llenado de paso[]
        for(int indice=pasos-1;indice>=0;indice--)
            ajuste->paso[indice]=indice;

        long double complex *qPropTan = (long double complex*)calloc(pasos,sizeof(long double complex));
        long double complex *qPropSag = (long double complex*)calloc(pasos,sizeof(long double complex));
        long double complex At, Bt, Ct, Dt;
        long double complex As, Bs, Cs, Ds;
        assert(qPropTan);
        assert(qPropSag);
        matriz *MTan, *MSag;
        MTan=nuevaMatriz(2,2);
        MSag=nuevaMatriz(2,2);
        // Cálculo inicial
        qPropTan[0]=qInTan;
        qPropSag[0]=qInSag;
        ajuste->paso[0]=0;
        long double n0Prop=nLin;
        if(termico==1)
        {
            for (int i=0;i<pasos;i++)
            {
                // Cálculo de spots
                long double w_pump_t=spot_q(qPropTan[i],n0Prop,lambdaIn);
                long double w_pump_s=spot_q(qPropSag[i],n0Prop,lambdaIn);
                // Cálculo de plano
                distTempPlano *plano=distTempPlanoCristal(N,ancho,alto,deltaZ,Power,chi,L,kth,Cp,rho,w_pump_s,w_pump_t,alpha,ajuste->paso[i],iteraciones,umbral,lambda_relax);
                // Ajuste de plano
                ajusteTemperaturaPlano *ajustePlano=ajusteCuadraticoPlanoPonderado(plano,iteraciones,umbral,lambda_relax);
                ajuste->a0[i]=ajustePlano->a0;
                ajuste->a1[i]=ajustePlano->a1;
                ajuste->a2[i]=ajustePlano->a2;
                ajuste->paso[i]=ajustePlano->paso;
                n0Prop=nLin+ajuste->a0[i]*dn_dT;
                // Formación de matrices
                // Factor de parábola=1/h^2
                long double parabola1=2*n0Prop*ajuste->a1[i]*dn_dT;
                long double parabola2=2*n0Prop*ajuste->a2[i]*dn_dT;
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
                free(ajustePlano);
                plano=borraPlano(plano);
            }
        }
        else
        {
            for (int i=0;i<pasos;i++)
            {
                // Cálculo de spots

                ajuste->a0[i]=0;
                ajuste->a1[i]=0;
                ajuste->a2[i]=0;
                ajuste->paso[i]=i;
                n0Prop=nLin+ajuste->a0[i]*dn_dT;
                // Formación de matrices
                // Factor de parábola=1/h^2
                long double parabola1=2*n0Prop*ajuste->a1[i]*dn_dT;
                long double parabola2=2*n0Prop*ajuste->a2[i]*dn_dT;
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
        //Limpieza
        borraMatriz(MTan);
        borraMatriz(MSag);
        free(qPropTan);
        free(qPropSag);
        //printf("Terminó ajuste térmico\n");
        return ajuste;
    }


void escribePlanoArchivoCompleto(distTempPlano * plano, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int i,j;
    fprintf(archivo, " ,Eje tangencial [m],");
    for(i=0;i<=plano->N-1;i++)
        fprintf(archivo, "%.15Le,",plano->X[i]);
    fprintf(archivo, "\nEje sagital [m],\n");
    for (j=0;j<=plano->N-1;j++)
    {
        fprintf(archivo,"%.15Le, ,",plano->Y[j]);
        for(i=0;i<=plano->N-1;i++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le,", plano->T[i*plano->N+j]);
        }
        fprintf(archivo,"\n");
    }
    fclose(archivo);
}

void escribePlanoArchivo(distTempPlano * plano, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int i,j;
    for (j=1;j<plano->N-1;j++)
    {
        for(i=1;i<plano->N-1;i++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le,", plano->T[i*plano->N+j]);
        }
        fprintf(archivo,"\n");
    }
    fclose(archivo);
}

void escribePlanoArchivoOrigin(distTempPlano * plano, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int i,j;
    for (j=1;j<plano->N-1;j++)
    {
        for(i=1;i<plano->N-1;i++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le,%.15Le,%.15Le\n", plano->X[i], plano->Y[j],plano->T[i*plano->N+j]);
        }
    }
    fclose(archivo);
}

void escribeCoeficientesAjuste(ajusteTemperaturaCristal *ajuste, int pasos, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    fprintf(archivo, "Paso,a0,a1,a2\n");
    int i;
    for(i=0;i<pasos;i++)
    {
        //Imprime el elemento de punto flotante
        fprintf(archivo, "%i,%.15Le,%.15Le,%.15Le\n", ajuste->paso[i],ajuste->a0[i],ajuste->a2[i],ajuste->a2[i]);
    }
    fclose(archivo);
}

ajusteTemperaturaCristal * borraVectorAjuste(ajusteTemperaturaCristal * vector)
{
    assert(vector);
    assert(vector->a0);
    free(vector->a0);
    vector->a0=NULL;
    assert(vector->a1);
    free(vector->a1);
    vector->a1=NULL;
    assert(vector->a2);
    free(vector->a2);
    vector->a2=NULL;
    assert(vector->paso);
    free(vector->paso);
    vector->paso=NULL;
    free(vector);
    vector=NULL;
    return vector;

}



// Rutina que genera vector térmico para cada variación en el spot de entrada.

ajusteTemperaturaCristal ** matrizAjusteCristal(spotsPump * spots, matriz * epsilon1, long double alpha, long double lambdaIn, long double nLin, \
    long double Power, long double chi, long double L, long double kth, long double Cp, long double rho, long double dn_dT, \
    long double ancho, long double alto, int iteraciones, int pasos, long double umbral, long double lambda_relax, int N, _Bool termico,
                                                CtoCppQProgressBar *bar)
    {
        CtoCppQProgressBar_setMaximum(bar,spots->numSpots);
        CtoCppQProgressBar_setMinimum(bar,0);
        ajusteTemperaturaCristal **matrizVectores=(ajusteTemperaturaCristal**)calloc(spots->numSpots,sizeof(ajusteTemperaturaCristal*));

        for(int i=0;i<spots->numSpots;i++)
        {
            matrizVectores[i]=vectorAjusteCristal(spots->qOutTan[i],spots->qOutSag[i], alpha, lambdaIn, nLin, \
                            Power, chi, L, kth, Cp, rho, dn_dT, ancho, alto, iteraciones, pasos, umbral, lambda_relax, N, termico);
            CtoCppQProgressBar_setValue(bar,i+1);
            CtoCppQProgressBar_repaint(bar);
        }
        return matrizVectores;
    }

ajusteTemperaturaCristal ** borraMatrizAjusteCristal(ajusteTemperaturaCristal ** matriz, int numVectores)
{
    for(int i=0;i<numVectores;i++)
    {
        matriz[i]=borraVectorAjuste(matriz[i]);
    }
    free(matriz);
    return NULL;
}

/*
 * ajusteTemperaturaCristal ** matrizAjusteCristal(spotsPump * spots, matriz * epsilon1, long double alpha, long double lambdaIn, long double nLin, \
    long double Power, long double chi, long double L, long double kth, long double Cp, long double rho, long double dn_dT, \
    long double ancho, long double alto, int iteraciones, int pasos, long double umbral, long double lambda_relax, int N, _Bool termico)
    {
        ajusteTemperaturaCristal **matrizVectores=(ajusteTemperaturaCristal**)calloc(spots->numSpots,sizeof(ajusteTemperaturaCristal*));

        for(int i=0;i<spots->numSpots;i++)
        {
            matrizVectores[i]=vectorAjusteCristal(spots->qOutTan[i],spots->qOutSag[i], alpha, lambdaIn, nLin, \
                            Power, chi, L, kth, Cp, rho, dn_dT, ancho, alto, iteraciones, pasos, umbral, lambda_relax, N, termico);
        }
        return matrizVectores;
    }
    */
