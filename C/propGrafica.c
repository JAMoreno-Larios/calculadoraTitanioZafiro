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

void escribeData(vectorDato *datos, char* filenameCW, char* filenameML)
{
    FILE *archivoCW, *archivoML;
    archivoCW=fopen(filenameCW,"w");
    archivoML=fopen(filenameML,"w");
    int i=0;
    for(i=0;i<datos->numeroElementos-1;i++)
    {
        fprintf(archivoCW,"%Le %Le\n",datos->pos[i],datos->spotCW[i]);
        fprintf(archivoML,"%Le %Le\n",datos->pos[i],datos->spotML[i]);
    }
    fclose(archivoCW);
    fclose(archivoML);
}

/* Escribe a gnuplot - Todos los datos. */
void gnuplotEscribe(vectorDato *EM1tan, vectorDato *EM1sag, vectorDato *EM2tan, vectorDato *EM2sag)
{
    // Nombres de archivo
    char* scriptFile="script.plot";
    char* datosEM1tanCW="datosEM1tanCW.dat";
    char* datosEM1sagCW="datosEM1sagCW.dat";
    char* datosEM2tanCW="datosEM2tanCW.dat";
    char* datosEM2sagCW="datosEM2sagCW.dat";
    char* datosEM1tanML="datosEM1tanML.dat";
    char* datosEM1sagML="datosEM1sagML.dat";
    char* datosEM2tanML="datosEM2tanML.dat";
    char* datosEM2sagML="datosEM2sagML.dat";

    // Creación de archivos de datos.
    escribeData(EM1tan, datosEM1tanCW, datosEM1tanML);
    escribeData(EM1sag, datosEM1sagCW, datosEM1sagML);
    escribeData(EM2tan, datosEM2tanCW, datosEM2tanML);
    escribeData(EM2sag, datosEM2sagCW, datosEM2sagML);

    FILE *script, *gnuplotPipe;
    script=fopen(scriptFile,"w");
    //Encabezado
    fprintf(script,"set encoding utf8\n");
    fprintf(script,"set terminal pdf enhanced solid size 6in,4in font \"Bookman URW,16\"\n");
    fprintf(script,"set datafile separator \" \"\n");
    // Conf. EM1 tan
    fprintf(script,"set title \"Propagación tangencial desde EM_1\\n{/*0.8 ε_1=-1 [mm], ε_2=-0.9 [mm]}\" font \"Bookman URW,16\" offset 0,0\n");
    fprintf(script,"set output \"EM1tan.pdf\"\n");
    fprintf(script,"set multiplot\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set size 1,1\n");
    fprintf(script,"set bmargin at screen 0.2\n");
    fprintf(script,"set rmargin at screen 0.93\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set key font \"Bookman URW,16\"\n");
    fprintf(script,"set xlabel \"Posición en la cavidad [m]\"\n");
    fprintf(script,"set ylabel \"Radio de haz ω_t [μm]\"\n");
    fprintf(script,"set xrange [*:*]\n");
    fprintf(script,"set xtics autofreq font \"Bookman URW,16\"\n");
    fprintf(script,"set yrange [*:*]\n");
    fprintf(script,"set ytics 2e-4 font \"Bookman URW,16\"\n");
    //fprintf(script,"set format x \"%%.2g\"\n");
    //fprintf(script,"set format y \"%%.1e\"\n");
    fprintf(script, "set label \'A: Interfaz Ti:Zafiro\' at 0.2,4e-4\n");
    fprintf(script, "set label \'B: Sección media\' at 0.2,3e-4\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM1tanCW,datosEM1tanML);

    // Zoom
    fprintf(script,"set size 0.45,0.2\n");
    fprintf(script,"set bmargin at screen 0.3\n");
    fprintf(script,"set rmargin at screen 0.88\n");
    fprintf(script,"set origin 0.5,0.4\n");
    fprintf(script,"set title \"Ti:Zafiro\" offset 0,-1\n");
    fprintf(script,"set xrange [%f:%f]\n",1.045,1.065);
    fprintf(script,"set yrange [%f:%f]\n",0.0,100e-6*1e6);
    fprintf(script,"set xtics (%f,%f) font \"Bookman URW,16\"\n",1.045,1.065);
    fprintf(script,"set ytics (%f,%f) font \"Bookman URW,16\"\n",0.0,100e-6*1e6);
    fprintf(script,"set xlabel \"\"\n");
    fprintf(script,"set ylabel \"\"\n");
    //fprintf(script,"set format x \"%%.4g\"\n");
    fprintf(script,"unset key\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM1tanCW,datosEM1tanML);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM1tan->posKerrIn,EM1tan->posKerrIn);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM1tan->posKerrMid,EM1tan->posKerrMid);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM1tan->posKerrOut,EM1tan->posKerrOut);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM1tan->posKerrIn);
    fprintf(script, "set label \'B\' at %Le,GPVAL_Y_MAX*0.9\n",EM1tan->posKerrMid);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM1tan->posKerrOut);
    fprintf(script,"replot\n");
    fprintf(script,"unset multiplot\n");

    // Conf. em1 sagital
    fprintf(script,"unset arrow\n");
    fprintf(script,"unset label\n");
    fprintf(script,"set title \"Propagación sagital desde EM_1\\n{/*0.8 ε_1=-1 [mm], ε_2=-0.9 [mm]}\" font \"Bookman URW,16\" offset 0,0\n");
    fprintf(script,"set output \"EM1sag.pdf\"\n");
    fprintf(script,"set multiplot\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set size 1,1\n");
    fprintf(script,"set bmargin at screen 0.2\n");
    fprintf(script,"set rmargin at screen 0.93\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set key font \"Bookman URW,16\"\n");
    fprintf(script,"set xlabel \"Posición en la cavidad [m]\"\n");
    fprintf(script,"set ylabel \"Radio de haz ω_s [μm]\"\n");
    fprintf(script,"set xrange [*:*]\n");
    fprintf(script,"set xtics autofreq font \"Bookman URW,16\"\n");
    fprintf(script,"set yrange [*:*]\n");
    fprintf(script,"set ytics 2e-4 font \"Bookman URW,16\"\n");
    //fprintf(script,"set format x \"%%.2g\"\n");
    //fprintf(script,"set format y \"%%.1e\"\n");
    fprintf(script, "set label \'A: Interfaz Ti:Zafiro\' at 0.2,4e-4\n");
    fprintf(script, "set label \'B: Sección media\' at 0.2,3e-4\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM1sagCW,datosEM1sagML);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM1sag->posKerrIn,EM1sag->posKerrIn);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM1sag->posKerrMid,EM1sag->posKerrMid);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM1sag->posKerrOut,EM1sag->posKerrOut);
    // Zoom
    fprintf(script,"set size 0.45,0.2\n");
    fprintf(script,"set bmargin at screen 0.3\n");
    fprintf(script,"set rmargin at screen 0.88\n");
    fprintf(script,"set origin 0.5,0.4\n");
    fprintf(script,"set title \"Ti:Zafiro\" offset 0,-1\n");
    fprintf(script,"set xrange [%f:%f]\n",1.045,1.065);
    fprintf(script,"set yrange [%f:%f]\n",0.0,100e-6*1e6);
    fprintf(script,"set xtics (%f,%f) font \"Bookman URW,16\"\n",1.045,1.065);
    fprintf(script,"set ytics (%f,%f) font \"Bookman URW,16\"\n",0.0,100e-6*1e6);
    fprintf(script,"set xlabel \"\"\n");
    fprintf(script,"set ylabel \"\"\n");
    //fprintf(script,"set format x \"%%.4g\"\n");

    fprintf(script,"unset key\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM1sagCW,datosEM1sagML);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM1sag->posKerrIn,EM1sag->posKerrIn);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM1sag->posKerrMid,EM1sag->posKerrMid);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM1sag->posKerrOut,EM1sag->posKerrOut);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM1sag->posKerrIn);
    fprintf(script, "set label \'B\' at %Le,GPVAL_Y_MAX*0.9\n",EM1sag->posKerrMid);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM1sag->posKerrOut);
    fprintf(script,"replot\n");
    fprintf(script,"unset multiplot\n");

    // Conf. EM2 tan
    fprintf(script,"unset arrow\n");
    fprintf(script,"unset label\n");
    fprintf(script,"set title \"Propagación tangencial desde EM_2\\n{/*0.8 ε_1=-1 [mm], ε_2=-0.9 [mm]}\" font \"Bookman URW,16\" offset 0,0\n");
    fprintf(script,"set output \"EM2tan.pdf\"\n");
    fprintf(script,"set multiplot\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set size 1,1\n");
    fprintf(script,"set bmargin at screen 0.2\n");
    fprintf(script,"set rmargin at screen 0.93\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set key font \"Bookman URW,16\"\n");
    fprintf(script,"set xlabel \"Posición en la cavidad [m]\"\n");
    fprintf(script,"set ylabel \"Radio de haz ω_t [μm]\"\n");
    fprintf(script,"set xrange [*:*]\n");
    fprintf(script,"set xtics autofreq font \"Bookman URW,16\"\n");
    fprintf(script,"set yrange [*:*]\n");
    fprintf(script,"set ytics 2e-4 font \"Bookman URW,16\"\n");
    //fprintf(script,"set format x \"%%.2g\"\n");
    //fprintf(script,"set format y \"%%.1e\"\n");
    fprintf(script, "set label \'A: Interfaz Ti:Zafiro\' at 0.2,4e-4\n");
    fprintf(script, "set label \'B: Sección media\' at 0.2,3e-4\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM2tanCW,datosEM2tanML);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2tan->posKerrIn,EM2tan->posKerrIn);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM2tan->posKerrMid,EM2tan->posKerrMid);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2tan->posKerrOut,EM2tan->posKerrOut);
    // Zoom
    fprintf(script,"set size 0.45,0.2\n");
    fprintf(script,"set bmargin at screen 0.3\n");
    fprintf(script,"set rmargin at screen 0.88\n");
    fprintf(script,"set origin 0.5,0.4\n");
    fprintf(script,"set title \"Ti:Zafiro\" offset 0,-1\n");
    fprintf(script,"set xrange [%f:%f]\n",1.045,1.065);
    fprintf(script,"set yrange [%f:%f]\n",0.0,100e-6*1e6);
    fprintf(script,"set xtics (%f,%f) font \"Bookman URW,16\"\n",1.045,1.065);
    fprintf(script,"set ytics (%f,%f) font \"Bookman URW,16\"\n",0.0,100e-6*1e6);
    fprintf(script,"set xlabel \"\"\n");
    fprintf(script,"set ylabel \"\"\n");
    //fprintf(script,"set format x \"%%.4g\"\n");
    fprintf(script,"unset key\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM2tanCW,datosEM2tanML);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2tan->posKerrIn,EM2tan->posKerrIn);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM2tan->posKerrMid,EM2tan->posKerrMid);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2tan->posKerrOut,EM2tan->posKerrOut);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM2tan->posKerrIn);
    fprintf(script, "set label \'B\' at %Le,GPVAL_Y_MAX*0.9\n",EM2tan->posKerrMid);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM2tan->posKerrOut);
    fprintf(script,"replot\n");
    fprintf(script,"unset multiplot\n");

    // Conf. em2 sagital
    fprintf(script,"unset arrow\n");
    fprintf(script,"unset label\n");
    fprintf(script,"set title \"Propagación sagital desde EM_2\\n{/*0.8 ε_1=-1 [mm], ε_2=-0.9 [mm]}\" font \"Bookman URW,16\" offset 0,0\n");
    fprintf(script,"set output \"EM2sag.pdf\"\n");
    fprintf(script,"set multiplot\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set size 1,1\n");
    fprintf(script,"set bmargin at screen 0.2\n");
    fprintf(script,"set rmargin at screen 0.93\n");
    fprintf(script,"set origin 0,0\n");
    fprintf(script,"set key font \"Bookman URW,16\"\n");
    fprintf(script,"set xlabel \"Posición en la cavidad [m]\"\n");
    fprintf(script,"set ylabel \"Radio de haz ω_s [μm]\"\n");
    fprintf(script,"set xrange [*:*]\n");
    fprintf(script,"set xtics autofreq font \"Bookman URW,16\"\n");
    fprintf(script,"set yrange [*:*]\n");
    fprintf(script,"set ytics 2e-4 font \"Bookman URW,16\"\n");
    //fprintf(script,"set format x \"%%.2g\"\n");
    //fprintf(script,"set format y \"%%.1e\"\n");
    fprintf(script, "set label \'A: Interfaz Ti:Zafiro\' at 0.2,4e-4\n");
    fprintf(script, "set label \'B: Sección media\' at 0.2,3e-4\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM2sagCW,datosEM2sagML);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2sag->posKerrIn,EM2sag->posKerrIn);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM2sag->posKerrMid,EM2sag->posKerrMid);
    //fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2sag->posKerrOut,EM2sag->posKerrOut);
    // Zoom
    fprintf(script,"set size 0.45,0.2\n");
    fprintf(script,"set bmargin at screen 0.3\n");
    fprintf(script,"set rmargin at screen 0.88\n");
    fprintf(script,"set origin 0.5,0.4\n");
    fprintf(script,"set title \"Ti:Zafiro\" offset 0,-1\n");
    fprintf(script,"set xrange [%f:%f]\n",1.045,1.065);
    fprintf(script,"set yrange [%f:%f]\n",0.0,100e-6*1e6);
    fprintf(script,"set xtics (%f,%f) font \"Bookman URW,16\"\n",1.045,1.065);
    fprintf(script,"set ytics (%f,%f) font \"Bookman URW,16\"\n",0.0,100e-6*1e6);
    fprintf(script,"set xlabel \"\"\n");
    fprintf(script,"set ylabel \"\"\n");
    //fprintf(script,"set format x \"%%.4g\"\n");
    fprintf(script,"unset key\n");
    fprintf(script,"plot \"%s\" using 1:($2*1e6) with lines title \"Emisión continua\", \"%s\" using 1:($2*1e6) with lines title \"Emisión pulsada\"\n",datosEM2sagCW,datosEM2sagML);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2sag->posKerrIn,EM2sag->posKerrIn);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'blue\'\n",EM2sag->posKerrMid,EM2sag->posKerrMid);
    fprintf(script,"set arrow from %Le,0 to %Le,GPVAL_Y_MAX*0.75 nohead lc rgb \'red\'\n",EM2sag->posKerrOut,EM2sag->posKerrOut);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM2sag->posKerrIn);
    fprintf(script, "set label \'B\' at %Le,GPVAL_Y_MAX*0.9\n",EM2sag->posKerrMid);
    fprintf(script, "set label \'A\' at %Le,GPVAL_Y_MAX*0.9\n",EM2sag->posKerrOut);
    fprintf(script,"replot\n");
    fprintf(script,"unset multiplot\n");
    fclose(script);
    printf("Fin\n");

}

/* Rutina para gráficar propagación lineal - Conlleva un error por la discreti-
zación de la propagación de rayos (en lugar de una matriz para propagar en espacio,
se utilizan muchas matrices */

void propagacionKerrGrafica(int pasos, long double deltaZ, long double n0,long double n2,
                             long double complex qInTan, long double complex qInSag, long double complex *qTan, \
                                            long double complex *qSag, long double *wTan, long double *wSag, long double chi, long double kth, long double Cp, \
                                            long double rho, long double dn_dv,long double P_laser, long double lambda0, \
                                            ajusteTemperaturaCristal *vectorPlano, _Bool ladoBombeo)
{
        long double wInTan=spot_q(qInTan,n0,lambda0);
        long double wInSag=spot_q(qInSag,n0,lambda0);
        long double complex At, Bt, Ct, Dt;
        long double complex As, Bs, Cs, Ds;
        matriz *MTan, *MSag;
        MTan=nuevaMatriz(2,2);
        MSag=nuevaMatriz(2,2);
        // Cálculo inicial
        long double n0Prop=n0;
        // Identificar si se propaga desde el lado bombeado o no, cambiar coeficientes
        // de ajuste acorde a esto.
        if(ladoBombeo==1)
        {
            // Primer propagación
            wInTan=spot_q(qInTan,n0,lambda0);
            wInSag=spot_q(qInSag,n0,lambda0);
            n0Prop=n0; // Descartados términos lineales por ser pequeños.
            // Formación de matrices
            // Factor de parábola=1/h^2
            long double parabola1=-2.0*vectorPlano->a1[0]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInTan,3)*wInSag);
            long double parabola2=-2.0*vectorPlano->a2[0]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInSag,3)*wInTan);
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
            qTan[0]=prop_q(MTan,qInTan,n0Prop,n0Prop);
            qSag[0]=prop_q(MSag,qInSag,n0Prop,n0Prop);
            wTan[0]=spot_q(qTan[0],n0,lambda0);
            wSag[0]=spot_q(qSag[0],n0,lambda0);
            for (int i=0;i<pasos;i++)
            {
                // Cálculo de spots
                wInTan=spot_q(qTan[i],n0Prop,lambda0);
                wInSag=spot_q(qSag[i],n0Prop,lambda0);
                // Cálculo de plano
                //n0Prop=n0+vectorPlano->a0[i]*dn_dv+n2*P_laser/(M_PI*wTan*wSag)*3.0/4.0; // Factores lineal, térmico y Kerr
                n0Prop=n0; // Descartados términos lineales por ser pequeños.
                // Formación de matrices
                // Factor de parábola=1/h^2
                parabola1=-2.0*vectorPlano->a1[i]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInTan,3)*wInSag);
                parabola2=-2.0*vectorPlano->a2[i]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInSag,3)*wInTan);
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
                    qTan[i+1]=prop_q(MTan,qTan[i],n0Prop,n0Prop);
                    qSag[i+1]=prop_q(MSag,qSag[i],n0Prop,n0Prop);
                    wTan[i+1]=spot_q(qTan[i],n0,lambda0);
                    wSag[i+1]=spot_q(qSag[i],n0,lambda0);
                }
            }
        }
        else
        {
            int indiceTermico=pasos-1;
            long double wInTan=spot_q(qInTan,n0Prop,lambda0);
            long double wInSag=spot_q(qInSag,n0Prop,lambda0);
            // Cálculo de plano
            n0Prop=n0;
            // Formación de matrices
            // Factor de parábola=1/h^2
            long double parabola1=-2.0*vectorPlano->a1[indiceTermico]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInTan,3)*wInSag);
            long double parabola2=-2.0*vectorPlano->a2[indiceTermico]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInSag,3)*wInTan);
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
            qTan[0]=prop_q(MTan,qInTan,n0Prop,n0Prop);
            qSag[0]=prop_q(MSag,qInSag,n0Prop,n0Prop);
            wTan[0]=spot_q(qTan[0],n0,lambda0);
            wSag[0]=spot_q(qSag[0],n0,lambda0);
            indiceTermico--;
            for (int i=0;i<pasos;i++)
            {
                // Cálculo de spots
                wInTan=spot_q(qTan[i],n0Prop,lambda0);
                wInSag=spot_q(qSag[i],n0Prop,lambda0);
                // Cálculo de plano
                n0Prop=n0;
                // Formación de matrices
                // Factor de parábola=1/h^2
                parabola1=-2.0*vectorPlano->a1[indiceTermico]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInTan,3)*wInSag);
                parabola2=-2.0*vectorPlano->a2[indiceTermico]*dn_dv/n0Prop+(n2*P_laser)/(n0Prop*M_PI*powl(wInSag,3)*wInTan);
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
                    qTan[i+1]=prop_q(MTan,qTan[i],n0Prop,n0Prop);
                    qSag[i+1]=prop_q(MSag,qSag[i],n0Prop,n0Prop);
                    wTan[i+1]=spot_q(qTan[i],n0,lambda0);
                    wSag[i+1]=spot_q(qSag[i],n0,lambda0);
                }
                indiceTermico--;
            }
        }
        //Limpieza
        borraMatriz(MTan);
        borraMatriz(MSag);
    }

void graficaPropagacion(char *conjugado_corto, char *conjugado_largo, int iteraciones, long double umbral, int N, long double Epsilon1, long double Epsilon2,
                     long double lambdaPump, long double n0, long double n2, long double P_laser, long double P_pump, long double chi, long double L,
                     long double kth, long double Cp, long double rho, long double dn_dv, long double ancho, long double alto, long double w_pump_t,
                     long double w_pump_s, long double nPump, long double lambda0, long double L1, long double L2,
                     long double f1, long double f2)
{
    // Lente térmica - en el futuro se podrá propagar desde la fuente láser
    long double complex pT,pS,qT,qS;
     long double angulos[2],delta1,delta2,f1t,f1s,f2t,f2s;
    pT=1/100e-3-I*535e-9/(nPump*M_PI*powl(w_pump_t,2));
    qT=1/pT;
    pS=1/100e-3-I*535e-9/(nPump*M_PI*powl(w_pump_s,2));
    qS=1/pS;
    long double alpha=1/1e-2;
    int pasos=1000;
    int pasosKerr=1000;
    long double lambda_relax=1.5; // Valor de relajación para el método Liebmann (Gauss-Seidel) de sol. de EDP.
    ajusteTemperaturaCristal *ajuste=vectorAjusteCristal(qT,qS,alpha,lambdaPump,n0,P_pump,chi,L,kth,Cp,rho,dn_dv,ancho,alto,iteraciones,pasos,umbral,lambda_relax,N);
    printf("Cálculo térmico terminado\n");
    long double paso=(long double)L/pasosKerr;

    // Cálculo CW
    matriz *wt1,*ws1,*wt2,*ws2,*qt1,*qs1,*qt2,*qs2,*epsilon1,*epsilon2;
    matriz *wt1_ml, *ws1_ml, *wt2_ml, *ws2_ml, *astigML1, *astigML2;

    // Llenando epsilon1 y epsilon2
    epsilon1=nuevaMatriz(1,1);
    epsilon2=nuevaMatriz(1,1); // Columnas
    fijaElemento(epsilon1,1,1,Epsilon1);
    fijaElemento(epsilon2,1,1,Epsilon2);

    // Creando matrices receptoras
    wt1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    ws1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    wt2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    ws2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qt1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qs1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qt2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qs2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);

    matriz **qAstigEM1_ML;
    matriz **qAstigEM2_ML;

    // Cálculo lineal
    calculoLineal(conjugado_corto,conjugado_largo,n0,lambda0,L1,L2,L,f1,f2,wt1,ws1,wt2,ws2,epsilon1,epsilon2,qt1,qs1,qt2,qs2);

    // Cálculo no lineal

    qAstigEM1_ML=propNoLinealEM1Astigmatico(qt1,qs1,ajuste,conjugado_corto,conjugado_largo,L1,L2,L,f1,f2,n0,n2,chi,kth,Cp,rho,dn_dv,P_laser,lambda0,epsilon1,epsilon2,iteraciones,pasosKerr,umbral);
    qAstigEM2_ML=propNoLinealEM2Astigmatico(qt2,qs2,ajuste,conjugado_corto,conjugado_largo,L1,L2,L,f1,f2,n0,n2,chi,kth,Cp,rho,dn_dv,P_laser,lambda0,epsilon1,epsilon2,iteraciones,pasosKerr,umbral);

    // Cálculo de spots
    wt1_ml=spot_q_matriz(qAstigEM1_ML[0],1,lambda0);
    ws1_ml=spot_q_matriz(qAstigEM1_ML[1],1,lambda0);
    wt2_ml=spot_q_matriz(qAstigEM2_ML[0],1,lambda0);
    ws2_ml=spot_q_matriz(qAstigEM2_ML[1],1,lambda0);

    // Cálculo de distancias
    anguloLineal("fin","inf",L,n0,L1,L2,f1,f2,angulos);
    delta1=distanciaCristal("fin",f1,L1,L,n0,angulos[0]);
    delta2=distanciaCristal("inf",f2,L2,L,n0,angulos[1]);

    // Cálculo de distancias focales
    f1t=f1*cosl(angulos[0]);
    f2t=f2*cosl(angulos[1]);
    f1s=f1/cosl(angulos[0]);
    f2s=f2/cosl(angulos[1]);

    delta1=delta1+Epsilon1;
    delta2=delta2+Epsilon2;

    // Búsqueda de configuración menos astigmática
    int fila1, columna1,fila2, columna2;
    long double complex valor1, valor2;

    /// Inicio Gráfica
    // matrices a usar
        matriz *propLibre, *espCurv1Tan, *espCurv1Sag, *espCurv2Tan, *espCurv2Sag,
            *brewsterEntradaTan, *brewsterSalidaTan, *brewsterEntradaSag, *brewsterSalidaSag,
            *pasoCristal;

    // Llenado de matrices
    propLibre=llenaMatriz(1,paso,0,1);
    espCurv1Tan=llenaMatriz(1,0,-1/f1t,1);
    espCurv1Sag=llenaMatriz(1,0,-1/f1s,1);
    espCurv2Tan=llenaMatriz(1,0,-1/f2t,1);
    espCurv2Sag=llenaMatriz(1,0,-1/f2s,1);
    brewsterEntradaTan=llenaMatriz(n0,0,0,1/n0);
    brewsterSalidaTan=llenaMatriz(1/n0,0,0,n0);
    brewsterEntradaSag=llenaMatriz(1,0,0,1);
    brewsterSalidaSag=llenaMatriz(1,0,0,1);
    pasoCristal=llenaMatriz(1,paso/n0,0,1);

    /// Espejo EM1

    /// Datos para propagación
    long double complex *q0EM1tanCW, *q1EM1tanCW, *q2EM1tanCW, *q3EM1tanCW, *q4EM1tanCW; // Vectores para almacenar q propagada.
    long double *w0EM1tanCW, *w1EM1tanCW, *w2EM1tanCW, *w3EM1tanCW, *w4EM1tanCW; // Vectores para spot
    long double complex *q0EM1sagCW, *q1EM1sagCW, *q2EM1sagCW, *q3EM1sagCW, *q4EM1sagCW; // Vectores para almacenar q propagada.
    long double *w0EM1sagCW, *w1EM1sagCW, *w2EM1sagCW, *w3EM1sagCW, *w4EM1sagCW; // Vectores para spotlong double complex *q0EM1tan, *q1EM1tan, *q2EM1tan, *q3EM1tan, *q4EM1tan; // Vectores para almacenar q propagada.

    long double complex *q0EM1tanML, *q1EM1tanML, *q2EM1tanML, *q3EM1tanML, *q4EM1tanML; // Vectores para almacenar q propagada.
    long double *w0EM1tanML, *w1EM1tanML, *w2EM1tanML, *w3EM1tanML, *w4EM1tanML; // Vectores para spot
    long double complex *q0EM1sagML, *q1EM1sagML, *q2EM1sagML, *q3EM1sagML, *q4EM1sagML; // Vectores para almacenar q propagada.
    long double *w0EM1sagML, *w1EM1sagML, *w2EM1sagML, *w3EM1sagML, *w4EM1sagML; // Vectores para spotlong double complex *q0EM1tan, *q1EM1tan, *q2EM1tan, *q3EM1tan, *q4EM1tan; // Vectores para almacenar q propagada.
    // Estimación de tamaño y reserva de memoriade vector q0.
    int contador=0;
    int numeroPasosq0EM1 = L1/paso+1;  // Incluye reflexión en espejo curvo

    q0EM1tanCW=(long double complex*)calloc(numeroPasosq0EM1,sizeof(long double complex));
    w0EM1tanCW=(long double*)calloc(numeroPasosq0EM1,sizeof(long double));
    q0EM1sagCW=(long double complex*)calloc(numeroPasosq0EM1,sizeof(long double complex));
    w0EM1sagCW=(long double*)calloc(numeroPasosq0EM1,sizeof(long double));

    q0EM1tanML=(long double complex*)calloc(numeroPasosq0EM1,sizeof(long double complex));
    w0EM1tanML=(long double*)calloc(numeroPasosq0EM1,sizeof(long double));
    q0EM1sagML=(long double complex*)calloc(numeroPasosq0EM1,sizeof(long double complex));
    w0EM1sagML=(long double*)calloc(numeroPasosq0EM1,sizeof(long double));


    // Propagación en q0EM1

    q0EM1tanCW[0]=obtieneElemento(qt1,1,1);
    w0EM1tanCW[0]=spot_q(q0EM1tanCW[0],1,lambda0);
    q0EM1sagCW[0]=obtieneElemento(qs1,1,1);
    w0EM1sagCW[0]=spot_q(q0EM1sagCW[0],1,lambda0);

    q0EM1tanML[0]=obtieneElemento(qAstigEM1_ML[0],1,1);
    w0EM1tanML[0]=spot_q(q0EM1tanML[0],1,lambda0);
    q0EM1sagML[0]=obtieneElemento(qAstigEM1_ML[1],1,1);
    w0EM1sagML[0]=spot_q(q0EM1sagML[0],1,lambda0);

    for(contador=1;contador<numeroPasosq0EM1-1;contador++)
    {
        q0EM1tanCW[contador]=prop_q(propLibre,q0EM1tanCW[contador-1],1,1);
        w0EM1tanCW[contador]=spot_q(q0EM1tanCW[contador],1,lambda0);
        q0EM1sagCW[contador]=prop_q(propLibre,q0EM1sagCW[contador-1],1,1);
        w0EM1sagCW[contador]=spot_q(q0EM1sagCW[contador],1,lambda0);

        q0EM1tanML[contador]=prop_q(propLibre,q0EM1tanML[contador-1],1,1);
        w0EM1tanML[contador]=spot_q(q0EM1tanML[contador],1,lambda0);
        q0EM1sagML[contador]=prop_q(propLibre,q0EM1sagML[contador-1],1,1);
        w0EM1sagML[contador]=spot_q(q0EM1sagML[contador],1,lambda0);
    }
    contador--;
    q0EM1tanCW[contador+1]=prop_q(espCurv1Tan,q0EM1tanCW[contador],1,1);
    w0EM1tanCW[contador+1]=spot_q(q0EM1tanCW[contador+1],1,lambda0);
    q0EM1sagCW[contador+1]=prop_q(espCurv1Sag,q0EM1sagCW[contador],1,1);
    w0EM1sagCW[contador+1]=spot_q(q0EM1sagCW[contador+1],1,lambda0);
    q0EM1tanML[contador+1]=prop_q(espCurv1Tan,q0EM1tanML[contador],1,1);
    w0EM1tanML[contador+1]=spot_q(q0EM1tanML[contador+1],1,lambda0);
    q0EM1sagML[contador+1]=prop_q(espCurv1Sag,q0EM1sagML[contador],1,1);
    w0EM1sagML[contador+1]=spot_q(q0EM1sagML[contador+1],1,lambda0);

    // Creación de q1
    int numeroPasosq1EM1 = delta1/paso+1; // incluye refracción en cristal
    q1EM1tanCW=(long double complex*)calloc(numeroPasosq1EM1,sizeof(long double complex));
    w1EM1tanCW=(long double*)calloc(numeroPasosq1EM1,sizeof(long double));
    q1EM1sagCW=(long double complex*)calloc(numeroPasosq1EM1,sizeof(long double complex));
    w1EM1sagCW=(long double*)calloc(numeroPasosq1EM1,sizeof(long double));
    q1EM1tanML=(long double complex*)calloc(numeroPasosq1EM1,sizeof(long double complex));
    w1EM1tanML=(long double*)calloc(numeroPasosq1EM1,sizeof(long double));
    q1EM1sagML=(long double complex*)calloc(numeroPasosq1EM1,sizeof(long double complex));
    w1EM1sagML=(long double*)calloc(numeroPasosq1EM1,sizeof(long double));
    // Propagación en q1EM1tan
    q1EM1tanCW[0]=prop_q(propLibre,q0EM1tanCW[contador+1],1,1);
    w1EM1tanCW[0]=spot_q(q1EM1tanCW[0],1,lambda0);
    q1EM1sagCW[0]=prop_q(propLibre,q0EM1sagCW[contador+1],1,1);
    w1EM1sagCW[0]=spot_q(q1EM1sagCW[0],1,lambda0);
    q1EM1tanML[0]=prop_q(propLibre,q0EM1tanML[contador+1],1,1);
    w1EM1tanML[0]=spot_q(q1EM1tanML[0],1,lambda0);
    q1EM1sagML[0]=prop_q(propLibre,q0EM1sagML[contador+1],1,1);
    w1EM1sagML[0]=spot_q(q1EM1sagML[0],1,lambda0);


    for(contador=1;contador<numeroPasosq1EM1-1;contador++)
    {
        q1EM1tanCW[contador]=prop_q(propLibre,q1EM1tanCW[contador-1],1,1);
        w1EM1tanCW[contador]=spot_q(q1EM1tanCW[contador],1,lambda0);
        q1EM1sagCW[contador]=prop_q(propLibre,q1EM1sagCW[contador-1],1,1);
        w1EM1sagCW[contador]=spot_q(q1EM1sagCW[contador],1,lambda0);
        q1EM1tanML[contador]=prop_q(propLibre,q1EM1tanML[contador-1],1,1);
        w1EM1tanML[contador]=spot_q(q1EM1tanML[contador],1,lambda0);
        q1EM1sagML[contador]=prop_q(propLibre,q1EM1sagML[contador-1],1,1);
        w1EM1sagML[contador]=spot_q(q1EM1sagML[contador],1,lambda0);
    }
    contador--;
    q1EM1tanCW[contador+1]=prop_q(brewsterEntradaTan,q1EM1tanCW[contador],1,n0);
    w1EM1tanCW[contador+1]=spot_q(q1EM1tanCW[contador+1],n0,lambda0);
    q1EM1sagCW[contador+1]=prop_q(brewsterEntradaSag,q1EM1sagCW[contador],1,n0);
    w1EM1sagCW[contador+1]=spot_q(q1EM1sagCW[contador+1],n0,lambda0);
    q1EM1tanML[contador+1]=prop_q(brewsterEntradaTan,q1EM1tanML[contador],1,n0);
    w1EM1tanML[contador+1]=spot_q(q1EM1tanML[contador+1],n0,lambda0);
    q1EM1sagML[contador+1]=prop_q(brewsterEntradaSag,q1EM1sagML[contador],1,n0);
    w1EM1sagML[contador+1]=spot_q(q1EM1sagML[contador+1],n0,lambda0);
    int posProp=contador+1;

    // Creación de q2
    int numeroPasosq2EM1 = pasos+1; // incluye refracción hacia afuera del cristal
    q2EM1tanCW=(long double complex*)calloc(numeroPasosq2EM1,sizeof(long double complex));
    w2EM1tanCW=(long double *)calloc(numeroPasosq2EM1,sizeof(long double));
    q2EM1sagCW=(long double complex*)calloc(numeroPasosq2EM1,sizeof(long double complex));
    w2EM1sagCW=(long double *)calloc(numeroPasosq2EM1,sizeof(long double));
    q2EM1tanML=(long double complex*)calloc(numeroPasosq2EM1,sizeof(long double complex));
    w2EM1tanML=(long double *)calloc(numeroPasosq2EM1,sizeof(long double));
    q2EM1sagML=(long double complex*)calloc(numeroPasosq2EM1,sizeof(long double complex));
    w2EM1sagML=(long double *)calloc(numeroPasosq2EM1,sizeof(long double));


    /// Propagación en cristal

    q2EM1tanCW[0]=prop_q(pasoCristal,q1EM1tanCW[contador+1],n0,n0);
    w2EM1tanCW[0]=spot_q(q2EM1tanCW[0],n0,lambda0);
    q2EM1sagCW[0]=prop_q(pasoCristal,q1EM1sagCW[contador+1],n0,n0);
    w2EM1sagCW[0]=spot_q(q2EM1sagCW[0],n0,lambda0);
    //q2EM1tanML[0]=prop_q(pasoCristal,q1EM1tanML[contador+1],n0,n0);
    //w2EM1tanML[0]=spot_q(q2EM1tanML[0],n0,lambda0);
    //q2EM1sagML[0]=prop_q(pasoCristal,q1EM1sagML[contador+1],n0,n0);
    //w2EM1sagML[0]=spot_q(q2EM1sagML[0],n0,lambda0);
    for(contador=1;contador<numeroPasosq2EM1-1;contador++)
    {
        q2EM1tanCW[contador]=prop_q(pasoCristal,q2EM1tanCW[contador-1],n0,n0);
        w2EM1tanCW[contador]=spot_q(q2EM1tanCW[contador],n0,lambda0);
        q2EM1sagCW[contador]=prop_q(pasoCristal,q2EM1sagCW[contador-1],n0,n0);
        w2EM1sagCW[contador]=spot_q(q2EM1sagCW[contador],n0,lambda0);
    }
    propagacionKerrGrafica(numeroPasosq2EM1-1,paso,n0,n2,q1EM1tanML[posProp],q1EM1sagML[posProp],q2EM1tanML,q2EM1sagML,w2EM1tanML,w2EM1sagML,chi,kth,Cp,rho,dn_dv,P_laser,lambda0,
                           ajuste,"1");
    //contador=999; // numeroPasosq2EM1-1
    contador--;
    q2EM1tanCW[contador+1]=prop_q(brewsterSalidaTan,q2EM1tanCW[contador],n0,1);
    w2EM1tanCW[contador+1]=spot_q(q2EM1tanCW[contador+1],1,lambda0);
    q2EM1sagCW[contador+1]=prop_q(brewsterSalidaSag,q2EM1sagCW[contador],n0,1);
    w2EM1sagCW[contador+1]=spot_q(q2EM1sagCW[contador+1],1,lambda0);
    q2EM1tanML[contador+1]=prop_q(brewsterSalidaTan,q2EM1tanML[contador],n0,1);
    w2EM1tanML[contador+1]=spot_q(q2EM1tanML[contador+1],1,lambda0);
    q2EM1sagML[contador+1]=prop_q(brewsterSalidaSag,q2EM1sagML[contador],n0,1);
    w2EM1sagML[contador+1]=spot_q(q2EM1sagML[contador+1],1,lambda0);

    // Creación de q3
    int numeroPasosq3EM1 = delta2/paso+1; // Incluye reflexión en espejo curvo.
    q3EM1tanCW=(long double complex*)calloc(numeroPasosq3EM1, sizeof(long double complex));
    w3EM1tanCW=(long double*)calloc(numeroPasosq3EM1,sizeof(long double));
    q3EM1sagCW=(long double complex*)calloc(numeroPasosq3EM1, sizeof(long double complex));
    w3EM1sagCW=(long double*)calloc(numeroPasosq3EM1,sizeof(long double));q3EM1tanML=(long double complex*)calloc(numeroPasosq3EM1, sizeof(long double complex));
    w3EM1tanML=(long double*)calloc(numeroPasosq3EM1,sizeof(long double));
    q3EM1sagML=(long double complex*)calloc(numeroPasosq3EM1, sizeof(long double complex));
    w3EM1sagML=(long double*)calloc(numeroPasosq3EM1,sizeof(long double));

    // Propagación de q3
    q3EM1tanCW[0]=prop_q(propLibre,q2EM1tanCW[contador+1],1,1);
    w3EM1tanCW[0]=spot_q(q3EM1tanCW[0],1,lambda0);
    q3EM1sagCW[0]=prop_q(propLibre,q2EM1sagCW[contador+1],1,1);
    w3EM1sagCW[0]=spot_q(q3EM1sagCW[0],1,lambda0);
    q3EM1tanML[0]=prop_q(propLibre,q2EM1tanML[contador+1],1,1);
    w3EM1tanML[0]=spot_q(q3EM1tanML[0],1,lambda0);
    q3EM1sagML[0]=prop_q(propLibre,q2EM1sagML[contador+1],1,1);
    w3EM1sagML[0]=spot_q(q3EM1sagML[0],1,lambda0);

    for(contador=1;contador<numeroPasosq3EM1-1;contador++)
    {
        q3EM1tanCW[contador]=prop_q(propLibre,q3EM1tanCW[contador-1],1,1);
        w3EM1tanCW[contador]=spot_q(q3EM1tanCW[contador],1,lambda0);
        q3EM1sagCW[contador]=prop_q(propLibre,q3EM1sagCW[contador-1],1,1);
        w3EM1sagCW[contador]=spot_q(q3EM1sagCW[contador],1,lambda0);
        q3EM1tanML[contador]=prop_q(propLibre,q3EM1tanML[contador-1],1,1);
        w3EM1tanML[contador]=spot_q(q3EM1tanML[contador],1,lambda0);
        q3EM1sagML[contador]=prop_q(propLibre,q3EM1sagML[contador-1],1,1);
        w3EM1sagML[contador]=spot_q(q3EM1sagML[contador],1,lambda0);
    }
    contador--;
    q3EM1tanCW[contador+1]=prop_q(espCurv2Tan,q3EM1tanCW[contador],1,1);
    w3EM1tanCW[contador+1]=spot_q(q3EM1tanCW[contador+1],1,lambda0);
    q3EM1sagCW[contador+1]=prop_q(espCurv2Sag,q3EM1sagCW[contador],1,1);
    w3EM1sagCW[contador+1]=spot_q(q3EM1sagCW[contador+1],1,lambda0);
    q3EM1tanML[contador+1]=prop_q(espCurv2Tan,q3EM1tanML[contador],1,1);
    w3EM1tanML[contador+1]=spot_q(q3EM1tanML[contador+1],1,lambda0);
    q3EM1sagML[contador+1]=prop_q(espCurv2Sag,q3EM1sagML[contador],1,1);
    w3EM1sagML[contador+1]=spot_q(q3EM1sagML[contador+1],1,lambda0);

    // Creación de q4
    int numeroPasosq4EM1= L2/paso;
    q4EM1tanCW=(long double complex*)calloc(numeroPasosq4EM1, sizeof(long double complex));
    w4EM1tanCW=(long double *)calloc(numeroPasosq4EM1,sizeof(long double));
    q4EM1sagCW=(long double complex*)calloc(numeroPasosq4EM1, sizeof(long double complex));
    w4EM1sagCW=(long double *)calloc(numeroPasosq4EM1,sizeof(long double));
    q4EM1tanML=(long double complex*)calloc(numeroPasosq4EM1, sizeof(long double complex));
    w4EM1tanML=(long double *)calloc(numeroPasosq4EM1,sizeof(long double));
    q4EM1sagML=(long double complex*)calloc(numeroPasosq4EM1, sizeof(long double complex));
    w4EM1sagML=(long double *)calloc(numeroPasosq4EM1,sizeof(long double));
    // propagación de q4
    q4EM1tanCW[0]=prop_q(propLibre,q3EM1tanCW[contador+1],1,1);
    w4EM1tanCW[0]=spot_q(q4EM1tanCW[0],1,lambda0);
    q4EM1sagCW[0]=prop_q(propLibre,q3EM1sagCW[contador+1],1,1);
    w4EM1sagCW[0]=spot_q(q4EM1sagCW[0],1,lambda0);
    q4EM1tanML[0]=prop_q(propLibre,q3EM1tanML[contador+1],1,1);
    w4EM1tanML[0]=spot_q(q4EM1tanML[0],1,lambda0);
    q4EM1sagML[0]=prop_q(propLibre,q3EM1sagML[contador+1],1,1);
    w4EM1sagML[0]=spot_q(q4EM1sagML[0],1,lambda0);
    for(contador=1;contador<numeroPasosq4EM1-1;contador++)
    {
        q4EM1tanCW[contador]=prop_q(propLibre,q4EM1tanCW[contador-1],1,1);
        w4EM1tanCW[contador]=spot_q(q4EM1tanCW[contador],1,lambda0);
        q4EM1sagCW[contador]=prop_q(propLibre,q4EM1sagCW[contador-1],1,1);
        w4EM1sagCW[contador]=spot_q(q4EM1sagCW[contador],1,lambda0);
        q4EM1tanML[contador]=prop_q(propLibre,q4EM1tanML[contador-1],1,1);
        w4EM1tanML[contador]=spot_q(q4EM1tanML[contador],1,lambda0);
        q4EM1sagML[contador]=prop_q(propLibre,q4EM1sagML[contador-1],1,1);
        w4EM1sagML[contador]=spot_q(q4EM1sagML[contador],1,lambda0);
    }
    contador--;
    q4EM1tanCW[contador+1]=prop_q(propLibre,q4EM1tanCW[contador],1,1);
    w4EM1tanCW[contador+1]=spot_q(q4EM1tanCW[contador+1],1,lambda0);
    q4EM1sagCW[contador+1]=prop_q(propLibre,q4EM1sagCW[contador],1,1);
    w4EM1sagCW[contador+1]=spot_q(q4EM1sagCW[contador+1],1,lambda0);
    q4EM1tanML[contador+1]=prop_q(propLibre,q4EM1tanML[contador],1,1);
    w4EM1tanML[contador+1]=spot_q(q4EM1tanML[contador+1],1,lambda0);
    q4EM1sagML[contador+1]=prop_q(propLibre,q4EM1sagML[contador],1,1);
    w4EM1sagML[contador+1]=spot_q(q4EM1sagML[contador+1],1,lambda0);


    /// Escritura EM1

    /// Estructura para guardar
    vectorDato *EM1tan=(vectorDato *)calloc(1,sizeof(vectorDato));
    EM1tan->numeroElementos=numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1+numeroPasosq4EM1-4; // Ultima propagación no hay refracción.
    EM1tan->pos=(long double *)calloc(EM1tan->numeroElementos,sizeof(long double));
    EM1tan->spotCW=(long double *)calloc(EM1tan->numeroElementos,sizeof(long double));
    EM1tan->spotML=(long double *)calloc(EM1tan->numeroElementos,sizeof(long double));
    vectorDato *EM1sag=(vectorDato *)calloc(1,sizeof(vectorDato));
    EM1sag->numeroElementos=numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1+numeroPasosq4EM1-4; // Ultima propagación no hay refracción.
    EM1sag->pos=(long double *)calloc(EM1sag->numeroElementos,sizeof(long double));
    EM1sag->spotCW=(long double *)calloc(EM1sag->numeroElementos,sizeof(long double));
    EM1sag->spotML=(long double *)calloc(EM1sag->numeroElementos,sizeof(long double));


    for(int i=0;i<numeroPasosq0EM1-1;i++)
    {
        EM1tan->pos[i]=i*paso;
        EM1tan->spotCW[i]=w0EM1tanCW[i];
        EM1tan->spotML[i]=w0EM1tanML[i];
        EM1sag->pos[i]=i*paso;
        EM1sag->spotCW[i]=w0EM1sagCW[i];
        EM1sag->spotML[i]=w0EM1sagML[i];
    }
    for(int i=0;i<numeroPasosq1EM1-1;i++)
    {
        EM1tan->pos[i+numeroPasosq0EM1-1]=(i+numeroPasosq0EM1-1)*paso;
        EM1tan->spotCW[i+numeroPasosq0EM1-1]=w1EM1tanCW[i];
        EM1tan->spotML[i+numeroPasosq0EM1-1]=w1EM1tanML[i];
        EM1sag->pos[i+numeroPasosq0EM1-1]=(i+numeroPasosq0EM1-1)*paso;
        EM1sag->spotCW[i+numeroPasosq0EM1-1]=w1EM1sagCW[i];
        EM1sag->spotML[i+numeroPasosq0EM1-1]=w1EM1sagML[i];
    }
    for(int i=0;i<numeroPasosq2EM1-1;i++)
    {
        if(i==0)
        {
            EM1tan->posKerrIn=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso;
            EM1sag->posKerrIn=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso; // Inicio de cristal
        }
        if(i==(int)(numeroPasosq2EM1-1)/2)
        {
            EM1tan->posKerrMid=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso; // Mitad de cristal
            EM1sag->posKerrMid=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso; // Mitad de cristal
        }
        if(i==(numeroPasosq2EM1-2))
        {
            EM1tan->posKerrOut=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso; // Fin de cristal
            EM1sag->posKerrOut=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso; // Fin de cristal
        }
        EM1tan->pos[i+numeroPasosq0EM1+numeroPasosq1EM1-2]=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso;
        EM1tan->spotCW[i+numeroPasosq0EM1+numeroPasosq1EM1-2]=w2EM1tanCW[i];
        EM1tan->spotML[i+numeroPasosq0EM1+numeroPasosq1EM1-2]=w2EM1tanML[i];
        EM1sag->pos[i+numeroPasosq0EM1+numeroPasosq1EM1-2]=(i+numeroPasosq0EM1+numeroPasosq1EM1-2)*paso;
        EM1sag->spotCW[i+numeroPasosq0EM1+numeroPasosq1EM1-2]=w2EM1sagCW[i];
        EM1sag->spotML[i+numeroPasosq0EM1+numeroPasosq1EM1-2]=w2EM1sagML[i];
    }
    for(int i=0;i<numeroPasosq3EM1-1;i++)
    {
        EM1tan->pos[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3]=(i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3)*paso;
        EM1tan->spotCW[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3]=w3EM1tanCW[i];
        EM1tan->spotML[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3]=w3EM1tanML[i];
        EM1sag->pos[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3]=(i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3)*paso;
        EM1sag->spotCW[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3]=w3EM1sagCW[i];
        EM1sag->spotML[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1-3]=w3EM1sagML[i];
    }
    for(int i=0;i<numeroPasosq4EM1-1;i++)
    {
        EM1tan->pos[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4]=(i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4)*paso;
        EM1tan->spotCW[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4]=w4EM1tanCW[i];
        EM1tan->spotML[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4]=w4EM1tanML[i];
        EM1sag->pos[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4]=(i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4)*paso;
        EM1sag->spotCW[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4]=w4EM1sagCW[i];
        EM1sag->spotML[i+numeroPasosq0EM1+numeroPasosq1EM1+numeroPasosq2EM1+numeroPasosq3EM1-4]=w4EM1sagML[i];
    }

    // Limpieza de vectores
    free(q0EM1tanCW);
    free(q1EM1tanCW);
    free(q2EM1tanCW);
    free(q3EM1tanCW);
    free(q4EM1tanCW);
    free(w0EM1tanCW);
    free(w1EM1tanCW);
    free(w2EM1tanCW);
    free(w3EM1tanCW);
    free(w4EM1tanCW);
    free(q0EM1tanML);
    free(q1EM1tanML);
    free(q2EM1tanML);
    free(q3EM1tanML);
    free(q4EM1tanML);
    free(w0EM1tanML);
    free(w1EM1tanML);
    free(w2EM1tanML);
    free(w3EM1tanML);
    free(w4EM1tanML);
    free(q0EM1sagCW);
    free(q1EM1sagCW);
    free(q2EM1sagCW);
    free(q3EM1sagCW);
    free(q4EM1sagCW);
    free(w0EM1sagCW);
    free(w1EM1sagCW);
    free(w2EM1sagCW);
    free(w3EM1sagCW);
    free(w4EM1sagCW);
    free(q0EM1sagML);
    free(q1EM1sagML);
    free(q2EM1sagML);
    free(q3EM1sagML);
    free(q4EM1sagML);
    free(w0EM1sagML);
    free(w1EM1sagML);
    free(w2EM1sagML);
    free(w3EM1sagML);
    free(w4EM1sagML);

    /// FIN EM1
    /// INICIO EM2

    long double complex *q0EM2tanCW, *q1EM2tanCW, *q2EM2tanCW, *q3EM2tanCW, *q4EM2tanCW; // Vectores para almacenar q propagada.
    long double *w0EM2tanCW, *w1EM2tanCW, *w2EM2tanCW, *w3EM2tanCW, *w4EM2tanCW; // Vectores para spot
    long double complex *q0EM2sagCW, *q1EM2sagCW, *q2EM2sagCW, *q3EM2sagCW, *q4EM2sagCW; // Vectores para almacenar q propagada.
    long double *w0EM2sagCW, *w1EM2sagCW, *w2EM2sagCW, *w3EM2sagCW, *w4EM2sagCW; // Vectores para spot
    long double complex *q0EM2tanML, *q1EM2tanML, *q2EM2tanML, *q3EM2tanML, *q4EM2tanML; // Vectores para almacenar q propagada.
    long double *w0EM2tanML, *w1EM2tanML, *w2EM2tanML, *w3EM2tanML, *w4EM2tanML; // Vectores para spot
    long double complex *q0EM2sagML, *q1EM2sagML, *q2EM2sagML, *q3EM2sagML, *q4EM2sagML; // Vectores para almacenar q propagada.
    long double *w0EM2sagML, *w1EM2sagML, *w2EM2sagML, *w3EM2sagML, *w4EM2sagML; // Vectores para spot

    // Estimación de tamaño y reserva de memoriade vector q0.

    int numeroPasosq0EM2 = L2/paso+1;  // Incluye reflexión en espejo curvo

    q0EM2tanCW=(long double complex*)calloc(numeroPasosq0EM2,sizeof(long double complex));
    w0EM2tanCW=(long double*)calloc(numeroPasosq0EM2,sizeof(long double));
    q0EM2sagCW=(long double complex*)calloc(numeroPasosq0EM2,sizeof(long double complex));
    w0EM2sagCW=(long double*)calloc(numeroPasosq0EM2,sizeof(long double));
    q0EM2tanML=(long double complex*)calloc(numeroPasosq0EM2,sizeof(long double complex));
    w0EM2tanML=(long double*)calloc(numeroPasosq0EM2,sizeof(long double));
    q0EM2sagML=(long double complex*)calloc(numeroPasosq0EM2,sizeof(long double complex));
    w0EM2sagML=(long double*)calloc(numeroPasosq0EM2,sizeof(long double));


    // Propagación en q0EM2

    q0EM2tanCW[0]=obtieneElemento(qt2,1,1);
    w0EM2tanCW[0]=spot_q(q0EM2tanCW[0],1,lambda0);
    q0EM2sagCW[0]=obtieneElemento(qs2,1,1);
    w0EM2sagCW[0]=spot_q(q0EM2sagCW[0],1,lambda0);
    q0EM2tanML[0]=obtieneElemento(qAstigEM2_ML[0],1,1);
    w0EM2tanML[0]=spot_q(q0EM2tanML[0],1,lambda0);
    q0EM2sagML[0]=obtieneElemento(qAstigEM2_ML[1],1,1);
    w0EM2sagML[0]=spot_q(q0EM2sagML[0],1,lambda0);

    for(contador=1;contador<numeroPasosq0EM2-1;contador++)
    {
        q0EM2tanCW[contador]=prop_q(propLibre,q0EM2tanCW[contador-1],1,1);
        w0EM2tanCW[contador]=spot_q(q0EM2tanCW[contador],1,lambda0);
        q0EM2sagCW[contador]=prop_q(propLibre,q0EM2sagCW[contador-1],1,1);
        w0EM2sagCW[contador]=spot_q(q0EM2sagCW[contador],1,lambda0);
        q0EM2tanML[contador]=prop_q(propLibre,q0EM2tanML[contador-1],1,1);
        w0EM2tanML[contador]=spot_q(q0EM2tanML[contador],1,lambda0);
        q0EM2sagML[contador]=prop_q(propLibre,q0EM2sagML[contador-1],1,1);
        w0EM2sagML[contador]=spot_q(q0EM2sagML[contador],1,lambda0);
    }
    contador--;
    q0EM2tanCW[contador+1]=prop_q(espCurv1Tan,q0EM2tanCW[contador],1,1);
    w0EM2tanCW[contador+1]=spot_q(q0EM2tanCW[contador+1],1,lambda0);
    q0EM2sagCW[contador+1]=prop_q(espCurv1Sag,q0EM2sagCW[contador],1,1);
    w0EM2sagCW[contador+1]=spot_q(q0EM2sagCW[contador+1],1,lambda0);
    q0EM2tanML[contador+1]=prop_q(espCurv1Tan,q0EM2tanML[contador],1,1);
    w0EM2tanML[contador+1]=spot_q(q0EM2tanML[contador+1],1,lambda0);
    q0EM2sagML[contador+1]=prop_q(espCurv1Sag,q0EM2sagML[contador],1,1);
    w0EM2sagML[contador+1]=spot_q(q0EM2sagML[contador+1],1,lambda0);

    // Creación de q1
    int numeroPasosq1EM2 = delta2/paso+1; // incluye refracción en cristal
    q1EM2tanCW=(long double complex*)calloc(numeroPasosq1EM2,sizeof(long double complex));
    w1EM2tanCW=(long double*)calloc(numeroPasosq1EM2,sizeof(long double));
    q1EM2sagCW=(long double complex*)calloc(numeroPasosq1EM2,sizeof(long double complex));
    w1EM2sagCW=(long double*)calloc(numeroPasosq1EM2,sizeof(long double));
    q1EM2tanML=(long double complex*)calloc(numeroPasosq1EM2,sizeof(long double complex));
    w1EM2tanML=(long double*)calloc(numeroPasosq1EM2,sizeof(long double));
    q1EM2sagML=(long double complex*)calloc(numeroPasosq1EM2,sizeof(long double complex));
    w1EM2sagML=(long double*)calloc(numeroPasosq1EM2,sizeof(long double));
    // Propagación en q1EM2tanML
    q1EM2tanCW[0]=prop_q(propLibre,q0EM2tanCW[contador+1],1,1);
    w1EM2tanCW[0]=spot_q(q1EM2tanCW[0],1,lambda0);
    q1EM2sagCW[0]=prop_q(propLibre,q0EM2sagCW[contador+1],1,1);
    w1EM2sagCW[0]=spot_q(q1EM2sagCW[0],1,lambda0);
    q1EM2tanML[0]=prop_q(propLibre,q0EM2tanML[contador+1],1,1);
    w1EM2tanML[0]=spot_q(q1EM2tanML[0],1,lambda0);
    q1EM2sagML[0]=prop_q(propLibre,q0EM2sagML[contador+1],1,1);
    w1EM2sagML[0]=spot_q(q1EM2sagML[0],1,lambda0);


    for(contador=1;contador<numeroPasosq1EM2-1;contador++)
    {
        q1EM2tanCW[contador]=prop_q(propLibre,q1EM2tanCW[contador-1],1,1);
        w1EM2tanCW[contador]=spot_q(q1EM2tanCW[contador],1,lambda0);
        q1EM2sagCW[contador]=prop_q(propLibre,q1EM2sagCW[contador-1],1,1);
        w1EM2sagCW[contador]=spot_q(q1EM2sagCW[contador],1,lambda0);
        q1EM2tanML[contador]=prop_q(propLibre,q1EM2tanML[contador-1],1,1);
        w1EM2tanML[contador]=spot_q(q1EM2tanML[contador],1,lambda0);
        q1EM2sagML[contador]=prop_q(propLibre,q1EM2sagML[contador-1],1,1);
        w1EM2sagML[contador]=spot_q(q1EM2sagML[contador],1,lambda0);
    }
    contador--;
    q1EM2tanCW[contador+1]=prop_q(brewsterEntradaTan,q1EM2tanCW[contador],1,n0);
    w1EM2tanCW[contador+1]=spot_q(q1EM2tanCW[contador+1],n0,lambda0);
    q1EM2sagCW[contador+1]=prop_q(brewsterEntradaSag,q1EM2sagCW[contador],1,n0);
    w1EM2sagCW[contador+1]=spot_q(q1EM2sagCW[contador+1],n0,lambda0);
    q1EM2tanML[contador+1]=prop_q(brewsterEntradaTan,q1EM2tanML[contador],1,n0);
    w1EM2tanML[contador+1]=spot_q(q1EM2tanML[contador+1],n0,lambda0);
    q1EM2sagML[contador+1]=prop_q(brewsterEntradaSag,q1EM2sagML[contador],1,n0);
    w1EM2sagML[contador+1]=spot_q(q1EM2sagML[contador+1],n0,lambda0);
    posProp=contador+1;

    // Creación de q2
    int numeroPasosq2EM2 = pasos+1; // incluye refracción hacia afuera del cristal
    q2EM2tanCW=(long double complex*)calloc(numeroPasosq2EM2,sizeof(long double complex));
    w2EM2tanCW=(long double *)calloc(numeroPasosq2EM2,sizeof(long double));
    q2EM2sagCW=(long double complex*)calloc(numeroPasosq2EM2,sizeof(long double complex));
    w2EM2sagCW=(long double *)calloc(numeroPasosq2EM2,sizeof(long double));
    q2EM2tanML=(long double complex*)calloc(numeroPasosq2EM2,sizeof(long double complex));
    w2EM2tanML=(long double *)calloc(numeroPasosq2EM2,sizeof(long double));
    q2EM2sagML=(long double complex*)calloc(numeroPasosq2EM2,sizeof(long double complex));
    w2EM2sagML=(long double *)calloc(numeroPasosq2EM2,sizeof(long double));


    /// Propagación en cristal
    q2EM2tanCW[0]=prop_q(pasoCristal,q1EM2tanCW[contador+1],n0,n0);
    w2EM2tanCW[0]=spot_q(q2EM2tanCW[0],n0,lambda0);
    q2EM2sagCW[0]=prop_q(pasoCristal,q1EM2sagCW[contador+1],n0,n0);
    w2EM2sagCW[0]=spot_q(q2EM2sagCW[0],n0,lambda0);
    //q2EM2tanML[0]=prop_q(pasoCristal,q1EM2tanML[contador+1],n0,n0);
    //w2EM2tanML[0]=spot_q(q2EM2tanML[0],n0,lambda0);
    //q2EM2sagML[0]=prop_q(pasoCristal,q1EM2sagML[contador+1],n0,n0);
    //w2EM2sagML[0]=spot_q(q2EM2sagML[0],n0,lambda0);
    for(contador=1;contador<numeroPasosq2EM2-1;contador++)
    {
        q2EM2tanCW[contador]=prop_q(pasoCristal,q2EM2tanCW[contador-1],n0,n0);
        w2EM2tanCW[contador]=spot_q(q2EM2tanCW[contador],n0,lambda0);
        q2EM2sagCW[contador]=prop_q(pasoCristal,q2EM2sagCW[contador-1],n0,n0);
        w2EM2sagCW[contador]=spot_q(q2EM2sagCW[contador],n0,lambda0);
    }
    propagacionKerrGrafica(numeroPasosq2EM2-1,paso,n0,n2,q1EM2tanML[posProp],q1EM2sagML[posProp],q2EM2tanML,q2EM2sagML,w2EM2tanML,w2EM2sagML,chi,kth,Cp,rho,dn_dv,P_laser,lambda0,
                           ajuste,"0");
    //contador=999; // numeroPasosq2EM2-1
    contador--;
    q2EM2tanCW[contador+1]=prop_q(brewsterSalidaTan,q2EM2tanCW[contador],n0,1);
    w2EM2tanCW[contador+1]=spot_q(q2EM2tanCW[contador+1],1,lambda0);
    q2EM2sagCW[contador+1]=prop_q(brewsterSalidaSag,q2EM2sagCW[contador],n0,1);
    w2EM2sagCW[contador+1]=spot_q(q2EM2sagCW[contador+1],1,lambda0);
    q2EM2tanML[contador+1]=prop_q(brewsterSalidaTan,q2EM2tanML[contador],n0,1);
    w2EM2tanML[contador+1]=spot_q(q2EM2tanML[contador+1],1,lambda0);
    q2EM2sagML[contador+1]=prop_q(brewsterSalidaSag,q2EM2sagML[contador],n0,1);
    w2EM2sagML[contador+1]=spot_q(q2EM2sagML[contador+1],1,lambda0);

    // Creación de q3
    int numeroPasosq3EM2 = delta1/paso+1; // Incluye reflexión en espejo curvo.
    q3EM2tanCW=(long double complex*)calloc(numeroPasosq3EM2, sizeof(long double complex));
    w3EM2tanCW=(long double*)calloc(numeroPasosq3EM2,sizeof(long double));
    q3EM2sagCW=(long double complex*)calloc(numeroPasosq3EM2, sizeof(long double complex));
    w3EM2sagCW=(long double*)calloc(numeroPasosq3EM2,sizeof(long double));
    q3EM2tanML=(long double complex*)calloc(numeroPasosq3EM2, sizeof(long double complex));
    w3EM2tanML=(long double*)calloc(numeroPasosq3EM2,sizeof(long double));
    q3EM2sagML=(long double complex*)calloc(numeroPasosq3EM2, sizeof(long double complex));
    w3EM2sagML=(long double*)calloc(numeroPasosq3EM2,sizeof(long double));

    // Propagación de q3
    q3EM2tanCW[0]=prop_q(propLibre,q2EM2tanCW[contador+1],1,1);
    w3EM2tanCW[0]=spot_q(q3EM2tanCW[0],1,lambda0);
    q3EM2sagCW[0]=prop_q(propLibre,q2EM2sagCW[contador+1],1,1);
    w3EM2sagCW[0]=spot_q(q3EM2sagCW[0],1,lambda0);
    q3EM2tanML[0]=prop_q(propLibre,q2EM2tanML[contador+1],1,1);
    w3EM2tanML[0]=spot_q(q3EM2tanML[0],1,lambda0);
    q3EM2sagML[0]=prop_q(propLibre,q2EM2sagML[contador+1],1,1);
    w3EM2sagML[0]=spot_q(q3EM2sagML[0],1,lambda0);

    for(contador=1;contador<numeroPasosq3EM2-1;contador++)
    {
        q3EM2tanCW[contador]=prop_q(propLibre,q3EM2tanCW[contador-1],1,1);
        w3EM2tanCW[contador]=spot_q(q3EM2tanCW[contador],1,lambda0);
        q3EM2sagCW[contador]=prop_q(propLibre,q3EM2sagCW[contador-1],1,1);
        w3EM2sagCW[contador]=spot_q(q3EM2sagCW[contador],1,lambda0);
        q3EM2tanML[contador]=prop_q(propLibre,q3EM2tanML[contador-1],1,1);
        w3EM2tanML[contador]=spot_q(q3EM2tanML[contador],1,lambda0);
        q3EM2sagML[contador]=prop_q(propLibre,q3EM2sagML[contador-1],1,1);
        w3EM2sagML[contador]=spot_q(q3EM2sagML[contador],1,lambda0);
    }
    contador--;
    q3EM2tanCW[contador+1]=prop_q(espCurv2Tan,q3EM2tanCW[contador],1,1);
    w3EM2tanCW[contador+1]=spot_q(q3EM2tanCW[contador+1],1,lambda0);
    q3EM2sagCW[contador+1]=prop_q(espCurv2Sag,q3EM2sagCW[contador],1,1);
    w3EM2sagCW[contador+1]=spot_q(q3EM2sagCW[contador+1],1,lambda0);
    q3EM2tanML[contador+1]=prop_q(espCurv2Tan,q3EM2tanML[contador],1,1);
    w3EM2tanML[contador+1]=spot_q(q3EM2tanML[contador+1],1,lambda0);
    q3EM2sagML[contador+1]=prop_q(espCurv2Sag,q3EM2sagML[contador],1,1);
    w3EM2sagML[contador+1]=spot_q(q3EM2sagML[contador+1],1,lambda0);

    // Creación de q4
    int numeroPasosq4EM2= L1/paso;
    q4EM2tanCW=(long double complex*)calloc(numeroPasosq4EM2, sizeof(long double complex));
    w4EM2tanCW=(long double *)calloc(numeroPasosq4EM2,sizeof(long double));
    q4EM2sagCW=(long double complex*)calloc(numeroPasosq4EM2, sizeof(long double complex));
    w4EM2sagCW=(long double *)calloc(numeroPasosq4EM2,sizeof(long double));
    q4EM2tanML=(long double complex*)calloc(numeroPasosq4EM2, sizeof(long double complex));
    w4EM2tanML=(long double *)calloc(numeroPasosq4EM2,sizeof(long double));
    q4EM2sagML=(long double complex*)calloc(numeroPasosq4EM2, sizeof(long double complex));
    w4EM2sagML=(long double *)calloc(numeroPasosq4EM2,sizeof(long double));
    // propagación de q4
    q4EM2tanCW[0]=prop_q(propLibre,q3EM2tanCW[contador+1],1,1);
    w4EM2tanCW[0]=spot_q(q4EM2tanCW[0],1,lambda0);
    q4EM2sagCW[0]=prop_q(propLibre,q3EM2sagCW[contador+1],1,1);
    w4EM2sagCW[0]=spot_q(q4EM2sagCW[0],1,lambda0);
    q4EM2tanML[0]=prop_q(propLibre,q3EM2tanML[contador+1],1,1);
    w4EM2tanML[0]=spot_q(q4EM2tanML[0],1,lambda0);
    q4EM2sagML[0]=prop_q(propLibre,q3EM2sagML[contador+1],1,1);
    w4EM2sagML[0]=spot_q(q4EM2sagML[0],1,lambda0);
    for(contador=1;contador<numeroPasosq4EM2-1;contador++)
    {
        q4EM2tanCW[contador]=prop_q(propLibre,q4EM2tanCW[contador-1],1,1);
        w4EM2tanCW[contador]=spot_q(q4EM2tanCW[contador],1,lambda0);
        q4EM2sagCW[contador]=prop_q(propLibre,q4EM2sagCW[contador-1],1,1);
        w4EM2sagCW[contador]=spot_q(q4EM2sagCW[contador],1,lambda0);
        q4EM2tanML[contador]=prop_q(propLibre,q4EM2tanML[contador-1],1,1);
        w4EM2tanML[contador]=spot_q(q4EM2tanML[contador],1,lambda0);
        q4EM2sagML[contador]=prop_q(propLibre,q4EM2sagML[contador-1],1,1);
        w4EM2sagML[contador]=spot_q(q4EM2sagML[contador],1,lambda0);
    }
    contador--;
    q4EM2tanCW[contador+1]=prop_q(propLibre,q4EM2tanCW[contador],1,1);
    w4EM2tanCW[contador+1]=spot_q(q4EM2tanCW[contador+1],1,lambda0);
    q4EM2sagCW[contador+1]=prop_q(propLibre,q4EM2sagCW[contador],1,1);
    w4EM2sagCW[contador+1]=spot_q(q4EM2sagCW[contador+1],1,lambda0);
    q4EM2tanML[contador+1]=prop_q(propLibre,q4EM2tanML[contador],1,1);
    w4EM2tanML[contador+1]=spot_q(q4EM2tanML[contador+1],1,lambda0);
    q4EM2sagML[contador+1]=prop_q(propLibre,q4EM2sagML[contador],1,1);
    w4EM2sagML[contador+1]=spot_q(q4EM2sagML[contador+1],1,lambda0);


    /// Escritura EM2tanML

    /// Estructura para guardar
    vectorDato *EM2tan=(vectorDato *)calloc(1,sizeof(vectorDato));
    EM2tan->numeroElementos=numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-4; // Ultima propagación no hay refracción.
    EM2tan->pos=(long double *)calloc(EM2tan->numeroElementos,sizeof(long double));
    EM2tan->spotCW=(long double *)calloc(EM2tan->numeroElementos,sizeof(long double));
    EM2tan->spotML=(long double *)calloc(EM2tan->numeroElementos,sizeof(long double));
    vectorDato *EM2sag=(vectorDato *)calloc(1,sizeof(vectorDato));
    EM2sag->numeroElementos=numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-4; // Ultima propagación no hay refracción.
    EM2sag->pos=(long double *)calloc(EM2sag->numeroElementos,sizeof(long double));
    EM2sag->spotCW=(long double *)calloc(EM2sag->numeroElementos,sizeof(long double));
    EM2sag->spotML=(long double *)calloc(EM2sag->numeroElementos,sizeof(long double));



    for(int i=0;i<numeroPasosq0EM2-1;i++)
    {
        EM2tan->pos[i]=(numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-5-i)*paso;
        EM2tan->spotCW[i]=w0EM2tanCW[i];
        EM2tan->spotML[i]=w0EM2tanML[i];
        EM2sag->pos[i]=(numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-5-i)*paso;
        EM2sag->spotCW[i]=w0EM2sagCW[i];
        EM2sag->spotML[i]=w0EM2sagML[i];
    }
    for(int i=0;i<numeroPasosq1EM2-1;i++)
    {
        EM2tan->pos[i+numeroPasosq0EM2-1]=(numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-4-i)*paso;
        EM2tan->spotCW[i+numeroPasosq0EM2-1]=w1EM2tanCW[i];
        EM2tan->spotML[i+numeroPasosq0EM2-1]=w1EM2tanML[i];
        EM2sag->pos[i+numeroPasosq0EM2-1]=(numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-4-i)*paso;
        EM2sag->spotCW[i+numeroPasosq0EM2-1]=w1EM2sagCW[i];
        EM2sag->spotML[i+numeroPasosq0EM2-1]=w1EM2sagML[i];
    }
    for(int i=0;i<numeroPasosq2EM2-1;i++)
    {
        if(i==0)
        {
            EM2tan->posKerrIn=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso; // Inicio de cristal
            EM2sag->posKerrIn=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso; // Inicio de cristal
        }
        if(i==(numeroPasosq2EM2-1)/2)
        {
            EM2tan->posKerrMid=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso; // Mitad de cristal
            EM2sag->posKerrMid=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso; // Mitad de cristal
        }
        if(i==(numeroPasosq2EM2-2))
        {
            EM2tan->posKerrOut=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso; // Fin de cristal
            EM2sag->posKerrOut=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso; // Fin de cristal
        }
        EM2tan->pos[i+numeroPasosq0EM2+numeroPasosq1EM2-2]=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso;
        EM2tan->spotCW[i+numeroPasosq0EM2+numeroPasosq1EM2-2]=w2EM2tanCW[i];
        EM2tan->spotML[i+numeroPasosq0EM2+numeroPasosq1EM2-2]=w2EM2tanML[i];
        EM2sag->pos[i+numeroPasosq0EM2+numeroPasosq1EM2-2]=(numeroPasosq2EM2+numeroPasosq3EM2+numeroPasosq4EM2-3-i)*paso;
        EM2sag->spotCW[i+numeroPasosq0EM2+numeroPasosq1EM2-2]=w2EM2sagCW[i];
        EM2sag->spotML[i+numeroPasosq0EM2+numeroPasosq1EM2-2]=w2EM2sagML[i];
    }
    for(int i=0;i<numeroPasosq3EM2-1;i++)
    {
        EM2tan->pos[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2-3]=(numeroPasosq3EM2+numeroPasosq4EM2-2-i)*paso;
        EM2tan->spotCW[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2-3]=w3EM2tanCW[i];
        EM2tan->spotML[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2-3]=w3EM2tanML[i];
        EM2sag->pos[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2-3]=(numeroPasosq3EM2+numeroPasosq4EM2-2-i)*paso;
        EM2sag->spotCW[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2-3]=w3EM2sagCW[i];
        EM2sag->spotML[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2-3]=w3EM2sagML[i];
    }
    for(int i=0;i<numeroPasosq4EM2-1;i++)
    {
        EM2tan->pos[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2-4]=(numeroPasosq4EM2-1-i)*paso;
        EM2tan->spotCW[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2-4]=w4EM2tanCW[i];
        EM2tan->spotML[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2-4]=w4EM2tanML[i];
        EM2sag->pos[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2-4]=(numeroPasosq4EM2-1-i)*paso;
        EM2sag->spotCW[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2-4]=w4EM2sagCW[i];
        EM2sag->spotML[i+numeroPasosq0EM2+numeroPasosq1EM2+numeroPasosq2EM2+numeroPasosq3EM2-4]=w4EM2sagML[i];
    }

    // Limpieza de vectores
    free(q0EM2tanCW);
    free(q1EM2tanCW);
    free(q2EM2tanCW);
    free(q3EM2tanCW);
    free(q4EM2tanCW);
    free(w0EM2tanCW);
    free(w1EM2tanCW);
    free(w2EM2tanCW);
    free(w3EM2tanCW);
    free(w4EM2tanCW);
    free(q0EM2tanML);
    free(q1EM2tanML);
    free(q2EM2tanML);
    free(q3EM2tanML);
    free(q4EM2tanML);
    free(w0EM2tanML);
    free(w1EM2tanML);
    free(w2EM2tanML);
    free(w3EM2tanML);
    free(w4EM2tanML);
    free(q0EM2sagCW);
    free(q1EM2sagCW);
    free(q2EM2sagCW);
    free(q3EM2sagCW);
    free(q4EM2sagCW);
    free(w0EM2sagCW);
    free(w1EM2sagCW);
    free(w2EM2sagCW);
    free(w3EM2sagCW);
    free(w4EM2sagCW);
    free(q0EM2sagML);
    free(q1EM2sagML);
    free(q2EM2sagML);
    free(q3EM2sagML);
    free(q4EM2sagML);
    free(w0EM2sagML);
    free(w1EM2sagML);
    free(w2EM2sagML);
    free(w3EM2sagML);
    free(w4EM2sagML);


    /// FIN EM2

    /// Escribe script
    gnuplotEscribe(EM1tan,EM1sag,EM2tan,EM2sag);


    /// Limpieza final

    free(EM1tan->pos);
    free(EM1tan->spotCW);
    free(EM1tan->spotML);
    free(EM1tan);
    free(EM1sag->pos);
    free(EM1sag->spotCW);
    free(EM1sag->spotML);
    free(EM1sag);
    free(EM2tan->pos);
    free(EM2tan->spotCW);
    free(EM2tan->spotML);
    free(EM2tan);
    free(EM2sag->pos);
    free(EM2sag->spotCW);
    free(EM2sag->spotML);
    free(EM2sag);

    wt1=borraMatriz(wt1);
    ws1=borraMatriz(ws1);
    wt2=borraMatriz(wt2);
    ws2=borraMatriz(ws2);
    qt1=borraMatriz(qt1);
    qs1=borraMatriz(qs1);
    qt2=borraMatriz(qt2);
    qs2=borraMatriz(qs2);


    wt1_ml=borraMatriz(wt1_ml);
    ws1_ml=borraMatriz(ws1_ml);
    wt2_ml=borraMatriz(wt2_ml);
    ws2_ml=borraMatriz(ws2_ml);
    free(qAstigEM1_ML[0]);
    qAstigEM1_ML[0]=NULL;
    free(qAstigEM1_ML[1]);
    qAstigEM1_ML[1]=NULL;
    free(qAstigEM1_ML);
    qAstigEM1_ML=NULL;
    free(qAstigEM2_ML[0]);
    qAstigEM2_ML[0]=NULL;
    free(qAstigEM2_ML[1]);
    qAstigEM2_ML[1]=NULL;
    free(qAstigEM2_ML);
    qAstigEM1_ML=NULL;
    free(ajuste->a0);
    ajuste->a0=NULL;
    free(ajuste->a1);
    ajuste->a1=NULL;
    free(ajuste->a2);
    ajuste->a2=NULL;
    free(ajuste->paso);
    ajuste->paso=NULL;
    free(ajuste);
    ajuste=NULL;
    epsilon1=borraMatriz(epsilon1);
    epsilon2=borraMatriz(epsilon2);
}

void matrizALista(matriz * mtz, matriz * epsilon1, matriz * epsilon2, char * ruta) //Sólo imprime reales
{
    FILE *archivo;
    archivo = fopen(ruta,"w");
    fprintf(archivo,"#Epsilon1 [mm] Epsilon2 [mm] Valor\n");
    for(int i=0;i<epsilon1->columnas;i++)
    {
        for(int j=0;j<epsilon2->columnas;j++)
        {
            fprintf(archivo,"%Le %Le %Le\n",obtieneElemento(epsilon1,1,i),obtieneElemento(epsilon2,1,j),creall(obtieneElemento(mtz,i,j)));
        }
    }
    fclose(archivo);
}
