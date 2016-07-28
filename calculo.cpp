#include <cstdlib>
#include <cstdio>
#include <ccomplex>
#include <QMessageBox>
#include "C/lineal.h"
#include "C/matrices.h"
#include "C/no_lineal.h"
#include <cstring>
#include "C/termico.h"
#include "C/error_iteraciones.h"
#include "C/acoplamiento.h"
#include "ctocppqprogressbar.h"
#include "parser.h"
#include "C/no_linealMatAstTerm.h"
#include <QThread>
#include "calculo.h"


// Macros

#define RutaEM1CW parsePtr->mapaChar["RutaEM1CW"]
#define RutaEM2CW parsePtr->mapaChar["RutaEM2CW"]
#define RutaEM1ML parsePtr->mapaChar["RutaEM1ML"]
#define RutaEM2ML parsePtr->mapaChar["RutaEM2ML"]
#define RutaIteracionesEM1 parsePtr->mapaChar["RutaIteracionesEM1"]
#define RutaIteracionesEM2 parsePtr->mapaChar["RutaIteracionesEM2"]
#define RutaIterMax parsePtr->mapaChar["RutaIterMax"]
#define conjugado_largo parsePtr->mapaChar["conjugado_largo"]
#define conjugado_corto parsePtr->mapaChar["conjugado_corto"]
#define EM1CW parsePtr->mapaBool["EM1CW"]
#define EM2CW parsePtr->mapaBool["EM2CW"]
#define EM1ML parsePtr->mapaBool["EM1ML"]
#define EM2ML parsePtr->mapaBool["EM2ML"]
#define termico parsePtr->mapaBool["termico"]
#define guardaSpotsIter parsePtr->mapaBool["guardaSpotsIter"]
#define guardaIterMax parsePtr->mapaBool["guardaIterMax"]
#define numeroLaminas parsePtr->mapaEnteros["numeroLaminas"]
#define iteraciones parsePtr->mapaEnteros["iteraciones"]
#define NTerm parsePtr->mapaEnteros["NTerm"]
#define L1 parsePtr->mapaDobles["L1"]
#define L2 parsePtr->mapaDobles["L2"]
#define L parsePtr->mapaDobles["L"]
#define n parsePtr->mapaDobles["n"]
#define f1 parsePtr->mapaDobles["f1"]
#define f2 parsePtr->mapaDobles["f2"]
#define lambda parsePtr->mapaDobles["lambda"]
#define PLaser parsePtr->mapaDobles["PLaser"]
#define umbral parsePtr->mapaDobles["umbral"]
#define ladoX parsePtr->mapaDobles["ladoX"]
#define ladoY parsePtr->mapaDobles["ladoY"]
#define alpha parsePtr->mapaDobles["alpha"]
#define n2 parsePtr->mapaDobles["n2"]
#define nPump parsePtr->mapaDobles["nPump"]
#define kth parsePtr->mapaDobles["kth"]
#define chi parsePtr->mapaDobles["chi"]
#define Cp parsePtr->mapaDobles["Cp"]
#define rho parsePtr->mapaDobles["rho"]
#define dn_dT parsePtr->mapaDobles["dn_Dt"]
#define wPumpTan parsePtr->mapaDobles["wPumpTan"]
#define wPumpSag parsePtr->mapaDobles["wPumpSag"]
#define divPump parsePtr->mapaDobles["divPump"]
#define lambdaPump parsePtr->mapaDobles["lambdaPump"]
#define PPump parsePtr->mapaDobles["PPump"]
#define RA parsePtr->mapaDobles["RA"]
#define RB parsePtr->mapaDobles["RB"]
#define tLente parsePtr->mapaDobles["tLente"]
#define nLente parsePtr->mapaDobles["nLente"]
#define inclinacionLente parsePtr->mapaDobles["inclinacionLente"]
#define tEspejo parsePtr->mapaDobles["tEspejo"]
#define nEspejo parsePtr->mapaDobles["nEspejo"]
#define separacionFuenteALente parsePtr->mapaDobles["separacionFuenteALente"]
#define separacionLenteAEspejo parsePtr->mapaDobles["separacionLenteAEspejo"]
#define epsilon1Min parsePtr->mapaDobles["epsilon1Min"]
#define epsilon1Max parsePtr->mapaDobles["epsilon1Max"]
#define epsilon2Min parsePtr->mapaDobles["epsilon2Min"]
#define epsilon2Max parsePtr->mapaDobles["epsilon2Max"]
#define pasoEpsilon1 parsePtr->mapaDobles["pasoEpsilon1"]
#define pasoEpsilon2 parsePtr->mapaDobles["pasoEpsilon2"]
#define numPasoEpsilon1 parsePtr->mapaDobles["numPasoEpsilon1"]
#define numPasoEpsilon2 parsePtr->mapaDobles["numPasoEpsilon2"]


calculo::calculo(parser *Par, CtoCppQProgressBar *avanceNoLinealEM1,CtoCppQProgressBar *avanceNoLinealEM2,\
                 CtoCppQProgressBar *avanceTermico, CtoCppQProgressBar *avanceSpots,QObject *parent) : QObject(parent)
{
    parsePtr=Par;
    barraNoLinealEM1=avanceNoLinealEM1;
    barraNoLinealEM2=avanceNoLinealEM2;
    barraTermico=avanceTermico;
    barraSpots=avanceSpots;
}

void calculo::run()
{
    spotsPump * spots;
    matriz *wt1,*ws1,*wt2,*ws2,*qt1,*qs1,*qt2,*qs2,*epsilon1,*epsilon2;
    // Llenando epsilon1 y epsilon2
    epsilon1=nuevaMatriz(1,numPasoEpsilon1);
    epsilon2=nuevaMatriz(1,numPasoEpsilon2); // Columnas
    fijaElemento(epsilon1,1,1,epsilon1Min);
    fijaElemento(epsilon2,1,1,epsilon2Min);
    for(int index=2;index<=numPasoEpsilon1;index++)
        fijaElemento(epsilon1,1,index,obtieneElemento(epsilon1,1,index-1)+pasoEpsilon1);
    for(int index=2;index<=numPasoEpsilon2;index++)
        fijaElemento(epsilon2,1,index,obtieneElemento(epsilon2,1,index-1)+pasoEpsilon2);

    // Creando matrices receptoras
    wt1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    ws1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    wt2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    ws2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qt1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qs1=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qt2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);
    qs2=nuevaMatriz(epsilon1->columnas,epsilon2->columnas);

    // Lineal

    calculoLineal(conjugado_corto,conjugado_largo,n,lambda,L1,L2,L,f1,f2,EM1CW,EM2CW,\
                  RutaEM1CW,RutaEM2CW,wt1,ws1,wt2,ws2,epsilon1,epsilon2,qt1,qs1,qt2,qs2);
    matriz *astigEM1CW,*astigEM2CW;
    long double _Complex solEM1CW, solEM2CW;
    long double spotTanEM1CW, spotSagEM1CW, spotTanEM2CW, spotSagEM2CW;
    int e1PosEM1CW, e2PosEM1CW, e1PosEM2CW, e2PosEM2CW;
    astigEM1CW=divMatricesElemAElem(wt1,ws1);
    astigEM2CW=divMatricesElemAElem(wt2,ws2);
    confNoAstigmatica(astigEM1CW,&solEM1CW,&e1PosEM1CW,&e2PosEM1CW);
    confNoAstigmatica(astigEM2CW,&solEM2CW,&e1PosEM2CW,&e2PosEM2CW);
    spotTanEM1CW=creall(obtieneElemento(wt1,e1PosEM1CW,e2PosEM1CW));
    spotSagEM1CW=creall(obtieneElemento(ws1,e1PosEM1CW,e2PosEM1CW));
    spotTanEM2CW=creall(obtieneElemento(wt2,e1PosEM2CW,e2PosEM2CW));
    spotSagEM2CW=creall(obtieneElemento(ws2,e1PosEM2CW,e2PosEM2CW));
    long double deltaFin1EM1, deltaFin2EM1, deltaFin1EM2, deltaFin2EM2;
    long double angulos[2];
    anguloLineal(conjugado_corto,conjugado_largo,L,n,L1,L2,f1,f2,angulos);
    deltaFin1EM1=distanciaCristal(conjugado_corto,f1,L1,L,n,angulos[0])+creall(obtieneElemento(epsilon1,1,e1PosEM1CW));
    deltaFin2EM1=distanciaCristal(conjugado_largo,f2,L2,L,n,angulos[1])+creall(obtieneElemento(epsilon2,1,e2PosEM1CW));
    //salida.append(QString("Para EM1 CW:\n\nTheta1= %1 [°]\nTheta2= %2 [°]\nDelta1= %3 [mm]\nDelta2= %4 [mm]\n").arg((double)angulos[0]*180.0/M_PI).arg((double)angulos[1]*180.0/M_PI).arg((double)deltaFin1EM1*1e3).arg((double)deltaFin2EM1*1e3));
    //salida.append(QString("\nSpot EM1 tangencial= %1 [um]\nSpot EM1 Sagital= %1 [um]\n").arg((double)spotTanEM1CW*1e6).arg((double)spotSagEM1CW*1e6));
    deltaFin1EM2=distanciaCristal(conjugado_corto,f1,L1,L,n,angulos[0])+creall(obtieneElemento(epsilon1,1,e1PosEM2CW));
    deltaFin2EM2=distanciaCristal(conjugado_largo,f2,L2,L,n,angulos[1])+creall(obtieneElemento(epsilon2,1,e2PosEM2CW));
    //salida.append(QString("Para EM2 CW:\n\nTheta1= %1 [°]\nTheta2= %2 [°]\nDelta1= %3 [mm]\nDelta2= %4 [mm]\n").arg((double)angulos[0]*180.0/M_PI).arg((double)angulos[1]*180.0/M_PI).arg((double)deltaFin1EM2*1e3).arg((double)deltaFin2EM2*1e3));
    //salida.append(QString("\nSpot EM2 tangencial= %1 [um]\nSpot EM2 Sagital= %1 [um]\n").arg((double)spotTanEM2CW*1e6).arg((double)spotSagEM2CW*1e6));
    //ui->salidaTexto->setText(salida);

    // No lineal

    if(this->EM1ML||this->EM2ML)
    {
        spots=acopleOpticoBarra(this->wPumpTan,this->wPumpSag,this->divPump,this->lambdaPump,this->separacionFuenteALente,this->tLente,this->RA,this->RB,this->nLente,\
                           this->inclinacionLente,this->separacionLenteAEspejo,this->tEspejo,INFINITY,f1*2,this->nEspejo,angulos[0],n,deltaFin1EM1,epsilon1,barraSpots); // Verificar si delta1 difiere

        ajusteTemperaturaCristal **matrizVectores = matrizAjusteCristal(spots,epsilon1,alpha,lambdaPump,n,PPump,chi,L,kth,Cp,rho,dn_dT,ladoX,ladoY,iteraciones,\
                                                                        numeroLaminas,umbral,1.3,NTerm,termico,barraTermico);
        // Acopladas termicas
        matriz **qAstigEM1_MLterm;
        matriz **qAstigEM2_MLterm;
        matriz *wt1MLmatAcopTerm, *ws1MLmatAcopTerm, *wt2MLmatAcopTerm, *ws2MLmatAcopTerm, *astigMatAcopTermEM1, *astigMatAcopTermEM2;

        if(this->EM1ML)
        {
            qAstigEM1_MLterm=propNoLinealEM1Astigmatico(qt1,qs1,matrizVectores,conjugado_corto,conjugado_largo,L1,L2,L,f1,f2,n,n2,chi,kth,Cp,rho,dn_dT,PLaser,\
                                                        lambda,epsilon1,epsilon2,iteraciones,numeroLaminas,umbral,guardaIterMax,guardaSpotsIter,\
                                                        RutaIterMax,RutaIteracionesEM1,RutaIteracionesEM1,barraNoLinealEM1);
            wt1MLmatAcopTerm=spot_q_matriz(qAstigEM1_MLterm[0],1,lambda);
            ws1MLmatAcopTerm=spot_q_matriz(qAstigEM1_MLterm[1],1,lambda);
            astigMatAcopTermEM1=divMatricesElemAElem(wt1MLmatAcopTerm,ws1MLmatAcopTerm);
            imprimeListaReal(epsilon1,epsilon2,wt1MLmatAcopTerm,ws1MLmatAcopTerm,RutaEM1ML);
            double delta1=distanciaCristal(conjugado_corto,f1,L1,L,n,angulos[0]);
            double delta2=distanciaCristal(conjugado_largo,f2,L2,L,n,angulos[1]);
            int fila, columna;
            long double _Complex sol;
            confNoAstigmatica(astigMatAcopTermEM1,&sol,&fila,&columna);
            //escribeMatrizArchivo_completo(astigMatAcopTermEM1,epsilon1,epsilon2,angulos,delta1,delta2,conjugado_corto,conjugado_largo,RutaEM1ML,fila,columna,creall(sol));
            imprimeListaReal(epsilon1,epsilon2,wt1MLmatAcopTerm,ws1MLmatAcopTerm,RutaEM1ML);
        }
        if(this->EM2ML)
        {
            qAstigEM2_MLterm=propNoLinealEM2Astigmatico(qt2,qs2,matrizVectores,conjugado_corto,conjugado_largo,L1,L2,L,f1,f2,n,n2,chi,kth,Cp,rho,dn_dT,PLaser,\
                                                        lambda,epsilon1,epsilon2,iteraciones,numeroLaminas,umbral,guardaIterMax,guardaSpotsIter,\
                                                        RutaIterMax,RutaIteracionesEM2,RutaIteracionesEM2,barraNoLinealEM2);
            wt2MLmatAcopTerm=spot_q_matriz(qAstigEM2_MLterm[0],1,lambda);
            ws2MLmatAcopTerm=spot_q_matriz(qAstigEM2_MLterm[1],1,lambda);
            astigMatAcopTermEM2=divMatricesElemAElem(wt2MLmatAcopTerm,ws2MLmatAcopTerm);
            imprimeListaReal(epsilon1,epsilon2,wt2MLmatAcopTerm,ws2MLmatAcopTerm,RutaEM2ML);
            double delta1=distanciaCristal(conjugado_corto,f1,L1,L,n,angulos[0]);
            double delta2=distanciaCristal(conjugado_largo,f2,L2,L,n,angulos[1]);
            int fila, columna;
            long double _Complex sol;
            confNoAstigmatica(astigMatAcopTermEM2,&sol,&fila,&columna);
            //escribeMatrizArchivo_completo(astigMatAcopTermEM2,epsilon1,epsilon2,angulos,delta1,delta2,conjugado_corto,conjugado_largo,RutaEM2ML,fila,columna,creall(sol));
            imprimeListaReal(epsilon1,epsilon2,wt2MLmatAcopTerm,ws2MLmatAcopTerm,RutaEM2ML);
        }


        printf("No lineal matrices acopladas con lente térmica terminado\n");

        wt1MLmatAcopTerm=borraMatriz(wt1MLmatAcopTerm);
        ws1MLmatAcopTerm=borraMatriz(ws1MLmatAcopTerm);
        wt2MLmatAcopTerm=borraMatriz(wt2MLmatAcopTerm);
        ws2MLmatAcopTerm=borraMatriz(ws2MLmatAcopTerm);
        astigMatAcopTermEM1=borraMatriz(astigMatAcopTermEM1);
        astigMatAcopTermEM2=borraMatriz(astigMatAcopTermEM2);

        // Borra térmico

        matrizVectores=borraMatrizAjusteCristal(matrizVectores,spots->numSpots);
        spots=borraSpotsPump(spots);
    }


    // Limpieza de matrices
    parsePtr->escribeDatos(parsePtr->mapaChar,parsePtr->mapaBool,parsePtr->mapaEnteros,parsePtr->mapaDobles);

    wt1=borraMatriz(wt1);
    ws1=borraMatriz(ws1);
    wt2=borraMatriz(wt2);
    ws2=borraMatriz(ws2);
    qt1=borraMatriz(qt1);
    qs1=borraMatriz(qs1);
    qt2=borraMatriz(qt2);
    qs2=borraMatriz(qs2);
    astigEM1CW=borraMatriz(astigEM1CW);
    astigEM2CW=borraMatriz(astigEM2CW);


    // Limpieza de rutas de texto
    if(RutaEM1CW!=NULL) {free(RutaEM1CW); RutaEM1CW=NULL;}
    if(RutaEM2CW!=NULL) {free(RutaEM2CW); RutaEM2CW=NULL;}
    if(RutaEM1ML!=NULL) {free(RutaEM1ML); RutaEM1ML=NULL;}
    if(RutaEM2ML!=NULL) {free(RutaEM2ML); RutaEM2ML=NULL;}
    if(RutaIterMax!=NULL) {free(RutaIterMax); RutaIterMax=NULL;}
    if(RutaIteracionesEM1!=NULL) {free(RutaIteracionesEM1); RutaIteracionesEM1=NULL;}
    if(RutaIteracionesEM2!=NULL) {free(RutaIteracionesEM2); RutaIteracionesEM2=NULL;}
    parsePtr->cargaVariables();

    //emit finished;

}
