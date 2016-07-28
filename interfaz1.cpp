#include "interfaz1.h"
#include "ui_interfaz1.h"
#include <cstdlib>
#include <cstdio>
#include <ccomplex>
#include <QMessageBox>
#include "C/lineal.h"
#include "C/matrices.h"
#include <cmath>
#include <QFile>
#include <QFileDialog>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include "C/no_lineal.h"
#include "opcionesLaser.h"
#include <cstring>
#include "C/termico.h"
#include "C/error_iteraciones.h"
#include "C/acoplamiento.h"
#include "ctocppqprogressbar.h"
#include "parser.h"
//#include "macros.h"
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

// Programa

interfaz1::interfaz1(parser * Par, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::interfaz1)
{
    parsePtr=Par;
    ui->setupUi(this);
}

interfaz1::~interfaz1()
{
    if(strlen(RutaEM1CW)!=0) {free(RutaEM1CW); RutaEM1CW=NULL;}
    if(strlen(RutaEM2CW)!=0) {free(RutaEM2CW); RutaEM2CW=NULL;}
    if(strlen(RutaEM1ML)!=0) {free(RutaEM1ML); RutaEM1ML=NULL;}
    if(strlen(RutaEM2ML)!=0) {free(RutaEM2ML); RutaEM2ML=NULL;}
    if(strlen(RutaIteracionesEM1)!=0) {free(RutaIteracionesEM1); RutaIteracionesEM1=NULL;}
    if(strlen(RutaIteracionesEM2)!=0) {free(RutaIteracionesEM2); RutaIteracionesEM2=NULL;}
    if(strlen(RutaIterMax)!=0) {free(RutaIterMax); RutaIterMax=NULL;}
    delete ui;
}

void interfaz1::on_initButton_clicked()
{
    this->controlToggle(true); // Desactiva controles
    // Relee el archivo de configuración
    if(parsePtr->archivoVacio())
        parsePtr->escribeDatos(parsePtr->mapaChar,parsePtr->mapaBool,parsePtr->mapaEnteros,parsePtr->mapaDobles);
    parsePtr->cargaVariables();
    // Prepara caja de mensaje
    QMessageBox aviso;
    aviso.setWindowTitle("Aviso");
    aviso.addButton(QMessageBox::Ok);

    // Forza las barras de progreso a cero
    ui->avanceNoLinealEM1->setValue(0);
    ui->avanceNoLinealEM2->setValue(0);
    ui->avanceTermico->setValue(0);
    ui->avanceSpots->setValue(0);

    // Guarda valores
    L1=ui->L1SpinBox->value();
    L2=ui->L2SpinBox->value();
    L=ui->LSpinBox->value()*1e-3;
    n=ui->nSpinBox->value();
    f1=ui->f1SpinBox->value()*1e-3;
    f2=ui->f2SpinBox->value()*1e-3;
    lambda=ui->lambdaLSpinBox->value()*1e-9;
    epsilon1Min=(long double)ui->epsilon1_min->value()*1e-3;
    epsilon1Max=(long double)ui->epsilon1_max->value()*1e-3;
    epsilon2Min=(long double)ui->epsilon2_min->value()*1e-3;
    epsilon2Max=(long double)ui->epsilon2_max->value()*1e-3;
    pasoEpsilon1=(long double)ui->pasoEpsilon1Box->value()*1e-3;
    pasoEpsilon2=(long double)ui->pasoEpsilon2Box->value()*1e-3;
    bool errorEpsilon=interfaz1::verificaEpsilon((double)epsilon1Min,(double)epsilon1Max,(double)pasoEpsilon1,(double)epsilon2Min,(double)epsilon2Max,(double)pasoEpsilon2);
    if(errorEpsilon) return;
    numPasoEpsilon1=(epsilon1Max-epsilon1Min)/pasoEpsilon1+1+1; // Compensado por error de truncamiento y considerando la existencia del cero
    numPasoEpsilon2=(epsilon2Max-epsilon2Min)/pasoEpsilon2+1+1;
    char conjugado_largoLocal[4];
    char conjugado_cortoLocal[4];

    int flagL1=ui->conjugadoL1->checkedId(); // Primer botón es -2, segundo es -3.
    int flagL2=ui->conjugadoL2->checkedId();
    //spotsPump * spots;
    QString salida=QString("Configuraciones menos astigmáticas\n");
    switch(flagL1)
    {
        case -2:
            strcpy(conjugado_cortoLocal,"fin");
        break;
        case -3:
            strcpy(conjugado_cortoLocal,"inf");
        break;
        default:
            aviso.setText("Seleccione conjugado para L1");
            aviso.exec();
            return;
    }
    switch(flagL2)
    {
        case -2:
            strcpy(conjugado_largoLocal,"fin");
        break;
        case -3:
            strcpy(conjugado_largoLocal,"inf");
        break;
        default:
        aviso.setText("Seleccione conjugado para L2");
        aviso.exec();
        return;
    }
    setConjugados(conjugado_cortoLocal,conjugado_largoLocal);
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
                           this->inclinacionLente,this->separacionLenteAEspejo,this->tEspejo,INFINITY,f1*2,this->nEspejo,angulos[0],n,deltaFin1EM1,epsilon1,ui->avanceSpots); // Verificar si delta1 difiere

        ajusteTemperaturaCristal **matrizVectores = matrizAjusteCristal(spots,epsilon1,alpha,lambdaPump,n,PPump,chi,L,kth,Cp,rho,dn_dT,ladoX,ladoY,iteraciones,\
                                                                        numeroLaminas,umbral,1.3,NTerm,termico,ui->avanceTermico);
        // Acopladas termicas
        matriz **qAstigEM1_MLterm;
        matriz **qAstigEM2_MLterm;
        matriz *wt1MLmatAcopTerm, *ws1MLmatAcopTerm, *wt2MLmatAcopTerm, *ws2MLmatAcopTerm, *astigMatAcopTermEM1, *astigMatAcopTermEM2;

        if(this->EM1ML)
        {
            qAstigEM1_MLterm=propNoLinealEM1Astigmatico(qt1,qs1,matrizVectores,conjugado_corto,conjugado_largo,L1,L2,L,f1,f2,n,n2,chi,kth,Cp,rho,dn_dT,PLaser,\
                                                        lambda,epsilon1,epsilon2,iteraciones,numeroLaminas,umbral,guardaIterMax,guardaSpotsIter,\
                                                        RutaIterMax,RutaIteracionesEM1,RutaIteracionesEM1,ui->avanceNoLinealEM1);
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
                                                        RutaIterMax,RutaIteracionesEM2,RutaIteracionesEM2,ui->avanceNoLinealEM2);
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

    this->controlToggle(false); // Reactiva controles
}

void interfaz1::on_actionImporta_constantes_triggered()
{
    QStringList variables;
    variables<<"L1"<<"L2"<<"L"<<"n0"<<"n2"<<"nPump"<<"kth"<<"chi"<<"Cp"<<"rho"<<"P_pump"<<"P_laser"<<"dn_dv"<<"lambda0"<<"lambdaPump"<<"f1"<<"f2";
    //archivo=QFileDialog::getOpenFileName(this,"Importa constantes",".","Texto plano (*.txt)");
    QFileDialog dialogo(this);
    dialogo.setFileMode(QFileDialog::ExistingFile);
    dialogo.setNameFilter("Texto plano (*.txt)");
    dialogo.setViewMode(QFileDialog::List);
    QStringList nombreArchivo;
    if(dialogo.exec())
        nombreArchivo=dialogo.selectedFiles();
    QFile archivo(nombreArchivo[0]);
    if (!archivo.open(QIODevice::ReadOnly|QIODevice::Text))
        return;
    QTextStream texto(&archivo);
    while(!texto.atEnd())
    {
        QStringList linea=texto.readLine().split(' ');

        for(int i=0; i<variables.size(); i++)
        {
            if(QString::compare(linea[0],variables[i],Qt::CaseInsensitive)==0)
                switch(i)
                {
                    case 0:
                        ui->L1SpinBox->setValue(linea.at(1).toDouble());
                        L1=linea.at(1).toDouble();
                    break;
                    case 1:
                        ui->L2SpinBox->setValue(linea.at(1).toDouble());
                        L2=linea.at(1).toDouble();
                    break;
                    case 2:
                        ui->LSpinBox->setValue(linea.at(1).toDouble()*1e3);
                        L=linea.at(1).toDouble();
                    break;
                    case 3:
                        ui->nSpinBox->setValue(linea.at(1).toDouble());
                        n=linea.at(1).toDouble();
                    break;
                    case 4:
                        n2=linea.at(1).toDouble();
                    break;
                    case 5:
                        nPump=linea.at(1).toDouble();
                    break;
                    case 6:
                        kth=linea.at(1).toDouble();
                    break;
                    case 7:
                        chi=linea.at(1).toDouble();
                    break;
                    case 8:
                        Cp=linea.at(1).toDouble();
                    break;
                    case 9:
                        rho=linea.at(1).toDouble();
                    break;
                    case 10:
                        PPump=linea.at(1).toDouble();
                    break;
                    case 11:
                        PLaser=linea.at(1).toDouble();
                    break;
                    case 12:
                        dn_dT=linea.at(1).toDouble();
                    break;
                    case 13:
                       ui->lambdaLSpinBox->setValue(linea.at(1).toDouble()*1e9);
                       lambda=linea.at(1).toDouble();
                    break;
                    case 14:
                        lambdaPump=linea.at(1).toDouble();
                    break;
                    case 15:
                        ui->f1SpinBox->setValue(linea.at(1).toDouble()*1e3);
                        f1=linea.at(1).toDouble();
                    break;
                    case 16:
                        ui->f2SpinBox->setValue(linea.at(1).toDouble()*1e3);
                        f2=linea.at(1).toDouble();
                    break;
                }
        }
    }
    parsePtr->escribeDatos(parsePtr->mapaChar,parsePtr->mapaBool,parsePtr->mapaEnteros,parsePtr->mapaDobles);
}

void interfaz1::on_actionSalir_triggered()
{
   close();
}

// Verifica que epsilonMin<epsilonMax
bool interfaz1::verificaEpsilon(double Epsilon1Min, double Epsilon1Max, double Epsilon1Step, double Epsilon2Min, double Epsilon2Max, double Epsilon2Step)
{
    QMessageBox aviso;
    aviso.setWindowTitle("Aviso");
    aviso.addButton(QMessageBox::Ok);
    bool hayError=0;
    QString error=QString("Se detectaron los siguientes errores al definir epsilon1 y epsilon2\n\n");
    if(Epsilon1Min>Epsilon1Max){
        error.append("El valor mínimo de epsilon1 es menor que el máximo\n");
        hayError=1;}
    if(Epsilon2Min>Epsilon2Max){
        error.append("El valor mínimo de epsilon2 es menor que el máximo\n");
        hayError=1;}
    if(Epsilon1Step==0.0){
        error.append("El paso para epsilon1 es cero\n.");
        hayError=1;}
    if(Epsilon2Step==0.0){
        error.append("El paso para epsilon2 es cero\n.");
        hayError=1;}
    if(hayError)
    {
        aviso.setText(error);
        aviso.exec();
    }
    return hayError;
}

void interfaz1::on_actionOpcionesNoLineal_triggered()
{
    opcionesLaser w(parsePtr);
    w.A=this;
    w.exec();
}

void interfaz1::setRutas(string nuevaRutaEM1CW, string nuevaRutaEM2CW, string nuevaRutaEM1ML, string nuevaRutaEM2ML,\
                         string nuevaRutaiteracionesEM1, string nuevaRutaiteracionesEM2, string nuevaRutaiterMax, parser *parsePtr)
{

    int caracteres=0;

            caracteres=snprintf(NULL,0,nuevaRutaEM1CW.c_str());
            RutaEM1CW=(char*)malloc(caracteres+1);
            caracteres=snprintf(NULL,0,nuevaRutaEM2CW.c_str());
            RutaEM2CW=(char*)malloc(caracteres+1);
            caracteres=snprintf(NULL,0,nuevaRutaEM1ML.c_str());
            RutaEM1ML=(char*)malloc(caracteres+1);
            caracteres=snprintf(NULL,0,nuevaRutaEM2ML.c_str());
            RutaEM2ML=(char*)malloc(caracteres+1);
            caracteres=snprintf(NULL,0,nuevaRutaiteracionesEM1.c_str());
            RutaIteracionesEM1=(char*)malloc(caracteres+1);
            caracteres=snprintf(NULL,0,nuevaRutaiteracionesEM2.c_str());
            RutaIteracionesEM2=(char*)malloc(caracteres+1);
            caracteres=snprintf(NULL,0,nuevaRutaiterMax.c_str());
            RutaIterMax=(char*)malloc(caracteres+1);

            sprintf(RutaEM1CW,nuevaRutaEM1CW.c_str());
            sprintf(RutaEM2CW,nuevaRutaEM2CW.c_str());
            sprintf(RutaEM1ML,nuevaRutaEM1ML.c_str());
            sprintf(RutaEM2ML,nuevaRutaEM2ML.c_str());
            sprintf(RutaIteracionesEM1,nuevaRutaiteracionesEM1.c_str());
            sprintf(RutaIteracionesEM2,nuevaRutaiteracionesEM2.c_str());
            sprintf(RutaIterMax,nuevaRutaiterMax.c_str());
}

void interfaz1::setCasos(bool EM1MLflag, bool EM2MLflag, bool termicoflag, bool guardaSpotIteraciones, bool guardaIteracionesMax, parser * parsePtr)
{
    EM1ML=EM1MLflag;
    EM2ML=EM2MLflag;
    termico=termicoflag;
    guardaSpotsIter=guardaSpotIteraciones;
    guardaIterMax=guardaIteracionesMax;
}

void interfaz1::setVariables(int numeroLaminasNuevo,int iteracionesNuevo,double umbralNuevo,int NTermNuevo,double ladoXNuevo,double ladoYNuevo,double alphaNuevo, double wPumpTanNuevo,\
                             double wPumpSagNuevo,double divPumpNuevo,double lambdaPumpNuevo,double PPumpNuevo,double RANuevo,double RBNuevo,double tLenteNuevo,double nLenteNuevo,\
                             double inclinacionLenteNuevo,double tEspejoNuevo,double nEspejoNuevo,double separacionFuenteALenteNuevo,double separacionLenteAEspejoNuevo,\
                             double PLaserNuevo, parser * parsePtr)
{
    // Asignar a puntero
    numeroLaminas=numeroLaminasNuevo;
    iteraciones=iteracionesNuevo;
    umbral=umbralNuevo;
    PLaser=PLaserNuevo*1e3;
    // Variables de cálculo térmico
    NTerm=NTermNuevo;
    ladoX=ladoXNuevo*1e-3;
    ladoY=ladoYNuevo*1e-3;
    alpha=alphaNuevo;
    // Variables fuente láser
    wPumpTan=wPumpTanNuevo*1e-3;
    wPumpSag = wPumpSagNuevo*1e-3;
    divPump = divPumpNuevo*1e-3;
    lambdaPump = lambdaPumpNuevo*1e-9;
    PPump = PPumpNuevo;
    // Variables lente
    RA = RANuevo*1e-3;
    RB = RBNuevo*1e-3;
    tLente = tLenteNuevo*1e-3;
    nLente = nLenteNuevo;
    inclinacionLente = inclinacionLenteNuevo;
    // Variables espejo
    tEspejo = tEspejoNuevo*1e-3;
    nEspejo = nEspejoNuevo;
    // Variables separaciones
    separacionFuenteALente = separacionFuenteALenteNuevo*1e-3;
    separacionLenteAEspejo = separacionLenteAEspejoNuevo*1e-3;
}

void interfaz1::setDefaults()
{
    int caracteres=0;
    caracteres=snprintf(NULL,0,"./EM1CW.txt");
    parsePtr->mapaChar["RutaEM1CW"]=this->RutaEM1CW=(char*)malloc(caracteres+1);
    caracteres=snprintf(NULL,0,"./EM2CW.txt");
    parsePtr->mapaChar["RutaEM2CW"]=this->RutaEM2CW=(char*)malloc(caracteres+1);
    caracteres=snprintf(NULL,0,"./EM1ML.txt");
    parsePtr->mapaChar["RutaEM1ML"]=this->RutaEM1ML=(char*)malloc(caracteres+1);
    caracteres=snprintf(NULL,0,"./EM2ML.txt");
    this->RutaEM2ML=(char*)malloc(caracteres+1);
    caracteres=snprintf(NULL,0,".");
    this->RutaIteraciones=(char*)malloc(caracteres+1);
    sprintf(this->RutaEM1CW,"./EM1CW.txt");
    sprintf(this->RutaEM2CW,"./EM2CW.txt");
    sprintf(this->RutaEM1ML,"./EM1ML.txt");
    sprintf(this->RutaEM2ML,"./EM2ML.txt");
    sprintf(this->RutaIteraciones,".");

    this->EM1CW=true;
    this->EM2CW=true;
    this->EM1ML=false;
    this->EM2ML=false;
    this->termico=false;
}

void interfaz1::on_actionAcerca_de_triggered()
{
    QMessageBox acercaDe;
    acercaDe.about(this,"Calculadora Ti:Zafiro","Escrito por José Agustín Moreno Larios\n2016");
}

void interfaz1::controlToggle(bool desactiva)
{
    if(desactiva)
    {
        ui->menubar->setEnabled(false);
        // Botones de radio
        ui->EM1Fin->setEnabled(false);
        ui->EM1Inf->setEnabled(false);
        ui->EM2Fin->setEnabled(false);
        ui->EM2Inf->setEnabled(false);
        // Controles
        ui->L1SpinBox->setEnabled(false);
        ui->L2SpinBox->setEnabled(false);
        ui->LSpinBox->setEnabled(false);
        ui->f1SpinBox->setEnabled(false);
        ui->f2SpinBox->setEnabled(false);
        ui->lambdaLSpinBox->setEnabled(false);
        ui->nSpinBox->setEnabled(false);
        // epsilon
        ui->epsilon1_min->setEnabled(false);
        ui->epsilon1_max->setEnabled(false);
        ui->epsilon2_min->setEnabled(false);
        ui->epsilon2_max->setEnabled(false);
        ui->pasoEpsilon1Box->setEnabled(false);
        ui->pasoEpsilon2Box->setEnabled(false);
        // Boton Init
        ui->initButton->setEnabled(false);
        ui->initButton->setText("Calculando");
    }
    else
    {
        ui->menubar->setEnabled(true);
        // Botones de radio
        ui->EM1Fin->setEnabled(true);
        ui->EM1Inf->setEnabled(true);
        ui->EM2Fin->setEnabled(true);
        ui->EM2Inf->setEnabled(true);
        // Controles
        ui->L1SpinBox->setEnabled(true);
        ui->L2SpinBox->setEnabled(true);
        ui->LSpinBox->setEnabled(true);
        ui->f1SpinBox->setEnabled(true);
        ui->f2SpinBox->setEnabled(true);
        ui->lambdaLSpinBox->setEnabled(true);
        ui->nSpinBox->setEnabled(true);
        // epsilon
        ui->epsilon1_min->setEnabled(true);
        ui->epsilon1_max->setEnabled(true);
        ui->epsilon2_min->setEnabled(true);
        ui->epsilon2_max->setEnabled(true);
        ui->pasoEpsilon1Box->setEnabled(true);
        ui->pasoEpsilon2Box->setEnabled(true);
        // Boton Init
        ui->initButton->setEnabled(true);
        ui->initButton->setText("Corre");
    }
}

void interfaz1::setConjugados(char *corto, char *largo)
{
    parsePtr->mapaChar["conjugado_corto"]=(char*)calloc(4,sizeof(char));
    parsePtr->mapaChar["conjugado_largo"]=(char*)calloc(4,sizeof(char));
    strcpy(parsePtr->mapaChar["conjugado_corto"],corto);
    strcpy(parsePtr->mapaChar["conjugado_largo"],largo);
}
