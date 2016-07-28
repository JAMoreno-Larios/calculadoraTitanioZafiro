#ifndef INTERFAZ1_H
#define INTERFAZ1_H

#include <QMainWindow>
#include <string>
#include "opcionesLaser.h"
#include "parser.h"
#include "calculo.h"

using std::string;

namespace Ui {
class interfaz1;
}

class interfaz1 : public QMainWindow
{
    Q_OBJECT

public:
    explicit interfaz1(parser *Par, QWidget *parent = 0);
    ~interfaz1();

    void setDefaults();
    //opcionesLaser *A;

    void setRutas(string nuevaRutaEM1CW, string nuevaRutaEM2CW, string nuevaRutaEM1ML, string nuevaRutaEM2ML,\
                             string nuevaRutaiteracionesEM1, string nuevaRutaiteracionesEM2, string nuevaRutaiterMax, parser *parsePtr);
    void setCasos(bool EM1MLflag, bool EM2MLflag, bool termicoflag, bool guardaSpotIteraciones, bool guardaIteracionesMax, parser * parsePtr);
    void setVariables(int numeroLaminasNuevo,int iteracionesNuevo,double umbralNuevo,int NTermNuevo,double ladoXNuevo,double ladoYNuevo,double alphaNuevo, double wPumpTanNuevo,\
                                 double wPumpSagNuevo,double divPumpNuevo,double lambdaPumpNuevo,double PPumpNuevo,double RANuevo,double RBNuevo,double tLenteNuevo,double nLenteNuevo,\
                                 double inclinacionLenteNuevo,double tEspejoNuevo,double nEspejoNuevo,double separacionFuenteALenteNuevo,double separacionLenteAEspejoNuevo,\
                                 double PLaserNuevo, parser * parsePtr);
    void setConjugados(char *corto, char *largo);


private slots:
    void on_initButton_clicked();

    void on_actionImporta_constantes_triggered();

    void on_actionSalir_triggered();

    void on_actionOpcionesNoLineal_triggered();

    void on_actionAcerca_de_triggered();

private:
    parser *parsePtr;
    Ui::interfaz1 *ui;
    calculo *noLineal;
    bool verificaEpsilon(double epsilon1Min, double epsilon1Max, double epsilon1Step, double epsilon2Min, double epsilon2Max, double epsilon2Step);
    void controlToggle(bool desactiva);
    char* RutaEM1CW=NULL;
    char* RutaEM2CW=NULL;
    char* RutaEM1ML=NULL;
    char* RutaEM2ML=NULL;
    char* RutaIteraciones=NULL;
    char* RutaIterMax=NULL;
    bool EM1CW;
    bool EM2CW;
    bool EM1ML;
    bool EM2ML;
    bool termico;
    bool guardaSpotsIter;
    bool guardaIterMax;
    // Variables no lineales
    int numeroLaminas;
    int iteraciones;
    double umbral;
    double PLaser;
    // Variables de cálculo térmico
    int NTerm;
    double ladoX;
    double ladoY;
    double alpha;
    double n2;
    double nPump;
    double kth;
    double chi;
    double Cp;
    double rho;
    double dn_dT;

    // Variables fuente láser
    double wPumpTan;
    double wPumpSag;
    double divPump;
    double lambdaPump;
    double PPump;
    // Variables lente
    double RA;
    double RB;
    double tLente;
    double nLente;
    double inclinacionLente;
    // Variables espejo
    double tEspejo;
    double nEspejo;
    // Variables separaciones
    double separacionFuenteALente;
    double separacionLenteAEspejo;

};

#endif // INTERFAZ1_H
