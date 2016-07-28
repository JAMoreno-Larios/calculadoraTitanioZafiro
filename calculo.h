#ifndef CALCULO_H
#define CALCULO_H

#include <QThread>
#include "parser.h"
#include <ccomplex>
#include <cstdlib>
#include <cstdio>
#include "ctocppqprogressbar.h"


class calculo : public QThread
{
    Q_OBJECT
public:
    explicit calculo(parser *Par, CtoCppQProgressBar *avanceNoLinealEM1,CtoCppQProgressBar *avanceNoLinealEM2,\
                                     CtoCppQProgressBar *avanceTermico, CtoCppQProgressBar *avanceSpots,QObject *parent=0);
    void fijaPtrBarras(CtoCppQProgressBar *avanceNoLinealEM1,CtoCppQProgressBar *avanceNoLinealEM2,CtoCppQProgressBar *avanceTermico, CtoCppQProgressBar *avanceSpots);

protected:
    void run();

public slots:
   // void ejecutaCalculo();

private:
    parser *parsePtr;
    CtoCppQProgressBar *barraNoLinealEM1;
    CtoCppQProgressBar *barraNoLinealEM2;
    CtoCppQProgressBar *barraTermico;
    CtoCppQProgressBar *barraSpots;
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

#endif // CALCULO_H
