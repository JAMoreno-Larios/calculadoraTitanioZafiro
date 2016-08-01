#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include "parser.h"
#include "ctocppqprogressbar.h"
#include <ccomplex>
#include "C/matrices.h"
#include <QTextStream>
#include <QTextEdit>


class worker : public QObject
{
    Q_OBJECT
public:
    explicit worker(parser *Par, CtoCppQProgressBar *barraSpots, \
                    CtoCppQProgressBar *barraTerm, CtoCppQProgressBar *barraEM1, \
                    CtoCppQProgressBar *barraEM2 , QTextEdit *salida, QObject *parent = 0);

signals:
    void finished();
    void error (QString err);

public slots:
    void process();

private:
    parser *parsePtr;
    CtoCppQProgressBar *barraSpots;
    CtoCppQProgressBar *barraTerm;
    CtoCppQProgressBar *barraEM1;
    CtoCppQProgressBar *barraEM2;
    QTextEdit *salidaTexto;
};

#endif // WORKER_H
