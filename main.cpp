#include "interfaz1.h"
#include <QApplication>
#include "parser.h"
#include <iostream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    std::string ruta="config.data";
    parser *Par = new parser(ruta);
    interfaz1 w(Par);
    w.setDefaults();
    w.show();

    return a.exec();
}
