#include "salidatexto.h"
#include "ui_salidatexto.h"
#include <QString>
#include <ctime>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>

salidaTexto::salidaTexto(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::salidaTexto)
{
    ui->setupUi(this);
}

salidaTexto::~salidaTexto()
{
    delete ui;
}

// Escribe a control de texto. La cadena debe contener saltos de línea.
void salidaTexto::escribeLineaTexto(QString texto)
{
    QString cadenaSalida;
    time(&rawtime);
    timeinfo=localtime(&rawtime);
    strftime(fecha,30,"%D - %T:\n~~~\n",timeinfo);
    cadenaSalida.append(fecha);
    cadenaSalida.append(texto);
    cadenaSalida.append("~~~\n\n");
    ui->texto->appendPlainText(cadenaSalida);
}

void salidaTexto::on_actionGuardar_triggered()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::AnyFile);
    menu.setNameFilter("Texto (*.txt");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    QFile archivo(nombreArchivo[0]);
    if (!archivo.open(QIODevice::WriteOnly|QIODevice::Text))
        return;
    QTextStream salidaArchivo(&archivo);
    salidaArchivo<<ui->texto->toPlainText();
    archivo.close();
}

void salidaTexto::on_actionCerrar_triggered()
{
    this->close();
}
