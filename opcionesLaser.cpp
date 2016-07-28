#include "opcionesLaser.h"
#include "ui_opcionesLaser.h"
#include <QFile>
#include <QFileDialog>
#include <QString>
#include <QStringList>
#include <QTextStream>
//#include "interfaz1.h"
#include <string>
#include <math.h>
#include <QMessageBox>
#include "parser.h"

using std::string;

opcionesLaser::opcionesLaser(parser *Par, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::opcionesLaser)
{
    parsePtr=Par;
    ui->setupUi(this);
    agrupaBotonesNoLineal();
    agrupaBotonesRadios();
    configuraControles();
}

opcionesLaser::~opcionesLaser()
{
    delete ui;
}

void opcionesLaser::agrupaBotonesNoLineal()
{
    QButtonGroup* grupoNoLineal = new QButtonGroup(this);
    grupoNoLineal->addButton(ui->enableEM1ML,1);
    grupoNoLineal->addButton(ui->enableEM2ML,2);
    grupoNoLineal->addButton(ui->enableTermico,3);
    grupoNoLineal->addButton(ui->enableSpotIter,4);
    grupoNoLineal->setExclusive(false);
    connect(grupoNoLineal,SIGNAL(buttonClicked(int)),this,SLOT(on_grupoNoLineal_clicked(void)));
}

void opcionesLaser::agrupaBotonesRadios()
{
    QButtonGroup* R11 = new QButtonGroup(this);
    QButtonGroup* R12 = new QButtonGroup(this);
    R11->addButton(ui->R11Inf,1);
    R11->addButton(ui->R11Fin,2);
    R11->setExclusive(true);
    R12->addButton(ui->R12Inf,1);
    R12->addButton(ui->R12Fin,2);
    R12->setExclusive(true);
    connect(R11,SIGNAL(buttonClicked(int)),this,SLOT(on_R11_clicked(void)));
    connect(R12,SIGNAL(buttonClicked(int)),this,SLOT(on_R12_clicked(void)));
}

void opcionesLaser::on_R11_clicked(void)
{
  if(ui->R11Fin->isChecked())
      ui->R11SpinBox->setEnabled(true);
  else
  {
      ui->R11SpinBox->setEnabled(false);
      ui->R11SpinBox->setValue(INFINITY);
  }
}

void opcionesLaser::on_R12_clicked(void)
{
  if(ui->R12Fin->isChecked())
      ui->R12SpinBox->setEnabled(true);
  else
  {
      ui->R12SpinBox->setEnabled(false);
      ui->R12SpinBox->setValue(INFINITY);
  }
}


void opcionesLaser::on_grupoNoLineal_clicked(void)
{
    if(ui->enableEM1ML->isChecked())
    {
        ui->menuRutaEM1ML->setEnabled(true);
        ui->rutaEM1ML->setEnabled(true);
    }
    else
    {
        ui->menuRutaEM1ML->setEnabled(false);
        ui->rutaEM1ML->setEnabled(false);
    }

    if(ui->enableEM2ML->isChecked())
    {
        ui->menuRutaEM2ML->setEnabled(true);
        ui->rutaEM2ML->setEnabled(true);
    }
    else
    {
        ui->menuRutaEM2ML->setEnabled(false);
        ui->rutaEM2ML->setEnabled(false);
    }

    if(ui->enableEM1ML->isChecked()||ui->enableEM2ML->isChecked())
    {
        ui->enableSpotIter->setEnabled(true);
        ui->enableTermico->setEnabled(true);
        ui->iterSpinBox->setEnabled(true);
        ui->numLaminasSpinBox->setEnabled(true);
        ui->umbralSpinBox->setEnabled(true);
        ui->SpotIterMax->setEnabled(true);
        ui->PLaserSpinBox->setEnabled(true);
    }
    else
    {
        ui->enableSpotIter->setEnabled(false);
        ui->enableSpotIter->setChecked(false);
        ui->iterSpinBox->setEnabled(false);
        ui->numLaminasSpinBox->setEnabled(false);
        ui->umbralSpinBox->setEnabled(false);
        ui->enableTermico->setEnabled(false);
        ui->enableTermico->setChecked(false);
        ui->SpotIterMax->setEnabled(false);
        ui->SpotIterMax->setChecked(false);
        ui->PLaserSpinBox->setEnabled(false);
    }
}

void opcionesLaser::on_menuRutaEM1CW_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::AnyFile);
    menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaEM1CW->setText(QString(nombreArchivo[0]));
}

void opcionesLaser::on_menuRutaEM2CW_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::AnyFile);
    menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaEM2CW->setText(QString(nombreArchivo[0]));
}

void opcionesLaser::on_menuRutaEM1ML_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::AnyFile);
    menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaEM1ML->setText(QString(nombreArchivo[0]));
}

void opcionesLaser::on_menuRutaEM2ML_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::AnyFile);
    menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaEM2ML->setText(QString(nombreArchivo[0]));
}

void opcionesLaser::on_enableTermico_toggled(bool checked)
{
    if(checked&&(ui->enableEM1ML->isChecked()||ui->enableEM2ML->isChecked()))
    {
        ui->tamX->setEnabled(true);
        ui->tamY->setEnabled(true);
        ui->N->setEnabled(true);
        ui->Pump->setEnabled(true);
        ui->Distancias->setEnabled(true);
        ui->Lente->setEnabled(true);
        ui->Espejo->setEnabled(true);
        ui->Cristal->setEnabled(true);
    }
    else
    {
        ui->tamX->setEnabled(false);
        ui->tamY->setEnabled(false);
        ui->N->setEnabled(false);
        ui->Pump->setEnabled(false);
        ui->Distancias->setEnabled(false);
        ui->Lente->setEnabled(false);
        ui->Espejo->setEnabled(false);
        ui->Cristal->setEnabled(false);
        //ui->enableTermico->setChecked(0);
    }
}

void opcionesLaser::on_buttonBox_accepted()
{
    // Guardar caminos para otras rutinas.
    string em1CW=ui->rutaEM1CW->text().toStdString();
    string em2CW=ui->rutaEM2CW->text().toStdString();
    string em1ML=ui->rutaEM1ML->text().toStdString();
    string em2ML=ui->rutaEM2ML->text().toStdString();
    string iteracionesEM1=ui->rutaSpotIterEM1->text().toStdString();
    string iteracionesEM2=ui->rutaSpotIterEM2->text().toStdString();
    string iterMax=ui->rutaIterMax->text().toStdString();

    A->setRutas(em1CW,em2CW,em1ML,em2ML,iteracionesEM1,iteracionesEM2,iterMax,parsePtr);
    A->setCasos(ui->enableEM1ML->isChecked(),ui->enableEM2ML->isChecked(),ui->enableTermico->isChecked(),ui->enableSpotIter->isChecked(),\
                ui->SpotIterMax->isChecked(),parsePtr);
    A->setVariables(ui->numLaminasSpinBox->value(),ui->iterSpinBox->value(),ui->umbralSpinBox->value(),ui->N->value(),ui->tamX->value(),ui->tamY->value(),ui->alphaSpinBox->value(),ui->wPumpTan->value(),\
                    ui->wPumpSag->value(),ui->pumpDiv->value(),ui->lambdaPump->value(),ui->pumpP->value(),ui->R11SpinBox->value(),ui->R12SpinBox->value(),\
                    ui->tLSpinBox->value(),ui->nLSpinBox->value(),ui->inclinacionSpinBox->value(),ui->tESpinBox->value(),ui->nESpinBox->value(),ui->FuenteALenteSpinBox->value(),\
                    ui->LenteAEspejoSpinBox->value(),ui->PLaserSpinBox->value(),parsePtr);
    parsePtr->escribeDatos(parsePtr->mapaChar,parsePtr->mapaBool,parsePtr->mapaEnteros,parsePtr->mapaDobles);
}

void opcionesLaser::configuraControles()
{
    ui->numLaminasSpinBox->setMinimum(100);
    ui->numLaminasSpinBox->setMaximum(10000);
    ui->iterSpinBox->setMinimum(1);
    ui->iterSpinBox->setMaximum(100000);
    ui->umbralSpinBox->setMinimum(1e-3);
    ui->umbralSpinBox->setMaximum(100.000);
    ui->umbralSpinBox->setDecimals(3); // Termina de configurar los spinboxes
    ui->PLaserSpinBox->setMaximum(INFINITY);
    ui->PLaserSpinBox->setMinimum(0);
    ui->N->setMinimum(1);
    ui->N->setMaximum(99);
    ui->tamX->setMinimum(1e-3);
    ui->tamX->setMaximum(100);
    ui->tamX->setDecimals(3);
    ui->tamY->setMinimum(1e-3);
    ui->tamY->setMaximum(100);
    ui->tamY->setDecimals(3);
    ui->alphaSpinBox->setMinimum(0);
    ui->alphaSpinBox->setMaximum(1000);
    ui->alphaSpinBox->setDecimals(2);
    ui->FuenteALenteSpinBox->setMinimum(0);
    ui->FuenteALenteSpinBox->setMaximum(INFINITY);
    ui->FuenteALenteSpinBox->setDecimals(3);
    ui->LenteAEspejoSpinBox->setMinimum(0);
    ui->LenteAEspejoSpinBox->setMaximum(INFINITY);
    ui->LenteAEspejoSpinBox->setDecimals(3);
    ui->wPumpTan->setMinimum(0.001);
    ui->wPumpTan->setMaximum(10);
    ui->wPumpTan->setDecimals(3);
    ui->wPumpSag->setMinimum(0.001);
    ui->wPumpSag->setMaximum(10);
    ui->wPumpSag->setDecimals(3);
    ui->pumpDiv->setMinimum(0);
    ui->pumpDiv->setMaximum(1000);
    ui->pumpDiv->setDecimals(2);
    ui->pumpP->setMinimum(0);
    ui->pumpP->setMaximum(100);
    ui->pumpP->setDecimals(2);
    ui->lambdaPump->setMinimum(400);
    ui->lambdaPump->setMaximum(600);
    ui->lambdaPump->setDecimals(2);
    ui->R11SpinBox->setMinimum(0.1);
    ui->R11SpinBox->setMaximum(INFINITY);
    ui->R12SpinBox->setMinimum(0.1);
    ui->R12SpinBox->setMaximum(INFINITY);
    ui->nLSpinBox->setMinimum(1);
    ui->nLSpinBox->setMaximum(10);
    ui->nLSpinBox->setDecimals(5);
    ui->tLSpinBox->setMinimum(0);
    ui->tLSpinBox->setMaximum(50);
    ui->tLSpinBox->setDecimals(3);
    ui->inclinacionSpinBox->setMinimum(-180);
    ui->inclinacionSpinBox->setMaximum(180);
    ui->inclinacionSpinBox->setDecimals(2);
    ui->nESpinBox->setMinimum(1);
    ui->nESpinBox->setMaximum(10);
    ui->nESpinBox->setDecimals(5);
    ui->tESpinBox->setMinimum(0);
    ui->tESpinBox->setMaximum(50);
    ui->tESpinBox->setDecimals(3);

    // Defaults

    if(parsePtr->archivoAbierto())
    {
        parsePtr->cargaVariables();
        ui->numLaminasSpinBox->setValue(parsePtr->mapaEnteros["numeroLaminas"]);
        ui->iterSpinBox->setValue(parsePtr->mapaEnteros["iteraciones"]);
        ui->umbralSpinBox->setValue(parsePtr->mapaDobles["umbral"]);
        ui->PLaserSpinBox->setValue(parsePtr->mapaDobles["PLaser"]*1e-3);
        ui->N->setValue(parsePtr->mapaEnteros["NTerm"]);
        ui->tamX->setValue(parsePtr->mapaDobles["ladoX"]*1e3);
        ui->tamY->setValue(parsePtr->mapaDobles["ladoY"]*1e3);
        ui->alphaSpinBox->setValue(parsePtr->mapaDobles["alpha"]);
        ui->FuenteALenteSpinBox->setValue(parsePtr->mapaDobles["separacionFuenteALente"]);
        ui->LenteAEspejoSpinBox->setValue(parsePtr->mapaDobles["separacionLenteAEspejo"]);
        ui->wPumpTan->setValue(parsePtr->mapaDobles["wPumpTan"]*1e3);
        ui->wPumpSag->setValue(parsePtr->mapaDobles["wPumpSag"]*1e3);
        ui->pumpDiv->setValue(parsePtr->mapaDobles["divPump"]*1e3);
        ui->pumpP->setValue(parsePtr->mapaDobles["PPump"]);
        ui->lambdaPump->setValue(parsePtr->mapaDobles["lambdaPump"]*1e9);
        ui->R11SpinBox->setValue(parsePtr->mapaDobles["RA"]*1e3);
        ui->R12SpinBox->setValue(parsePtr->mapaDobles["RB"]*1e3);
        ui->tLSpinBox->setValue(parsePtr->mapaDobles["tLente"]*1e3);
        ui->nLSpinBox->setValue(parsePtr->mapaDobles["nLente"]);
        ui->tESpinBox->setValue(parsePtr->mapaDobles["tEspejo"]*1e3);
        ui->nESpinBox->setValue(parsePtr->mapaDobles["nEspejo"]);

        // Texto

        ui->rutaEM1CW->setText(parsePtr->mapaChar["RutaEM1CW"]);
        ui->rutaEM2CW->setText(parsePtr->mapaChar["RutaEM2CW"]);
        ui->rutaEM1ML->setText(parsePtr->mapaChar["RutaEM1ML"]);
        ui->rutaEM2ML->setText(parsePtr->mapaChar["RutaEM2ML"]);
        ui->rutaSpotIterEM1->setText(parsePtr->mapaChar["RutaIteracionesEM1"]);
        ui->rutaSpotIterEM2->setText(parsePtr->mapaChar["RutaIteracionesEM2"]);
        ui->rutaIterMax->setText(parsePtr->mapaChar["RutaIterMax"]);
    }
    else
    {
        ui->numLaminasSpinBox->setValue(1000);
        ui->iterSpinBox->setValue(20000);
        ui->umbralSpinBox->setValue(0.1);
        ui->PLaserSpinBox->setValue(100);
        ui->N->setValue(51);
        ui->tamX->setValue(3);
        ui->tamY->setValue(3);
        ui->alphaSpinBox->setValue(100);
        ui->FuenteALenteSpinBox->setValue(1000);
        ui->LenteAEspejoSpinBox->setValue(45);
        ui->wPumpTan->setValue(1);
        ui->wPumpSag->setValue(1);
        ui->pumpDiv->setValue(0.3);
        ui->pumpP->setValue(8);
        ui->lambdaPump->setValue(532);
        ui->R11SpinBox->setValue(INFINITY);
        ui->R12SpinBox->setValue(100);
        ui->tLSpinBox->setValue(2.5);
        ui->nLSpinBox->setValue(1.5195);
        ui->tESpinBox->setValue(6.35);
        ui->nESpinBox->setValue(1.5195);
        ui->rutaEM1CW->setText("./EM1CW.txt");
        ui->rutaEM2CW->setText("./EM2CW.txt");
        ui->rutaEM1ML->setText("./EM1ML.txt");
        ui->rutaEM2ML->setText("./EM2ML.txt");
        //ui->rutaSpotIterEM1->setText(".");
        //ui->rutaSpotIterEM2->setText();
        //ui->rutaIterMax->setText();
    }
}

void opcionesLaser::on_N_editingFinished()
{
    QMessageBox aviso;
    aviso.setWindowTitle("Aviso");
    aviso.addButton(QMessageBox::Ok);
    if(ui->N->value()%2==0)
    {
        aviso.setText("El número de puntos de cálculo del plano térmico debe ser impar.");
        aviso.exec();
    }
}

void opcionesLaser::on_SpotIterMax_toggled(bool checked)
{
    if(checked)
    {
        ui->menuRutaIterMax->setEnabled(true);
        ui->rutaIterMax->setEnabled(true);
    }
    else
    {
        ui->menuRutaIterMax->setEnabled(false);
        ui->rutaIterMax->setEnabled(false);
    }
}

void opcionesLaser::on_menuRutaIterMax_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::AnyFile);
    menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaIterMax->setText(QString(nombreArchivo[0]));
}

void opcionesLaser::on_enableSpotIter_toggled(bool checked)
{
    if(checked)
    {
        ui->menuRutaSpotIterEM1->setEnabled(true);
        ui->rutaSpotIterEM1->setEnabled(true);
        ui->menuRutaSpotIterEM2->setEnabled(true);
        ui->rutaSpotIterEM2->setEnabled(true);
    }
    else
    {
        ui->menuRutaSpotIterEM1->setEnabled(false);
        ui->rutaSpotIterEM1->setEnabled(false);
        ui->menuRutaSpotIterEM2->setEnabled(false);
        ui->rutaSpotIterEM2->setEnabled(false);
    }
}

void opcionesLaser::on_menuRutaSpotIterEM1_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::Directory);
    //menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaSpotIterEM1->setText(QString(nombreArchivo[0]));
}

void opcionesLaser::on_menuRutaSpotIterEM2_clicked()
{
    QFileDialog menu(this);
    menu.setFileMode(QFileDialog::Directory);
    //menu.setNameFilter("All files (*.*)");
    menu.setViewMode(QFileDialog::Detail);
    QStringList nombreArchivo;
    if(menu.exec())
        nombreArchivo=menu.selectedFiles();
    ui->rutaSpotIterEM2->setText(QString(nombreArchivo[0]));
}
