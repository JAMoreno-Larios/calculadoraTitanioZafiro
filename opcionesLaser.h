#ifndef OPCIONESLASER_H
#define OPCIONESLASER_H

#include <QDialog>
#include "interfaz1.h"
#include "parser.h"
using std::string;
class interfaz1;

namespace Ui {
class opcionesLaser;
}

class opcionesLaser : public QDialog
{
    Q_OBJECT

public:
    explicit opcionesLaser(parser *Par, QWidget *parent = 0);
    ~opcionesLaser();
    interfaz1 *A;

private slots:
    void agrupaBotonesNoLineal();

    void agrupaBotonesRadios();

    void on_R11_clicked(void);

    void on_R12_clicked(void);

    void on_grupoNoLineal_clicked(void);

    void on_menuRutaEM1CW_clicked();

    void on_menuRutaEM2CW_clicked();

    void on_menuRutaEM1ML_clicked();

    void on_menuRutaEM2ML_clicked();

    void on_enableSpotIter_toggled(bool checked);

    void on_menuRutaSpotIterEM1_clicked();

    void on_buttonBox_accepted();

    void on_enableTermico_toggled(bool checked);

    void configuraControles();

    void on_N_editingFinished();

    void on_SpotIterMax_toggled(bool checked);

    void on_menuRutaIterMax_clicked();

    void on_menuRutaSpotIterEM2_clicked();

private:
    Ui::opcionesLaser *ui;
    parser *parsePtr;
};

#endif // OPCIONESLASER_H
