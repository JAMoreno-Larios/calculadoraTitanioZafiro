#ifndef SALIDATEXTO_H
#define SALIDATEXTO_H

#include <QMainWindow>
#include <ctime>

namespace Ui {
class salidaTexto;
}

class salidaTexto : public QMainWindow
{
    Q_OBJECT

public:
    explicit salidaTexto(QWidget *parent = 0);
    ~salidaTexto();
    void escribeLineaTexto(QString texto);

private slots:
    void on_actionGuardar_triggered();

    void on_actionCerrar_triggered();

private:
    Ui::salidaTexto *ui;
    char fecha[12];
    time_t rawtime;
    struct tm * timeinfo;

};

#endif // SALIDATEXTO_H
