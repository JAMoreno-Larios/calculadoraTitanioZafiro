#ifndef CTOCPPQPROGRESSBAR_H
#define CTOCPPQPROGRESSBAR_H

#ifdef __cplusplus
#include <QWidget>
#include <QProgressBar>
class CtoCppQProgressBar : public QProgressBar
{
public:
    explicit CtoCppQProgressBar(QWidget *parent=0);
};

#else

typedef struct CtoCppQProgressBar CtoCppQProgressBar;

#endif

// Funciones de entrada.
// Macros
#ifdef __cplusplus
    #define EXPORT_C extern "C"
#else
    #define EXPORT_C
#endif

EXPORT_C CtoCppQProgressBar* CtoCppQProgressBar_new(void);
EXPORT_C void CtoCppQProgressBar_delete(CtoCppQProgressBar *bar);
EXPORT_C void CtoCppQProgressBar_setMinimum(CtoCppQProgressBar *bar, int minimum);
EXPORT_C void CtoCppQProgressBar_setMaximum(CtoCppQProgressBar *bar, int maximum);
EXPORT_C void CtoCppQProgressBar_setValue(CtoCppQProgressBar *bar, int value);
EXPORT_C void CtoCppQProgressBar_repaint(CtoCppQProgressBar *bar);

#endif // CTOCPPQPROGRESSBAR_H
