#include "ctocppqprogressbar.h"

CtoCppQProgressBar::CtoCppQProgressBar(QWidget *parent)
{

}

// access functions

EXPORT_C CtoCppQProgressBar* CtoCppQProgressBar_new()
{
    return new CtoCppQProgressBar();
}

EXPORT_C void CtoCppQProgressBar_delete(CtoCppQProgressBar *bar)
{
    delete bar;
}

EXPORT_C void CtoCppQProgressBar_setMinimum(CtoCppQProgressBar *bar, int minimum)
{
    return bar->setMinimum(minimum);
}

EXPORT_C void CtoCppQProgressBar_setMaximum(CtoCppQProgressBar *bar, int maximum)
{
    return bar->setMaximum(maximum);
}

EXPORT_C void CtoCppQProgressBar_setValue(CtoCppQProgressBar *bar, int value)
{
    return bar->setValue(value);
}

EXPORT_C void CtoCppQProgressBar_repaint(CtoCppQProgressBar *bar)
{
    return bar->parentWidget()->parentWidget()->repaint();
}
