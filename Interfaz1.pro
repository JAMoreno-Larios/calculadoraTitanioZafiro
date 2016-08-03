#-------------------------------------------------
#
# Project created by QtCreator 2016-06-06T21:08:38
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Interfaz1
TEMPLATE = app


SOURCES += main.cpp\
        interfaz1.cpp \
    C/error_iteraciones.c \
    C/lineal.c \
    C/matrices.c \
    C/no_lineal.c \
    C/no_linealMatAstTerm.c \
    opcionesLaser.cpp \
    C/acoplamiento.c \
    C/termico.c \
    ctocppqprogressbar.cpp \
    parser.cpp \
    worker.cpp \
    salidatexto.cpp

HEADERS  += interfaz1.h \
    C/error_iteraciones.h \
    C/lineal.h \
    C/matrices.h \
    C/no_lineal.h \
    C/no_linealMatAstTerm.h \
    opcionesLaser.h \
    C/acoplamiento.h \
    C/termico.h \
    ctocppqprogressbar.h \
    parser.h \
    worker.h \
    salidatexto.h

FORMS    += \
    interfaz1.ui \
    opcionesLaser.ui \
    salidatexto.ui

RESOURCES += \
    recursos.qrc
#CONFIG += c++11
