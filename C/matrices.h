#ifdef __cplusplus
extern "C" {
#endif

#ifndef MATRICES_H_INCLUDED
#define MATRICES_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
//#include <complex.h>
#include <tgmath.h>
#include <assert.h>
#include <string.h>

typedef struct
{
    int filas;
    int columnas;
    long double _Complex * datos;
}matriz;

#define elem(mtz,fila,columna) mtz->datos[(columna-1)*mtz->filas+(fila-1)]

#define length(x) (sizeof(x)/sizeof((x)[0]))

/* Crea una matriz de filas por columnas con valores igual a 0.
Regresa NULL si filas y columnas <=0 o un puntero a la nueva matriz
en caso contrario */

matriz * nuevaMatriz(int filas, int columnas);
/* Borra matriz, no elimina variable. */
matriz * borraMatriz(matriz * mtz);
/* Copia una matriz */
matriz * copiaMatriz(matriz * mtz);
/* Imprime la matriz a stdout */
int fijaElemento(matriz * mtz, int fila, int columna, long double _Complex valor);
/* Recupera elemento de la matriz */
long double _Complex obtieneElemento(matriz * mtz, int fila, int columna);
/* Obtiene el número de filas de la matriz */
int nFilas(matriz * mtz);
/* Obtiene el número de columnas de la matriz */
int nColumnas(matriz * mtz);
/* Imprime matriz a stdout */
int imprimeMatriz(matriz * mtz);
matriz * sumaMatriz(matriz * A, matriz * B);
matriz * multMatriz(matriz * A, matriz * B);
matriz * cteMultMatriz(long double _Complex Constante, matriz * mtz);
matriz * divMatricesElemAElem(matriz * A, matriz * B);
void escribeMatrizArchivo(matriz * mtz, char *nombre);
void escribeMatrizArchivo_ejes(matriz * mtz, matriz * epsilon1, matriz * epsilon2, char *nombre);
void escribeMatrizArchivo_completo(matriz * mtz, matriz * epsilon1, matriz * epsilon2, long double *angulos, long double delta1, long double delta2, char* conjugado_corto, char *conjugado_largo, char *nombre, int fila_noAstig, int columna_noAstig, long double astigmatismo);
void escribeMatrizDoubleArchivo(matriz * mtz, char *nombre);
// Sólo para matrices de 2x2
matriz * llenaMatriz(long double _Complex A, long double _Complex B,long double _Complex C, long double _Complex D);
void llenaMatrizSinCrear(matriz * mtz, long double _Complex A, long double _Complex B,long double _Complex C, long double _Complex D);
matriz * variasMultMatriciales(matriz ** arreglo_matrices, int num_matrices);
void imprimeListaReal(matriz * ejeX, matriz * ejeY, matriz * datosTan, matriz * datosSag, char* rutaArchivo);
void imprimeListaCompleja(matriz * ejeX, matriz * ejeY, matriz * datosTan, matriz * datosSag, char* rutaArchivo);

#endif // MATRICES_H_INCLUDED

#ifdef __cplusplus
}
#endif
