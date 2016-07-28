#ifdef __cplusplus
extern "C" {
#endif
#ifndef ERROR_ITERACIONES_H_INCLUDED
#define ERROR_ITERACIONES_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void guardaVariaciones(long double * spots, int * iteraciones, int iterMax, long double epsilon1, long double epsilon2, char * tipo, char * subCarpeta,  _Bool termino);
int existeArchivo(const char *nombre);
void guardaSpotIteraciones(long double spotTan, long double spotSag, int iterMax, long double epsilon1, long double epsilon2, char * ruta, _Bool termino);


#endif // ERROR_ITERACIONES_H_INCLUDED

#ifdef __cplusplus
}
#endif
