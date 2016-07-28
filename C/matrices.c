/* Librería de operaciones matriciales básicas para números complejos */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
//#include <tgmath.h>
#include <assert.h>
#include <string.h>
#include <math.h>

typedef struct
{
    int filas;
    int columnas;
    long double _Complex * datos;
}matriz;

/* Crea una matriz de filas por columnas con valores igual a 0.
Regresa NULL si filas y columnas <=0 o un puntero a la nueva matriz
en caso contrario */

/*int nuevaMatriz(int filas, int columnas, matriz * m) // REVISAR
{
    if (filas<=0||columnas<=0) return -2;
    // Reserva para la estructura de matriz.
    m=calloc(1,sizeof(matriz));
    //Fijar dimensiones
    m->filas=filas;
    m->columnas=columnas;
    // Reserva para un arreglo de longitud filas * columnas
    m->datos=(double _Complex*)calloc(filas*columnas,sizeof(double _Complex));
    return 0;
}*/


matriz * nuevaMatriz(int filas, int columnas)
{
    if (filas<=0||columnas<=0) return NULL;
    // Reserva para la estructura de matriz.
    matriz *m=(matriz *)calloc(1,sizeof(matriz));
    //Fijar dimensiones
    m->filas=filas;
    m->columnas=columnas;
    // Reserva para un arreglo de longitud filas * columnas
    m->datos=(long double _Complex*)calloc(filas*columnas,sizeof(long double _Complex));
    return m;
}

/* Borra matriz, no elimina variable. */

matriz * borraMatriz(matriz * mtz)
{
    assert(mtz);
    // Libera datos de la matriz
    assert (mtz->datos);
    free(mtz->datos);
    mtz->datos=NULL;
    // Libera a mtz
    free(mtz);
    mtz=NULL;
    return mtz;
}

#define elem(mtz,fila,columna) \
    mtz->datos[(columna-1)*mtz->filas+(fila-1)]

/* Copia una matriz */

matriz * copiaMatriz(matriz * mtz)
{
    if (!mtz) return NULL;
    // Crea nueva matriz para guardar copya
    matriz * cp= nuevaMatriz(mtz->filas,mtz->columnas);
    // Copia datos de la matriz a cp.
    memcpy(cp->datos, mtz->datos, mtz->filas*mtz->columnas*sizeof(long double _Complex));
    return cp;
}

/* Imprime la matriz a stdout */

int fijaElemento(matriz * mtz, int fila, int columna, long double _Complex valor)
{
    if (!mtz) return -1;
    assert(mtz->datos);
    if (fila<=0 || fila > mtz->filas || columna<=0 || columna > mtz->columnas)
        return -2;

    elem(mtz,fila,columna)=valor;
    return 0;
}

/* Recupera elemento de la matriz */

long double _Complex obtieneElemento(matriz * mtz, int fila, int columna)
{
    long double _Complex valor;
    if(!mtz) return -1;
    assert(mtz->datos);
    if(fila <=0 || fila > mtz->filas || columna<=0 ||columna>mtz->columnas) return -2;
    valor=elem(mtz,fila,columna);
    return  valor;
}

/* Obtiene el número de filas de la matriz */

int nFilas(matriz * mtz)
{
    int *n;
    if (!mtz) return -1;
    *n = mtz->filas;
    return *n;
    }
/* Obtiene el número de columnas de la matriz */
int nColumnas(matriz * mtz)
{
    int *n;
    if (!mtz) return -1;
    *n=mtz->columnas;
    return *n;
}
int imprimeMatriz(matriz * mtz)
{
    if (!mtz) return -1;

    int fila, columna;
    for (fila=1;fila<=mtz->filas;++fila)
    {
        for(columna=1;columna<=mtz->columnas;++columna)
        {
            //Imprime el elemento de punto flotante
            printf("%.15Le+%.15Lei ", creall(elem(mtz,fila,columna)),cimagl(elem(mtz,fila,columna)));
        }
        printf("\n");
    }
    return 0;
}

matriz * sumaMatriz(matriz * A, matriz * B)
{
    assert(A&&B);
    assert(A->filas==B->filas);
    assert(A->columnas==B->columnas);
    matriz * suma;
    suma=nuevaMatriz(A->filas,A->columnas); // A == B
    int fila, columna;
    for (columna=A->columnas;columna>=1;columna--)
        for (fila=A->filas;fila>=1;fila--)
            elem(suma,fila,columna)=elem(A,fila,columna)+elem(B,fila,columna);
    return suma;
}

matriz * multMatriz(matriz * A, matriz * B)
{
    assert(A&&B);
    assert(A->columnas==B->filas);
    int fila, columna, k;
    matriz * prod=nuevaMatriz(A->filas,B->columnas);
    for (columna=B->columnas;columna>=1;columna--)
        for (fila=A->filas;fila>=1;fila--)
            {
                long double _Complex valor=0+I*0;
                for (k=A->columnas;k>=1;k--)
                    valor+=elem(A,fila,k)*elem(B,k,columna);
                elem(prod,fila,columna)=valor;
            }
    return prod;
}

matriz * cteMultMatriz(long double _Complex Constante, matriz * mtz)
{
    assert(mtz);
    matriz *prod;
    prod=nuevaMatriz(mtz->filas,mtz->columnas);
    prod=copiaMatriz(mtz);
    int fila, columna;
    for(columna=1;columna<=mtz->columnas;columna++)
        for(fila=1;fila<=mtz->filas;fila++)
            elem(prod,fila,columna)=Constante*elem(prod,fila,columna);
    return prod;
}

matriz * divMatricesElemAElem(matriz * A, matriz * B) // División elemento a elemento de A/B: A11/B11, A12/B12, etc. Si número es no finito, lo hace cero.
{
    assert(A&&B);
    assert(A->filas==B->filas);
    assert(A->columnas==B->columnas);
    matriz * division;
    division=nuevaMatriz(A->filas,A->columnas); // A == B
    int fila, columna;
    for (columna=A->columnas;columna>=1;columna--)
        for (fila=A->filas;fila>=1;fila--){
           elem(division,fila,columna)=elem(A,fila,columna)/elem(B,fila,columna);
           long double divResultadoAbs=cabsl(elem(division,fila,columna));
           //if(isfinite(elem(division,fila,columna))==0)
           if(isfinite(divResultadoAbs)==0)
            elem(division,fila,columna)=0;
        }
    return division;
}

void escribeMatrizArchivo(matriz * mtz, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int fila, columna;
    for (fila=1;fila<=mtz->filas;fila++)
    {
        for(columna=1;columna<=mtz->columnas;columna++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le+%.15Lei,", creall(elem(mtz,fila,columna)),cimagl(elem(mtz,fila,columna)));
            //fprintf(archivo, "%6.5f\t", cimagl(elem(mtz,fila,columna)));
        }
        fprintf(archivo,"\n");
    }
    fclose(archivo);
}

void escribeMatrizArchivo_ejes(matriz * mtz, matriz * epsilon1, matriz * epsilon2, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int fila, columna;
    fprintf(archivo, " ,epsilon 2 [m],");
    for(columna=1;columna<=mtz->columnas;columna++)
        fprintf(archivo, "%.15Le,",creall(elem(epsilon2,1,columna)));
    fprintf(archivo, "\nepsilon 1 [m],\n");
    for (fila=1;fila<=mtz->filas;fila++)
    {
        fprintf(archivo,"%.15Le, ,",creall(elem(epsilon1,1,fila)));
        for(columna=1;columna<=mtz->columnas;columna++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le+%.15Lei,", creall(elem(mtz,fila,columna)),cimagl(elem(mtz,fila,columna)));
            //fprintf(archivo, "%6.5f\t", cimagl(elem(mtz,fila,columna)));
        }
        fprintf(archivo,"\n");
    }
    fclose(archivo);
}

void escribeMatrizArchivo_completo(matriz * mtz, matriz * epsilon1, matriz * epsilon2, long double *angulos, long double delta1, long double delta2, char* conjugado_corto, char *conjugado_largo, char *nombre, int fila_noAstig, int columna_noAstig, long double astigmatismo)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int fila, columna;
    fprintf(archivo, "Configuracion optima,\n\n");
    fprintf(archivo,"delta 1 [m],Theta 1 [rad],delta 2 [m],Theta 2 [rad],Conjugado corto, Conjugado largo, Astigmatismo (spot tan / spot sag)\n");
    fprintf(archivo,"%.15Le,%.15Le,%.15Le,%.15Le,%s,%s,%.15Le\n\n",delta1+creall(obtieneElemento(epsilon1,1,fila_noAstig)),angulos[0],delta2+creall(obtieneElemento(epsilon2,1,columna_noAstig)),angulos[1],conjugado_corto,conjugado_largo, creall(astigmatismo));
    fprintf(archivo, "epsilon 1 [m],epsilon2 [m]\n");
    fprintf(archivo, "%.15Le,%.15Le,\n\n",creall(obtieneElemento(epsilon1,1,fila_noAstig)),creall(obtieneElemento(epsilon2,1,columna_noAstig)));
    fprintf(archivo, " ,epsilon 2 [m],");
    for(columna=1;columna<=mtz->columnas;columna++)
        fprintf(archivo, "%.15Le,",creall(elem(epsilon2,1,columna)));
    fprintf(archivo, "\nepsilon 1 [m],\n");
    for (fila=1;fila<=mtz->filas;fila++)
    {
        fprintf(archivo,"%.15Le, ,",creall(elem(epsilon1,1,fila)));
        for(columna=1;columna<=mtz->columnas;columna++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le+%.15Lei,", creall(elem(mtz,fila,columna)),cimagl(elem(mtz,fila,columna)));
            //fprintf(archivo, "%6.5f\t", cimagl(elem(mtz,fila,columna)));
        }
        fprintf(archivo,"\n");
    }
    fclose(archivo);
}

void escribeMatrizDoubleArchivo(matriz * mtz, char *nombre)
{
    FILE *archivo;
    archivo = fopen(nombre,"w");
    int fila, columna;
    fprintf(archivo,"double datos={");
    for (fila=1;fila<=mtz->filas-1;fila++)
    {
        fprintf(archivo,"{");
        for(columna=1;columna<=mtz->columnas;columna++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le,", creall(elem(mtz,fila,columna)));
            //fprintf(archivo, "%6.5f\t", cimagl(elem(mtz,fila,columna)));
        }
        fprintf(archivo,"},\n");
    }
    fila=mtz->filas;
    fprintf(archivo,"{");
    for(columna=1;columna<=mtz->columnas;columna++)
        {
            //Imprime el elemento de punto flotante
            fprintf(archivo, "%.15Le,", creall(elem(mtz,fila,columna)));
            //fprintf(archivo, "%6.5f\t", cimagl(elem(mtz,fila,columna)));
        }
    fprintf(archivo,"}};");
    fclose(archivo);
}


// Funciones válidas para matrices de 2x2

matriz * llenaMatriz(long double _Complex A, long double _Complex B,long double _Complex C, long double _Complex D)
{
    matriz * mtz;
    mtz = nuevaMatriz(2,2);
    fijaElemento(mtz,1,1,A);
    fijaElemento(mtz,1,2,B);
    fijaElemento(mtz,2,1,C);
    fijaElemento(mtz,2,2,D);
    //borraMatriz(mtz);
    return mtz;
}


void llenaMatrizSinCrear(matriz * mtz, long double _Complex A, long double _Complex B,long double _Complex C, long double _Complex D)
{
    fijaElemento(mtz,1,1,A);
    fijaElemento(mtz,1,2,B);
    fijaElemento(mtz,2,1,C);
    fijaElemento(mtz,2,2,D);
}

matriz * variasMultMatriciales(matriz ** arreglo_matrices, int num_matrices) // A1*A2*A3*...*An
{
    matriz * temp=llenaMatriz(1,0,0,1);
    int contador=num_matrices;
    for(contador=num_matrices-1;contador>=0;contador--)
    {
        matriz * hold=multMatriz(arreglo_matrices[contador],temp);
        llenaMatrizSinCrear(temp,obtieneElemento(hold,1,1),obtieneElemento(hold,1,2),obtieneElemento(hold,2,1),obtieneElemento(hold,2,2));
        borraMatriz(hold);
    }

    return temp;
}

// Función para imprimir datos de cálculos para gnuplot

void imprimeListaReal(matriz * ejeX, matriz * ejeY, matriz * datosTan, matriz * datosSag, char* rutaArchivo)
{
    FILE *archivo;
    archivo = fopen(rutaArchivo,"w");
    fprintf(archivo,"#Datos a graficar con Gnuplot, e1, e2, tan, sag\n");
    for(int i=1;i<=ejeX->columnas;i++)
    {
        for(int j=1;j<=ejeY->columnas;j++)
        {
            fprintf(archivo,"%Le %Le %Le %Le\n",creall(obtieneElemento(ejeX,1,i)),creall(obtieneElemento(ejeY,1,j)),creall(obtieneElemento(datosTan,i,j)),creall(obtieneElemento(datosSag,i,j)));
        }
    }
    fclose(archivo);
}

void imprimeListaCompleja(matriz * ejeX, matriz * ejeY, matriz * datosTan, matriz * datosSag, char* rutaArchivo)
{
    FILE *archivo;
    archivo = fopen(rutaArchivo,"w");
    fprintf(archivo,"#Datos a graficar con Gnuplot, e1, e2, tan, sag\n");
    for(int i=1;i<=ejeX->columnas;i++)
    {
        for(int j=1;j<=ejeY->columnas;j++)
        {
            fprintf(archivo,"%Le %Le %Le+%Lei %Le+%Lei\n",creall(obtieneElemento(ejeX,1,i)),creall(obtieneElemento(ejeY,1,j)),creall(obtieneElemento(datosTan,i,j)),cimagl(obtieneElemento(datosTan,i,j)),creall(obtieneElemento(datosSag,i,j)),cimagl(obtieneElemento(datosSag,i,j)));
        }
    }
    fclose(archivo);
}
